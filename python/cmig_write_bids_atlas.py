#!/usr/bin/env python3
"""
Install the ABCD T1 atlas into a BIDS derivative dataset.

cmig_proc_bids_anat writes SpatialReference as bids::atlases/atlas-ABCD_T1w.nii.gz
for scans resampled into the atlas. A bids:: URI resolves within the dataset that
contains it, so the atlas has to live inside the derivative tree rather than being
referenced from the pipeline install. This is a dataset-setup step, run once per root rather than
from a per-session module, since parallel per-session jobs would race.

Conversion uses mri_convert, so no nibabel dependency is introduced; a sourced
FreeSurfer is therefore a prerequisite. The voxel grid must survive the
conversion exactly, since the file's whole purpose is to be the geometric
reference the URI names, so the output grid is verified against the input.
Each file is written to a temporary file beside it and renamed, so concurrent
jobs, or one killed during the conversion, never leave a partial file.

Also writes the root-level atlas-at_dseg.tsv lookup into the analysis root
(write_dseg_lookup), which the AtlasTrack dseg of every session inherits.

Created:  08/10/2026 by Don Hagler
Prev mod: 08/10/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/13/2026 by Don Hagler
"""

import argparse
import collections
import json
import os
import subprocess
import sys
import tempfile

from cmig_sing_utils import cmig_getenv

# None: $CMIG_PIPELINE_ATLASDIR/T1_Atlas, read when used
DEFAULT_ATLAS_DIR = None

SUBDIR = 'atlases'

# source file -> output file; the mask carries no _T1w, since its suffix is mask
ATLAS_FILES = [
    ('T1_atlas_ABCD.mgz',      'atlas-ABCD_T1w.nii.gz'),
    ('T1_atlas_ABCD_mask.mgz', 'atlas-ABCD_desc-brain_mask.nii.gz'),
]

# the DatasetLinks label of the atlas package (cmig_pipeline_atlases_<tag>.tar.gz)
ATLAS_DATASET_LABEL = 'cmig_pipeline_atlases'

FNAME_JSON = 'atlas-ABCD_T1w.json'
DESCRIPTION = ('ABCD T1-weighted atlas, the geometric reference for scans resampled '
               'into atlas space by cmig_proc_bids_anat')
FNAME_MASK_JSON = 'atlas-ABCD_desc-brain_mask.json'
MASK_DESCRIPTION = 'Brain mask of the ABCD T1-weighted atlas'

# the AtlasTrack dseg lookup in the analysis root, from CMIG_ColorLUT.txt: the
#   dseg labels are 21000 plus the fiber number, and 21202-21209 for voxels in
#   two to nine overlapping fibers
FNAME_DSEG_TSV = 'atlas-at_dseg.tsv'
DSEG_CODE_RANGES = [(21100, 21199), (21202, 21209)]


def _mri_info(fname, flag):
    """Return one mri_info field as a normalized string, or None on failure."""
    try:
        out = subprocess.run(['mri_info', flag, fname], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True)
    except OSError:
        return None
    if out.returncode != 0:
        return None
    return ' '.join(out.stdout.split())


def check_grid(fname_in, fname_out):
    """Verify the output voxel grid matches the input.

    Returns (errcode, errmsg). Comparison is on dimensions, voxel sizes and the
    vox2ras matrix; if mri_info cannot report them, that is reported rather than
    passed over, since an unverified grid is the failure this guards against.
    """
    for flag in ('--dim', '--res', '--vox2ras'):
        val_in = _mri_info(fname_in, flag)
        val_out = _mri_info(fname_out, flag)
        if val_in is None or val_out is None:
            return 1, 'could not read {0} via mri_info for grid check'.format(flag)
        if val_in != val_out:
            return 1, 'grid mismatch on {0}: input {1}, output {2}'.format(
                flag, val_in, val_out)
    return 0, None


def make_temp(fname_out, suffix='.tmp'):
    """Create an empty temporary file beside fname_out (tempfile.mkstemp, a
    unique name in the same directory, so os.replace can rename it onto
    fname_out), with a new file's mode under the umask rather than 0600.
    Returns its path; raises OSError."""
    fd, fname_tmp = tempfile.mkstemp(
        prefix='.{0}.'.format(os.path.basename(fname_out)), suffix=suffix,
        dir=os.path.dirname(fname_out) or '.')
    os.close(fd)
    umask = os.umask(0)
    os.umask(umask)
    os.chmod(fname_tmp, 0o666 & ~umask)
    return fname_tmp


def _remove(fname):
    if fname and os.path.isfile(fname):
        os.remove(fname)


def convert_atlas_file(fname_in, fname_out, forceflag=0):
    """Convert one atlas file with mri_convert into a temporary file, verify
    the grid, then rename it to fname_out.

    Returns (errcode, errmsg). An existing output is left alone unless forceflag
    is set, since the setup is idempotent and will be re-run.
    """
    if os.path.isfile(fname_out) and not forceflag:
        return 0, None
    if not os.path.isfile(fname_in):
        return 1, 'atlas source not found: {0}'.format(fname_in)

    # mri_convert picks the output format from the extension
    suffix = '.nii.gz' if fname_out.endswith('.nii.gz') else os.path.splitext(fname_out)[1]
    try:
        fname_tmp = make_temp(fname_out, suffix=suffix)
    except OSError as err:
        return 1, 'cannot create a temporary file for {0}: {1}'.format(fname_out, err)
    try:
        out = subprocess.run(['mri_convert', fname_in, fname_tmp],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)
    except OSError as err:
        _remove(fname_tmp)
        return 1, 'failed to run mri_convert, is FreeSurfer sourced: {0}'.format(err)
    if out.returncode != 0:
        _remove(fname_tmp)
        return 1, 'mri_convert failed for {0}: {1}'.format(
            os.path.basename(fname_in), out.stdout.strip().splitlines()[-1:] or '')
    if not os.path.getsize(fname_tmp):
        _remove(fname_tmp)
        return 1, 'mri_convert wrote no output for {0}'.format(os.path.basename(fname_in))

    errcode, errmsg = check_grid(fname_in, fname_tmp)
    if errcode:
        _remove(fname_tmp)
        return errcode, errmsg
    try:
        os.replace(fname_tmp, fname_out)
    except OSError as err:
        _remove(fname_tmp)
        return 1, 'failed to write {0}: {1}'.format(fname_out, err)
    return 0, None


def _write_file(fname_out, write_fn):
    """Write through a temporary file beside fname_out and rename.
    Returns (errcode, errmsg)."""
    fname_tmp = None
    try:
        fname_tmp = make_temp(fname_out)
        with open(fname_tmp, 'w') as fid:
            write_fn(fid)
        os.replace(fname_tmp, fname_out)
    except (OSError, IOError, TypeError, ValueError) as err:
        _remove(fname_tmp)
        return 1, 'failed to write {0}: {1}'.format(fname_out, err)
    return 0, None


def write_json(fname_json, info):
    """Write a JSON sidecar. Returns (errcode, errmsg)."""
    def write_fn(fid):
        json.dump(info, fid, indent=4)
        fid.write('\n')
    return _write_file(fname_json, write_fn)


def source_info(atlas_dir, fname_src):
    """Sources of an atlas file: a bids: URI in the atlas package (the
    cmig_pipeline_atlases DatasetLinks entry of the proc root, its release
    asset), the atlas folder name (T1_Atlas) and the source file name."""
    subdir = os.path.basename(os.path.normpath(atlas_dir))
    info = collections.OrderedDict()
    info['Sources'] = ['bids:{0}:{1}/{2}'.format(ATLAS_DATASET_LABEL, subdir,
                                                 fname_src)]
    return info


def write_atlas_json(fname_json, atlas_dir, fname_src, forceflag=0):
    """Write the sidecar for the atlas image.

    Returns (errcode, errmsg).
    """
    if os.path.isfile(fname_json) and not forceflag:
        return 0, None

    info = collections.OrderedDict()
    info['Description'] = DESCRIPTION
    info.update(source_info(atlas_dir, fname_src))
    # the atlas defines its own space
    info['SpatialReference'] = 'orig'
    return write_json(fname_json, info)


def write_mask_json(fname_json, atlas_dir, fname_src, forceflag=0):
    """Write the sidecar for the atlas brain mask, on the atlas image grid.

    Returns (errcode, errmsg).
    """
    if os.path.isfile(fname_json) and not forceflag:
        return 0, None

    info = collections.OrderedDict()
    info['Description'] = MASK_DESCRIPTION
    info.update(source_info(atlas_dir, fname_src))
    info['Type'] = 'Brain'
    info['SpatialReference'] = 'bids::{0}/{1}'.format(SUBDIR, ATLAS_FILES[0][1])
    info['SkullStripped'] = False
    return write_json(fname_json, info)


def read_color_lut(fname_lut):
    """Read a FreeSurfer-style color table: {code: (name, r, g, b)}.
    Comments and blank lines are skipped."""
    lut = collections.OrderedDict()
    with open(fname_lut) as fid:
        for line in fid:
            tok = line.split('#', 1)[0].split()
            if len(tok) < 5 or not tok[0].isdigit():
                continue
            lut[int(tok[0])] = (tok[1], int(tok[2]), int(tok[3]), int(tok[4]))
    return lut


def write_dseg_lookup(analysis_rootdir, fname_lut=None, verbose=0, forceflag=0):
    """Write {analysis_rootdir}/atlas-at_dseg.tsv (index, name, color) for the
    AtlasTrack dseg codes in CMIG_ColorLUT.txt (default
    $CMIG_PIPELINE_DIR/parms/CMIG_ColorLUT.txt).

    Returns (errcode, errmsg).
    """
    fname_out = os.path.join(analysis_rootdir, FNAME_DSEG_TSV)
    if os.path.isfile(fname_out) and not forceflag:
        return 0, None
    if fname_lut is None:
        fname_lut = os.path.join(cmig_getenv('CMIG_PIPELINE_DIR'), 'parms',
                                 'CMIG_ColorLUT.txt')
    try:
        lut = read_color_lut(fname_lut)
    except (OSError, IOError, ValueError) as err:
        return 1, 'failed to read {0}: {1}'.format(fname_lut, err)

    rows = []
    for code, (name, r, g, b) in sorted(lut.items()):
        if any(lo <= code <= hi for lo, hi in DSEG_CODE_RANGES):
            rows.append('{0}\t{1}\t#{2:02X}{3:02X}{4:02X}'.format(code, name, r, g, b))
    if not rows:
        return 1, 'no AtlasTrack codes in {0}'.format(fname_lut)

    try:
        os.makedirs(analysis_rootdir, exist_ok=True)
    except OSError as err:
        return 1, 'cannot create {0}: {1}'.format(analysis_rootdir, err)
    errcode, errmsg = _write_file(
        fname_out, lambda fid: fid.write('index\tname\tcolor\n' +
                                         ''.join(r + '\n' for r in rows)))
    if not errcode and verbose:
        print('wrote {0}'.format(fname_out))
    return errcode, errmsg


def write_bids_atlas(proc_rootdir, atlas_dir=DEFAULT_ATLAS_DIR, verbose=0, forceflag=0):
    """Install the atlas into {proc_rootdir}/atlases/.

    Only proc_rootdir gets a copy: cmig_proc_bids_anat is the only module writing
    the atlas URI, and a bids:: URI cannot reach across datasets.

    Returns (errcode, errmsg).
    """
    if atlas_dir is None:
        atlas_dir = os.path.join(cmig_getenv('CMIG_PIPELINE_ATLASDIR'), 'T1_Atlas')
    if not atlas_dir:
        errmsg = 'no atlas directory supplied'
        print('WARNING: {0}'.format(errmsg))
        return 1, errmsg

    outdir = os.path.join(proc_rootdir, SUBDIR)
    try:
        os.makedirs(outdir, exist_ok=True)
    except OSError as err:
        errmsg = 'cannot create {0}: {1}'.format(outdir, err)
        print('WARNING: {0}'.format(errmsg))
        return 1, errmsg

    for fname_src, fname_dst in ATLAS_FILES:
        fname_in = os.path.join(atlas_dir, fname_src)
        fname_out = os.path.join(outdir, fname_dst)
        existed = os.path.isfile(fname_out)
        errcode, errmsg = convert_atlas_file(fname_in, fname_out, forceflag)
        if errcode:
            print('WARNING: {0}'.format(errmsg))
            return errcode, errmsg
        if verbose:
            print('{0} {1}'.format('kept' if existed and not forceflag else 'wrote',
                                   fname_out))

    for fname, write_fn, fname_src in [
            (FNAME_JSON, write_atlas_json, ATLAS_FILES[0][0]),
            (FNAME_MASK_JSON, write_mask_json, ATLAS_FILES[1][0])]:
        fname_json = os.path.join(outdir, fname)
        existed = os.path.isfile(fname_json)
        errcode, errmsg = write_fn(fname_json, atlas_dir, fname_src, forceflag)
        if errcode:
            print('WARNING: {0}'.format(errmsg))
            return errcode, errmsg
        if verbose:
            print('{0} {1}'.format('kept' if existed and not forceflag else 'wrote',
                                   fname_json))

    return 0, None


def parse_args():
    parser = argparse.ArgumentParser(
        description='Install the ABCD T1 atlas into a BIDS derivative dataset, '
                    'and the atlas-at_dseg.tsv lookup into the analysis root.'
    )
    parser.add_argument('--proc_rootdir',
                        help='Root directory of processed BIDS derivatives, '
                             'for the atlas')
    parser.add_argument('--analysis_rootdir',
                        help='Root directory of analysis output, for the '
                             'atlas-at_dseg.tsv lookup')
    parser.add_argument('--atlas_dir', default=DEFAULT_ATLAS_DIR,
                        help='Directory holding the atlas source files '
                             '(default: $CMIG_PIPELINE_ATLASDIR/T1_Atlas)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='Display status messages (default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing output (default: 0)')
    return parser.parse_args()


def run(args):
    errcode, errmsg = 0, None
    if not args.proc_rootdir and not args.analysis_rootdir:
        errmsg = 'neither --proc_rootdir nor --analysis_rootdir supplied'
        print('ERROR: {0}'.format(errmsg))
        return 1, errmsg
    if args.proc_rootdir:
        errcode, errmsg = write_bids_atlas(
            args.proc_rootdir, atlas_dir=args.atlas_dir,
            verbose=args.verbose, forceflag=args.forceflag)
    if args.analysis_rootdir:
        tmp_errcode, tmp_errmsg = write_dseg_lookup(
            args.analysis_rootdir, verbose=args.verbose, forceflag=args.forceflag)
        if tmp_errcode:
            print('WARNING: {0}'.format(tmp_errmsg))
            errcode, errmsg = tmp_errcode, tmp_errmsg
    return errcode, errmsg


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
