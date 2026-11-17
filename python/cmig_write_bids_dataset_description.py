#!/usr/bin/env python3
"""
Write dataset_description.json and .bidsignore for cmig BIDS derivative
datasets.

Each derivative root is a separate BIDS dataset and requires its own
dataset_description.json. This is a dataset-setup step, run once per root
rather than from a per-session module, since parallel per-session jobs would
race to write the same file.

The file also anchors the bids:: URIs written as SpatialReference: the dataset
root is the nearest parent directory containing a valid dataset_description.json,
so without one a URI has no anchor. Its DatasetLinks resolve the bids:<label>:
URIs that name a file in another root (e.g. bids:cmig_proc:sub-X/...).

The .bidsignore beside it lists the root's known non-BIDS content, and the
analysis root also gets the atlas-at_dseg.tsv lookup (cmig_write_bids_atlas).

Created:  08/09/2026 by Don Hagler
Prev mod: 08/09/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/13/2026 by Don Hagler
"""

import argparse
import collections
import json
import os
import sys

import cmig_sing_utils
import cmig_write_bids_atlas

DEFAULT_BIDS_VERSION = '1.11.1'
# the first GeneratedBy Name, which BIDS requires to be a substring of the
#   root's directory name (check_root_name)
DEFAULT_PIPELINE_NAME = 'cmig'
# no longer written: the single cmig entry names the pipeline and carries its
#   version; kept so the callers' --software_name still parses
DEFAULT_SOFTWARE_NAME = 'cmig_pipeline'
# None: $CMIG_PIPELINE_VERSION, read when used
DEFAULT_SOFTWARE_VERSION = None
DEFAULT_FS_VERSION = 820
FS_NAME = 'freesurfer'
CODE_URL = 'https://github.com/cmig-research-group/cmig_pipeline'
# the CodeURL of a tagged release; an untagged build ('local') gets CODE_URL
TAG_CODE_URL = CODE_URL + '/tree/{0}'
# the release tag and the image's name:tag (Dockerfile ENV); versions.json
#   keys cmig_pipeline_tag and image_tag take precedence
PIPELINE_TAG_ENV = 'CMIG_PIPELINE_TAG'
IMAGE_TAG_ENV = cmig_sing_utils.IMAGE_TAG_ENV
AFNI_CODE_URL = 'https://github.com/afni/afni/tree/{0}'
# the atlas package of a tagged atlases version (versions' atlases), the
#   release asset the image fetches; an untagged one ('local') or none (the
#   host) gets no link
ATLAS_ASSET_URL = (CODE_URL + '/releases/download/{0}/'
                   + cmig_write_bids_atlas.ATLAS_DATASET_LABEL + '_{0}.tar.gz')

FNAME_DESC = 'dataset_description.json'
FNAME_IGNORE = '.bidsignore'

# session-root files outside any datatype directory: protocol compliance,
#   pcqc, session info, status, manual raw QC (new and old names), review locks
SESSION_QC_IGNORE = [
    '# session-root QC files, outside any datatype directory',
    'sub-*/ses-*/*_pc.tsv',
    'sub-*/ses-*/*_pcqc.tsv',
    'sub-*/ses-*/*_info.tsv',
    'sub-*/ses-*/*_status.tsv',
    'sub-*/ses-*/*_status.json',
    'sub-*/ses-*/*_manualrawqc.tsv',
    'sub-*/ses-*/*_desc-manualrawqc_*.tsv',
    'sub-*/ses-*/*.lock',
]

# per-datatype status files
STATUS_IGNORE = [
    '# per-datatype status files',
    'sub-*/ses-*/*/*_status.tsv',
]

# per-root .bidsignore lines; logs/ is allowed by BIDS and needs no entry;
#   a folder is named without a trailing / (bids-validator ignores "foo/"
#   and "foo/**", but honors "foo" and "**/foo")
BIDSIGNORE = {
    'proc_rootdir': [
        '# the ABCD T1 atlas, the SpatialReference of atlas-space scans',
        'atlases',
    ] + STATUS_IGNORE + [
        '# motion estimates',
        '*_desc-motion_timeseries.tsv',
    ],
    'fsurf_rootdir': [
        "# FreeSurfer's own flat layout",
        'sub-*_ses-*',
        'fsaverage',
    ],
    'analysis_rootdir': [
        '# MATLAB structs (saved -v7)',
        '*.mat',
        '# custom suffixes, not yet listed: status, behavioral, runs, deconv,',
        '#   fstat, var, corr, corrmat, tseries, beta, sem, baseline, volinfo,',
        '#   motionqc, and the measure suffixes (thk, area, sulc, t1gm, dtiFA,',
        '#   rsiRND, ...)',
        '# surface masks, with no valid name in func/',
        '*_hemi-*_mask.func.gii',
        '*_hemi-*_mask.json',
    ],
    'rawqc_rootdir': SESSION_QC_IGNORE + STATUS_IGNORE + [
        '# QC images',
        '*.jpg',
    ],
    'postqc_rootdir': [
        '# the FreeSurfer QC, outside any datatype directory',
        'sub-*/ses-*/fsurf',
        '# session-root manual post QC files',
        'sub-*/ses-*/*_manualpostqc.tsv',
        'sub-*/ses-*/*_desc-manualpostqc*.tsv',
        'sub-*/ses-*/*_status.tsv',
        'sub-*/ses-*/*.lock',
    ] + STATUS_IGNORE + [
        '# QC images',
        '*.jpg',
        '# session mean images',
        '*_desc-meanNorm_boldref.nii.gz',
        '*_desc-meanNorm_boldref.json',
    ],
    'pc_rootdir': SESSION_QC_IGNORE,
}

# the dataset label of raw_rootdir in DatasetLinks (bids:rawdata:...)
RAW_LABEL = 'rawdata'

# per-root table: rootdir argument, dataset label, GeneratedBy kind,
#   source rootdirs, GeneratedBy Description
#   'cmig'  -> [cmig]
#   'fsurf' -> [freesurfer, cmig]
ROOTS = [
    ('proc_rootdir',     'cmig_proc',     'cmig',  ['raw_rootdir'],
     'preprocessing: cmig_proc_bids_anat, cmig_proc_bids_func, cmig_proc_bids_dwi'),
    ('fsurf_rootdir',    'freesurfer',    'fsurf', ['proc_rootdir'],
     'FreeSurfer recon-all, run by cmig_fsurf_bids'),
    ('analysis_rootdir', 'cmig_analysis', 'cmig',  ['proc_rootdir', 'fsurf_rootdir'],
     'analysis: cmig_analyze_bids_anat, cmig_analyze_bids_func, '
     'cmig_analyze_bids_dwi, cmig_analyze_bids_behav'),
    ('rawqc_rootdir',    'cmig_rawqc',    'cmig',  ['raw_rootdir'],
     'raw data QC: cmig_rawqc_bids, cmig_pc_bids_abcd, cmig_check_event_pcqc, '
     'cmig_manual_rawqc_bids'),
    ('postqc_rootdir',   'cmig_postqc',   'cmig',  ['proc_rootdir', 'fsurf_rootdir'],
     'post-processing QC: cmig_postqc_bids_anat, cmig_postqc_bids_dwi, '
     'cmig_postqc_bids_fmap, cmig_postqc_bids_fsurf, cmig_postqc_bids_func'),
    ('pc_rootdir',       'cmig_pc',       'cmig',  ['raw_rootdir'],
     'protocol compliance: cmig_pc_bids_abcd'),
]


def format_fs_version(FS_version):
    """Render a FreeSurfer version as a dotted string: 810 -> '8.1.0'.

    Values that are already dotted, or not three digits, are passed through.
    """
    vstr = str(FS_version)
    if vstr.isdigit() and len(vstr) == 3:
        return '.'.join(vstr)
    return vstr


def get_code_url(versions=None):
    """The cmig CodeURL: the tree of the release tag (versions'
    cmig_pipeline_tag, else $CMIG_PIPELINE_TAG), the repo root when there is
    none or it is 'local'."""
    versions = versions or {}
    tag = (versions.get('cmig_pipeline_tag')
           or os.environ.get(PIPELINE_TAG_ENV, '')).strip()
    if not tag or tag == 'local':
        return CODE_URL
    return TAG_CODE_URL.format(tag)


def get_container(versions=None):
    """The GeneratedBy Container: run from the image ($CMIG_DOCKER set
    there), Type singularity when singularity or apptainer runs it, else
    docker; on the host, singularity when versions (a .sif's versions.json)
    has an image_tag, else None. Tag is versions' image_tag, else
    $CMIG_IMAGE_TAG, omitted when empty or n/a."""
    versions = versions or {}
    image_tag = (versions.get('image_tag')
                 or os.environ.get(IMAGE_TAG_ENV, '')).strip()
    if image_tag == cmig_sing_utils.NA:
        image_tag = ''
    if os.environ.get('CMIG_DOCKER'):
        if (os.environ.get('SINGULARITY_CONTAINER')
                or os.environ.get('APPTAINER_CONTAINER')):
            ctype = 'singularity'
        else:
            ctype = 'docker'
    elif versions.get('image_tag'):
        ctype = 'singularity'
    else:
        return None
    container = collections.OrderedDict([('Type', ctype)])
    if image_tag:
        container['Tag'] = image_tag
    return container


def make_package_entries(versions, kind='cmig'):
    """One GeneratedBy entry (Name, Version) per package of versions
    (cmig_sing_utils.PACKAGES order); AFNI also gets the CodeURL of its
    commit. The fsurf roots skip FreeSurfer, which has its own entry."""
    entries = []
    for key, _, name in cmig_sing_utils.PACKAGES:
        if not versions.get(key) or (kind == 'fsurf' and key == 'freesurfer'):
            continue
        entry = collections.OrderedDict([('Name', name),
                                         ('Version', versions[key])])
        if key == 'afni' and versions.get('afni_commit'):
            entry['CodeURL'] = AFNI_CODE_URL.format(versions['afni_commit'])
        entries.append(entry)
    return entries


def make_generated_by(kind, pipeline_name=DEFAULT_PIPELINE_NAME,
                      software_name=DEFAULT_SOFTWARE_NAME,
                      software_version=DEFAULT_SOFTWARE_VERSION,
                      FS_version=DEFAULT_FS_VERSION, description=None,
                      container=None, versions=None):
    """Build the GeneratedBy array for one root.

    kind 'cmig': one entry, pipeline_name (default cmig) with Version
    (software_version, default the versions' cmig_pipeline, else
    $CMIG_PIPELINE_VERSION, omitted when unset), CodeURL (get_code_url()),
    Container (default get_container()) and Description (the step).
    kind 'fsurf': a freesurfer entry with FS_version first, then that entry.
    Then one entry per package of versions (an image's versions.json;
    None: $CMIG_VERSIONS_FILE, none without it).
    software_name is accepted for the callers' options and not written.
    """
    if versions is None:
        versions = cmig_sing_utils.read_versions()
    if software_version is None:
        software_version = (versions.get('cmig_pipeline')
                            or os.environ.get('CMIG_PIPELINE_VERSION'))
    if container is None:
        container = get_container(versions)
    entries = []
    if kind == 'fsurf':
        entries.append(collections.OrderedDict([
            ('Name', FS_NAME),
            ('Version', format_fs_version(FS_version)),
        ]))
    cmig = collections.OrderedDict([('Name', pipeline_name)])
    if software_version:
        cmig['Version'] = str(software_version)
    if description:
        cmig['Description'] = description
    cmig['CodeURL'] = get_code_url(versions)
    if container:
        cmig['Container'] = container
    entries.append(cmig)
    entries += make_package_entries(versions, kind)
    return entries


def check_root_name(rootdir, kind='cmig', pipeline_name=DEFAULT_PIPELINE_NAME):
    """BIDS requires the first GeneratedBy Name to be a substring of the
    root's directory name; return a warning message when it is not, else
    None."""
    name = FS_NAME if kind == 'fsurf' else pipeline_name
    dirname = os.path.basename(os.path.normpath(os.path.abspath(rootdir)))
    if name in dirname:
        return None
    return ('the directory name of {0} does not contain "{1}", the first '
            'GeneratedBy Name in its dataset_description.json, as BIDS '
            'requires').format(rootdir, name)


def make_root_url(source_dir, rootdir, relpath_flag=0):
    """The path of source_dir as rootdir's description names it, with a
    trailing /: absolute, or relative to rootdir with relpath_flag=1."""
    if relpath_flag:
        url = os.path.relpath(os.path.abspath(source_dir),
                              os.path.abspath(rootdir))
    else:
        url = os.path.abspath(source_dir)
    if not url.endswith('/'):
        url += '/'
    return url


def make_source_datasets(source_dirs, rootdir, relpath_flag=0):
    """Build the SourceDatasets array for one root.

    relpath_flag=0 writes absolute paths, correct when the roots may be in
    arbitrary locations. relpath_flag=1 writes paths relative to rootdir,
    which works when the roots are siblings under one master root dir and the
    dataset is later shared as a unit.
    """
    sources = []
    for source_dir in source_dirs:
        if not source_dir:
            continue
        sources.append(collections.OrderedDict(
            [('URL', make_root_url(source_dir, rootdir, relpath_flag))]))
    return sources


def get_atlas_url(versions=None):
    """The URL of the atlas package (ATLAS_ASSET_URL) for versions'
    atlases tag; None when there is none or it is 'local'."""
    versions = versions or {}
    tag = versions.get('atlases', '').strip()
    if not tag or tag == 'local':
        return None
    return ATLAS_ASSET_URL.format(tag)


def make_dataset_links(rootdirs, argname, relpath_flag=0, versions=None,
                       written=()):
    """Build the DatasetLinks object for one root: the dataset label of
    every other supplied root that exists or is among written (the roots
    this call describes, created as it goes) (rawdata for raw_rootdir) to
    its path, as SourceDatasets writes it; the proc root also gets
    cmig_pipeline_atlases (get_atlas_url), for its atlas sidecars'
    Sources."""
    labels = [('raw_rootdir', RAW_LABEL)] + [(r[0], r[1]) for r in ROOTS]
    links = collections.OrderedDict()
    for name, label in labels:
        if (name == argname or not rootdirs.get(name)
                or not (name in written or os.path.isdir(rootdirs[name]))):
            continue
        links[label] = make_root_url(rootdirs[name], rootdirs[argname],
                                     relpath_flag)
    if argname == 'proc_rootdir':
        atlas_url = get_atlas_url(versions)
        if atlas_url:
            links[cmig_write_bids_atlas.ATLAS_DATASET_LABEL] = atlas_url
    return links


def _write_file(fname_out, write_fn):
    """Write through a temporary file beside fname_out and rename, so
    per-session jobs starting together on a new root never leave a partly
    written file. Returns (errcode, errmsg)."""
    fname_tmp = None
    try:
        dirname = os.path.dirname(fname_out) or '.'
        os.makedirs(dirname, exist_ok=True)
        fname_tmp = cmig_write_bids_atlas.make_temp(fname_out)
        with open(fname_tmp, 'w') as fid:
            write_fn(fid)
        os.replace(fname_tmp, fname_out)
    except (OSError, IOError, TypeError, ValueError) as err:
        if fname_tmp and os.path.isfile(fname_tmp):
            os.remove(fname_tmp)
        return 1, 'failed to write {0}: {1}'.format(fname_out, err)
    return 0, None


def write_bidsignore(rootdir, lines, forceflag=0):
    """Write {rootdir}/.bidsignore from a list of lines (patterns and #
    comments). An existing file is left alone unless forceflag is set.

    Returns (errcode, errmsg).
    """
    fname_out = os.path.join(rootdir, FNAME_IGNORE)
    if os.path.isfile(fname_out) and not forceflag:
        return 0, None
    return _write_file(fname_out,
                       lambda fid: fid.write(''.join(l + '\n' for l in lines)))


def write_dataset_description(rootdir, name, bids_version=DEFAULT_BIDS_VERSION,
                              generated_by=None, source_datasets=None,
                              bidsignore=None, forceflag=0,
                              dataset_links=None):
    """Write dataset_description.json for a single derivative root.

    Knows nothing about which root this is: generated_by and source_datasets
    arrive fully formed, so extra or differently versioned entries need no
    change here.

    bidsignore, a list of lines, writes .bidsignore beside it (None: none).
    dataset_links, a dict of dataset label to path, is written as
    DatasetLinks (None or empty: none).

    An existing file is left alone unless forceflag is set, since the setup is
    idempotent and will be re-run, and silently rewriting a hand-edited
    description would be worse than skipping. Each file is checked on its own,
    so a root described before .bidsignore existed gets one.

    Returns (errcode, errmsg).
    """
    errcode, errmsg = 0, None
    if bidsignore is not None:
        errcode, errmsg = write_bidsignore(rootdir, bidsignore, forceflag)

    fname_out = os.path.join(rootdir, FNAME_DESC)
    if os.path.isfile(fname_out) and not forceflag:
        return errcode, errmsg

    info = collections.OrderedDict()
    info['Name'] = name
    info['BIDSVersion'] = bids_version
    info['DatasetType'] = 'derivative'
    info['GeneratedBy'] = generated_by if generated_by else []
    if source_datasets:
        info['SourceDatasets'] = source_datasets
    if dataset_links:
        info['DatasetLinks'] = dataset_links

    def write_fn(fid):
        json.dump(info, fid, indent=4)
        fid.write('\n')
    tmp_errcode, tmp_errmsg = _write_file(fname_out, write_fn)
    if tmp_errcode:
        return tmp_errcode, tmp_errmsg
    return errcode, errmsg


def write_dataset_descriptions(name, write_roots=None,
                               raw_rootdir=None, proc_rootdir=None,
                               fsurf_rootdir=None, analysis_rootdir=None,
                               rawqc_rootdir=None, postqc_rootdir=None,
                               pc_rootdir=None,
                               bids_version=DEFAULT_BIDS_VERSION,
                               pipeline_name=DEFAULT_PIPELINE_NAME,
                               software_name=DEFAULT_SOFTWARE_NAME,
                               software_version=DEFAULT_SOFTWARE_VERSION,
                               FS_version=DEFAULT_FS_VERSION,
                               relpath_flag=0, verbose=0, forceflag=0,
                               versions=None):
    """Write dataset_description.json for derivative roots.

    versions, the package versions added to GeneratedBy (an image's
    versions.json, read by the host for an image batch; None:
    $CMIG_VERSIONS_FILE, as in the image).

    Everything supplied is available as a source for SourceDatasets;
    write_roots names which roots get a file. It is a list of rootdir argument
    names, e.g. ['proc_rootdir', 'fsurf_rootdir']; the default None means every
    supplied root. A caller describing only its own output still passes the
    rootdirs its output derives from, so they can be named as sources.

    raw_rootdir is never described even if named in write_roots: it is the
    rawdata input, has its own description, and is only a source here.

    Each root gets a distinct Name, suffixed from the dataset label, since each
    root is a separate BIDS dataset, and its .bidsignore; the analysis root
    also gets the atlas-at_dseg.tsv lookup. Its DatasetLinks map the label of
    every other supplied root that exists to that root's path.

    Writes what it can and reports at the end rather than returning on the first
    failure, so one unwritable root does not leave the others undescribed.

    Returns (errcode, errmsg).
    """
    rootdirs = {
        'raw_rootdir':      raw_rootdir,
        'proc_rootdir':     proc_rootdir,
        'fsurf_rootdir':    fsurf_rootdir,
        'analysis_rootdir': analysis_rootdir,
        'rawqc_rootdir':    rawqc_rootdir,
        'postqc_rootdir':   postqc_rootdir,
        'pc_rootdir':       pc_rootdir,
    }

    if write_roots is not None:
        for argname in write_roots:
            if argname == 'raw_rootdir':
                print('WARNING: raw_rootdir is a source, not a derivative root; not described')
            elif not rootdirs.get(argname):
                print('WARNING: {0} named in write_roots but not supplied'.format(argname))

    errcode = 0
    errmsg = None
    nwritten = 0
    if versions is None:
        versions = cmig_sing_utils.read_versions()
    # the roots described here, which exist once written
    written = [r[0] for r in ROOTS if rootdirs.get(r[0])
               and (write_roots is None or r[0] in write_roots)]

    for argname, label, kind, source_names, description in ROOTS:
        rootdir = rootdirs.get(argname)
        if not rootdir:
            continue
        if write_roots is not None and argname not in write_roots:
            continue
        generated_by = make_generated_by(kind, pipeline_name, software_name,
                                         software_version, FS_version,
                                         description=description,
                                         versions=versions)
        source_datasets = make_source_datasets(
            [rootdirs.get(sname) for sname in source_names],
            rootdir, relpath_flag)
        tmp_errcode, tmp_errmsg = write_dataset_description(
            rootdir, '{0} - {1}'.format(name, label), bids_version,
            generated_by, source_datasets, bidsignore=BIDSIGNORE.get(argname),
            forceflag=forceflag,
            dataset_links=make_dataset_links(rootdirs, argname, relpath_flag,
                                             versions, written))
        if not tmp_errcode and argname == 'analysis_rootdir':
            tmp_errcode, tmp_errmsg = cmig_write_bids_atlas.write_dseg_lookup(
                rootdir, verbose=verbose, forceflag=forceflag)
        if tmp_errcode:
            errcode = tmp_errcode
            errmsg = tmp_errmsg
            print('WARNING: {0}'.format(tmp_errmsg))
        else:
            nwritten += 1
            if verbose:
                print('wrote {0}'.format(os.path.join(rootdir, FNAME_DESC)))

    if not nwritten and not errcode:
        errmsg = 'no derivative roots supplied'
        print('WARNING: {0}'.format(errmsg))

    return errcode, errmsg


def parse_args():
    parser = argparse.ArgumentParser(
        description='Write dataset_description.json for cmig BIDS derivative datasets.'
    )
    parser.add_argument('--name', required=True,
                        help='Dataset name; each root is suffixed with its label '
                             "(e.g. 'ABCD Study' -> 'ABCD Study - cmig_proc')")
    parser.add_argument('--write_roots', nargs='+', default=None,
                        help='Rootdir argument names to describe, e.g. proc_rootdir '
                             'fsurf_rootdir; others supplied serve only as sources '
                             '(default: every supplied root)')
    parser.add_argument('--raw_rootdir',
                        help='Root directory of BIDS rawdata; named as a source, not written to')
    parser.add_argument('--proc_rootdir',
                        help='Root directory of processed BIDS derivatives')
    parser.add_argument('--fsurf_rootdir',
                        help='Root directory of FreeSurfer derivatives')
    parser.add_argument('--analysis_rootdir',
                        help='Root directory of analysis output')
    parser.add_argument('--rawqc_rootdir',
                        help='Root directory of raw QC output')
    parser.add_argument('--postqc_rootdir',
                        help='Root directory of post-processing QC output')
    parser.add_argument('--pc_rootdir',
                        help='Root directory of protocol compliance output')
    parser.add_argument('--bids_version', default=DEFAULT_BIDS_VERSION,
                        help='BIDS specification version (default: %(default)s)')
    parser.add_argument('--pipeline_name', default=DEFAULT_PIPELINE_NAME,
                        help='Name of the pipeline that produced the output, the '
                             'first GeneratedBy Name (default: %(default)s)')
    parser.add_argument('--software_name', default=DEFAULT_SOFTWARE_NAME,
                        help='Not written; accepted for compatibility '
                             '(default: %(default)s)')
    parser.add_argument('--software_version', default=DEFAULT_SOFTWARE_VERSION,
                        help='Version of the pipeline '
                             '(default: $CMIG_PIPELINE_VERSION)')
    parser.add_argument('--FS_version', type=int, default=DEFAULT_FS_VERSION,
                        help='FreeSurfer version number (default: %(default)s)')
    parser.add_argument('--relpath_flag', type=int, default=0, choices=[0, 1],
                        help='Write SourceDatasets and DatasetLinks paths '
                             'relative to each root, assuming one master root '
                             'dir (default: 0)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='Display status messages (default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing output (default: 0)')
    return parser.parse_args()


def run(args):
    return write_dataset_descriptions(
        args.name,
        write_roots=args.write_roots,
        raw_rootdir=args.raw_rootdir,
        proc_rootdir=args.proc_rootdir,
        fsurf_rootdir=args.fsurf_rootdir,
        analysis_rootdir=args.analysis_rootdir,
        rawqc_rootdir=args.rawqc_rootdir,
        postqc_rootdir=args.postqc_rootdir,
        pc_rootdir=args.pc_rootdir,
        bids_version=args.bids_version,
        pipeline_name=args.pipeline_name,
        software_name=args.software_name,
        software_version=args.software_version,
        FS_version=args.FS_version,
        relpath_flag=args.relpath_flag,
        verbose=args.verbose,
        forceflag=args.forceflag,
    )


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
