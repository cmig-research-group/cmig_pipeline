#!/usr/bin/env python3
"""
python wrapper for FreeSurfer's recon-all with input from BIDS anat

Created:  11/14/2025 by Don Hagler
Prev mod: 11/17/2025 by Don Hagler
Prev mod: 03/06/2026 by Don Hagler
Prev mod: 05/06/2026 by Don Hagler
Prev mod: 06/02/2026 by Don Hagler
Prev mod: 06/17/2026 by Don Hagler
Prev mod: 08/10/2026 by Don Hagler
Prev mod: 08/20/2026 by Don Hagler
Prev mod: 09/02/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 11/09/2026 by Don Hagler
Last mod: 11/11/2026 by Don Hagler
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

import cmig_sing_utils
import cmig_write_status as cmig_status

DEFAULT_FS_VERSION = 820
DEFAULT_OPTIONS = '-all'


def build_parser():
    parser = argparse.ArgumentParser(
        description='Run FreeSurfer recon-all on BIDS anat data.'
    )
    parser.add_argument('--participant_id', required=True,
                        help='Participant ID (e.g., sub-XXXXXXXX)')
    parser.add_argument('--session_id', required=True,
                        help='Session ID (e.g., ses-00A)')
    parser.add_argument('--proc_rootdir', required=True,
                        help='Root directory of processed anat data (read-only input)')
    parser.add_argument('--output_rootdir', required=True,
                        help='Root directory for FreeSurfer output (created if absent)')
    parser.add_argument('--T2_flag', type=int, default=0, choices=[0, 1, 2],
                        help='T2w handling: 0=ignore, 1=use if available, 2=require (default: 0)')
    parser.add_argument('--options', default=DEFAULT_OPTIONS,
                        help='recon-all options (default: %(default)s)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing FreeSurfer output (default: 0)')
    parser.add_argument('--FS_version', type=int, default=DEFAULT_FS_VERSION,
                        help='FreeSurfer version number (default: %(default)s)')
    parser.add_argument('--setup_script', default=None,
                        help='Path to SetUpFreeSurfer.csh (default: $PUBSW/bin/SetUpFreeSurfer.csh)')
    return parser


def parse_args():
    return build_parser().parse_args()


def run(args):
    if not os.path.isdir(args.proc_rootdir):
        sys.exit(f'Error: directory not found: {args.proc_rootdir}')

    os.makedirs(args.output_rootdir, exist_ok=True)

    return fsurf_recon(
        args.participant_id, args.session_id,
        args.proc_rootdir, args.output_rootdir,
        FS_version=args.FS_version,
        T2_flag=args.T2_flag,
        options=args.options,
        forceflag=bool(args.forceflag),
        setup_script=args.setup_script,
    )


def _log_errmsg(fsdir, n_keep=3):
    """Return the recon-all.log lines that name the cause of a failure.

    recon-all's stderr names the step that failed, e.g. 'threshold mask volume at
    0.001', while the cause is only in the log, e.g. an invalid license key read
    by mris_volmask. Matched at the line start: the log also carries benign lines
    with ERROR mid-line, such as the LBFGS 'ERROR RETURN OF LINE SEARCH' messages.
    """
    fname_log = os.path.join(fsdir, 'scripts', 'recon-all.log')
    if not os.path.isfile(fname_log):
        return ''
    lines = []
    try:
        with open(fname_log, errors='replace') as f:
            for line in f:
                if line.startswith('ERROR'):
                    line = line.strip()
                    if line not in lines:
                        lines.append(line)
    except OSError:
        print(f'WARNING: failed to read {fname_log}')
        return ''
    if not lines:
        return ''
    return 'recon-all.log: ' + ' | '.join(lines[-n_keep:])


def _write_status(fsdir, participant_id, session_id, errcode, errmsg,
                  FS_version=None, startdate=None, procflag=1):
    """Write the fsurf status row and print its block.

    The row and the printed block come from cmig_write_status.py, which the
    MATLAB modules' cmig_write_status.m mirrors, so the two stay comparable.
    procflag 0 leaves an existing status file as it is.
    """
    status_file = cmig_status.cmig_write_status(
        participant_id, session_id, fsdir, 'desc-fsurf', errcode, errmsg,
        startdate=startdate, procflag=procflag)
    # the finished pair cmig_check_batchdir.py reads to name a block that
    #   succeeded. it must sit directly after errmsg, and only a success prints
    #   it. for a MATLAB module the finished half comes from the module and the
    #   completed half from cmig_sing_utils.check_status; this module is its own
    #   wrapper, so it prints both
    if errcode == 0:
        print('cmig_fsurf_bids: finished')
        print('cmig_fsurf_bids: completed successfully')
    # recon-all did the work here, not a compiled app, so it is what gets named
    cmig_sing_utils.add_status_versions(
        status_file, app_name='recon-all',
        app_version=FS_version if FS_version else None)


def fsurf_recon(participant_id, session_id, proc_rootdir, fsurf_rootdir,
                FS_version=DEFAULT_FS_VERSION, T2_flag=0, options=DEFAULT_OPTIONS,
                forceflag=False, setup_script=None):
    """
    Runs FreeSurfer recon-all on BIDS data.
    Returns: (errcode, errmsg)
    """
    errcode = 0
    errmsg = None
    startdate = cmig_status.module_start()

    # setup script default
    if setup_script is None:
        pubsw = os.getenv('PUBSW', '')
        setup_script = f'{pubsw}/bin/SetUpFreeSurfer.csh'

    # get input files from BIDS
    fnames_t1w = get_bids_data(participant_id, session_id, proc_rootdir, 'anat', 'T1w')
    fnames_t2w = get_bids_data(participant_id, session_id, proc_rootdir, 'anat', 'T2w')

    subjid = f'{participant_id}_{session_id}'
    fsdir = os.path.join(fsurf_rootdir, subjid)

    if not fnames_t1w:
        errmsg = f'T1w input not found for {participant_id} {session_id} in {proc_rootdir}'
        print(f'ERROR: {errmsg}')
        _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                      startdate)
        return 1, errmsg

    if not fnames_t2w and T2_flag == 2:
        errmsg = f'T2w input not found for {participant_id} {session_id} in {proc_rootdir}'
        print(f'ERROR: {errmsg}')
        _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                      startdate)
        return 1, errmsg
    elif fnames_t2w and T2_flag == 1:
        print(f'NOTE: T2w input found for {participant_id} {session_id} in {proc_rootdir}')

    # check existing directory
    if os.path.isdir(fsdir):
        if not check_fsurf_files(fsdir):
            print('WARNING: removing incomplete FreeSurfer recon')
            try:
                shutil.rmtree(fsdir)
            except Exception as e:
                errmsg = f'failed to remove incomplete FreeSurfer recon: {e}'
                print(f'WARNING: {errmsg}')
                _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                              startdate)
                return 1, errmsg
        elif forceflag:
            print('WARNING: removing existing FreeSurfer recon')
            try:
                shutil.rmtree(fsdir)
            except Exception as e:
                errmsg = f'failed to remove existing FreeSurfer recon: {e}'
                print(f'WARNING: {errmsg}')
                _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                              startdate)
                return 1, errmsg
        else:
            print(f'NOTE: complete FreeSurfer recon dir:\n{fsdir}')
            # a skip keeps an errcode 0 status, so its procdate stays the time
            #   the recon ran and the steps that read the recon are not stale
            prior = cmig_status.read_status(cmig_status.status_fpath(
                fsdir, participant_id, session_id, 'desc-fsurf'))
            procflag = int(not (prior and (prior.get('errcode') or '').strip() == '0'))
            _write_status(fsdir, participant_id, session_id, 0, None, FS_version,
                          startdate, procflag)
            return 0, None

    # build command string
    cmd  = f'source {setup_script} {FS_version}\n'
    cmd += f'setenv SUBJECTS_DIR {fsurf_rootdir}\n'
    cmd += f'recon-all -subjid {subjid}'
    for fname in fnames_t1w:
        cmd += f' -i {fname}'
    if T2_flag and fnames_t2w:
        cmd += f' -T2 {fnames_t2w[0]}'
    cmd += f' {options}'

    # start process
    print('Running FreeSurfer...')
    process = subprocess.Popen(
        cmd,
        shell=True,
        executable='/bin/csh',
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        bufsize=1
    )

    # stream stdout line by line
    for line in process.stdout:
        print(line, end='')

    # collect stderr after process finishes
    stderr_output = process.stderr.read()

    # wait for the process to finish
    process.wait()

    # report error if non-zero return code
    if process.returncode != 0:
        errmsg = f'FreeSurfer command ended in error: {stderr_output}'
        # recon-all reports the step that failed, not why; the reason is in the log
        log_errmsg = _log_errmsg(fsdir)
        if log_errmsg:
            errmsg = f'{errmsg}\n{log_errmsg}'
        print(f'ERROR: {errmsg}')
        _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                      startdate)
        return 1, errmsg

    # check output
    if not os.path.isdir(fsdir):
        errmsg = 'no FreeSurfer output'
        print(f'ERROR: {errmsg}')
        _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                      startdate)
        return 1, errmsg

    if not check_fsurf_files(fsdir):
        errmsg = 'incomplete FreeSurfer output'
        print(f'ERROR: {errmsg}')
        _write_status(fsdir, participant_id, session_id, 1, errmsg, FS_version,
                      startdate)
        return 1, errmsg

    _write_status(fsdir, participant_id, session_id, 0, None, FS_version,
                  startdate)
    return 0, None


def check_fsurf_files(fsdir):
    """Check whether FreeSurfer output is complete."""
    # NOTE: the same 22 files cmig_check_fsurf_files.m requires. The two lists had
    #   drifted apart, 7 here against 21 there, and this is the check on the native
    #   path, so a recon could pass during processing and then be judged incomplete
    #   at analysis, where cmig_check_bids_fsurf sets fs_complete_flag from the
    #   MATLAB list. ribbon.mgz is needed by the dwi white matter bias correction
    volumes = ['aseg.mgz', 'ribbon.mgz', 'orig.mgz', 'nu.mgz']
    surfaces = ['orig', 'white', 'pial', 'sphere.reg', 'thickness', 'area', 'sulc']
    annots = ['aparc.annot', 'aparc.a2009s.annot']
    required_files = [f'{fsdir}/mri/{v}' for v in volumes]
    for hemi in ('lh', 'rh'):
        required_files += [f'{fsdir}/surf/{hemi}.{s}' for s in surfaces]
        required_files += [f'{fsdir}/label/{hemi}.{a}' for a in annots]
    return all(os.path.isfile(f) for f in required_files)


def get_bids_data(participant_id, session_id, proc_rootdir, modality, suffix):
    """
    Retrieve BIDS filenames.

    Parameters
    ----------
    participant_id : str
        BIDS subject ID (with 'sub-' prefix).
    session_id : str
        BIDS session ID (with 'ses-' prefix).
    proc_rootdir : str
        Root directory of the BIDS dataset.
    modality : str
        BIDS datatype (e.g., 'func', 'anat', 'dwi').
    suffix : str
        BIDS suffix (e.g., 'bold', 'T1w').

    Returns
    -------
    list of str
        Matching file paths.
    """
    subj_dir = os.path.join(proc_rootdir, participant_id)
    if session_id:
        subj_dir = os.path.join(subj_dir, session_id)

    modality_dir = os.path.join(subj_dir, modality)
    if not os.path.exists(modality_dir):
        return []

    files = []
    for pattern in [f'*{suffix}.nii', f'*{suffix}.nii.gz']:
        files.extend(glob.glob(os.path.join(modality_dir, pattern)))

    return files


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
