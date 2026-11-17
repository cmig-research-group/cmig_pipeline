#!/usr/bin/env python3
"""
cmig_proc_bids_dwi.py

dMRI processing of one session: runs the dwi processing pipeline
(bash/cmig_proc_dwi.bash: run and field map selection, gradwarp, TORTOISE
correction, Synb0-DISCO where there is no usable field map, registration to
the processed T1w), checks its outputs and writes the status file.

Usage:
  cmig_proc_bids_dwi.py --participant_id sub-X --session_id ses-Y
      --raw_rootdir DIR --output_rootdir DIR [--proc_anat_rootdir DIR]
      [--fs_license FILE] [--n_cpus N] [--output_res MM]
      [--output_size XxYxZ] [--output_orientation XYZ] [--max_rotation DEG]
      [--proc_dir DIR] [--cleanup_flag 0|1] [--forceflag 0|1]

Outputs, in {output_rootdir}/{participant_id}/{session_id}/:
  dwi/   the processed dwi series (run-NN_desc-preproc_dwi, with .bval, .bvec,
         .json, and run-NN_desc-motion_timeseries.tsv), and the status file
         {pid}_{sid}_desc-proc_status.tsv
  fmap/  the B0dx displacement field (.nii.gz, .json) and the corrected AP
         field map (desc-preproc_epi)
  dwi/proc_dwi/ (or {proc_dir}/{pid}_{sid}_proc_dwi with --proc_dir)
         intermediate files and logs; with cleanup_flag 1, deleted after a
         session whose outputs pass the check, once its benchmark/stages.tsv
         is printed (TORTOISE's version and the field map's volume and shell
         counts are in the dwi sidecar)

A session that cannot be processed (no usable dwi run, field map or sMRI) is
not run; its status file reports errcode 1 and the reason.

The tools (FSL, FreeSurfer, ANTs, c3d, TORTOISE, Synb0-DISCO) come from the
environment, as the CMIG BIDS app image sets it; GRADIENT_DISTORTION_FIELDS_DIR
and CMIG_PIPELINE_ATLASDIR must be set.

Created:  11/02/2026 by Don Hagler
Prev mod: 11/02/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Last mod: 11/12/2026 by Don Hagler

Based on bids_app.py, by Diliana Pecheva (abcd_preproc 81210b9)
"""

import argparse
import datetime
import os
import shutil
import signal
import subprocess
import sys
import threading

import cmig_sing_utils
import cmig_write_status as cmig_status
from cmig_check_bids_dwi_inputs import input_check, not_processed_msg, processable
from cmig_check_bids_dwi_outputs import MAX_ROTATION, PROC_DIR, check_session

HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE = os.path.join(HERE, '..', 'bash', 'cmig_proc_dwi.bash')
MODULE = 'desc-proc'
LOG_NAME = 'cmig_proc_bids_dwi.log'   # a session's complete log, in {proc_dir}/logs/
STAGES_HEADER = ('stage\tstart\ttool\texit\twall_s\tuser_s\tsys_s\tmaxrss_kb\tcpu_pct'
                 '\tmajflt\tminflt\tfs_in\tfs_out\tcommand\n')


def build_parser():
    parser = argparse.ArgumentParser(
        description='dMRI processing of one session (cmig_proc_dwi.bash)')
    parser.add_argument('--participant_id', required=True,
                        help='participant ID (e.g. sub-XXXXXXXX)')
    parser.add_argument('--session_id', required=True,
                        help='session ID (e.g. ses-00A)')
    parser.add_argument('--raw_rootdir', required=True,
                        help='BIDS rawdata root (dwi, field maps, scans.tsv)')
    parser.add_argument('--output_rootdir', required=True,
                        help='root of the outputs (the proc root)')
    parser.add_argument('--proc_anat_rootdir', default=None,
                        help='root of the processed T1w (anat/); default: output_rootdir')
    parser.add_argument('--fs_license', default=None,
                        help='FreeSurfer licence file; default: $FS_LICENSE')
    parser.add_argument('--n_cpus', type=int, default=None,
                        help='threads for the tools; default: the CPUs this process may '
                             'run on (its affinity); give it where the scheduler does not '
                             'confine a job to its CPUs (e.g. SGE)')
    parser.add_argument('--output_res', default='1.7',
                        help='output voxel size (mm) (default: %(default)s)')
    parser.add_argument('--output_size', default='140x140x140',
                        help='output grid (default: %(default)s)')
    parser.add_argument('--output_orientation', default='LAS',
                        help='direction each output axis points towards, as TORTOISE\'s '
                             '--output_orientation; the second letter A or P '
                             '(default: %(default)s)')
    parser.add_argument('--max_rotation', type=float, default=MAX_ROTATION,
                        help='output check: largest rotation (degrees) between a series\' '
                             'gradient directions and the raw data\'s (default: %(default)s)')
    parser.add_argument('--proc_dir', default=None,
                        help='folder for the intermediate files, in a session subfolder '
                             '{pid}_{sid}_proc_dwi (e.g. a job\'s local scratch); default: '
                             'dwi/proc_dwi/ in the session\'s output folder')
    parser.add_argument('--cleanup_flag', type=int, default=1, choices=[0, 1],
                        help='1: after a session whose outputs pass the check, print '
                             'benchmark/stages.tsv and delete the intermediate folder '
                             '(default: %(default)s)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='rerun even if the status file reports success and the '
                             'outputs pass the check (default: %(default)s)')
    return parser


def parse_args():
    return build_parser().parse_args()


def run(args):
    return proc_dwi(
        args.participant_id, args.session_id, args.raw_rootdir, args.output_rootdir,
        proc_anat_rootdir=args.proc_anat_rootdir, fs_license=args.fs_license,
        n_cpus=args.n_cpus, output_res=args.output_res, output_size=args.output_size,
        output_orientation=args.output_orientation, max_rotation=args.max_rotation,
        proc_dir=args.proc_dir, cleanup_flag=args.cleanup_flag, forceflag=args.forceflag)


def _now():
    now = datetime.datetime.now()
    return now.strftime('%Y-%m-%dT%H:%M:%S.') + '{0:03d}'.format(now.microsecond // 1000)


def _status_ok(session_dir, participant_id, session_id):
    """True if the status file exists and reports errcode 0."""
    fpath = cmig_status.status_fpath(os.path.join(session_dir, 'dwi'),
                                     participant_id, session_id, MODULE)
    try:
        with open(fpath) as f:
            lines = f.read().splitlines()
        row = dict(zip(lines[0].split('\t'), lines[1].split('\t')))
        return row.get('errcode') == '0'
    except (OSError, IndexError):
        return False


def _write_status(session_dir, participant_id, session_id, errcode, errmsg, startdate):
    """Write dwi/{pid}_{sid}_desc-proc_status.tsv and print its block."""
    status_file = cmig_status.cmig_write_status(
        participant_id, session_id, os.path.join(session_dir, 'dwi'), MODULE,
        errcode, errmsg, startdate=startdate)
    # this module is its own wrapper, so it prints both halves of the
    #   finished pair cmig_check_batchdir.py reads, as cmig_fsurf_bids does
    if errcode == 0:
        print('cmig_proc_bids_dwi: finished')
        print('cmig_proc_bids_dwi: completed successfully')
    cmig_sing_utils.add_status_versions(status_file, app_name='cmig_proc_bids_dwi')


class _Tee:
    """In place of sys.stdout or sys.stderr while a session runs: everything
    written goes both to the job's output and to the session's log."""

    def __init__(self, stream, log, lock):
        self.stream, self.log, self.lock = stream, log, lock

    def write(self, text):
        with self.lock:
            self.stream.write(text)
            self.stream.flush()
            self.log.write(text)
            self.log.flush()
        return len(text)

    def flush(self):
        self.stream.flush()


def _copy_lines(src, dst):
    for line in src:
        dst.write(line)


def _print_stages(bench_tsv):
    """Print the step timings (benchmark/stages.tsv) to stdout, as they go with
    the intermediate folder."""
    try:
        with open(bench_tsv) as f:
            text = f.read()
    except OSError:
        print('cmig_proc_bids_dwi: no {0}'.format(bench_tsv))
        return
    print('cmig_proc_bids_dwi: stages.tsv:')
    sys.stdout.write(text)
    sys.stdout.flush()


def proc_dwi(participant_id, session_id, raw_rootdir, output_rootdir,
             proc_anat_rootdir=None, fs_license=None, n_cpus=None,
             output_res='1.7', output_size='140x140x140', output_orientation='LAS',
             max_rotation=MAX_ROTATION, proc_dir=None, cleanup_flag=1, forceflag=0):
    """Run cmig_proc_dwi.bash for one session, check its outputs and write the
    status file. Everything printed, the pipeline's output and errors included,
    is also written to {proc_dir}/logs/cmig_proc_bids_dwi.log; with cleanup_flag
    1, the intermediate folder (proc_dir) is deleted after a session whose
    outputs pass the check.

    Returns: (errcode, errmsg)
    """
    if proc_anat_rootdir is None:
        proc_anat_rootdir = output_rootdir
    raw_rootdir = os.path.abspath(raw_rootdir)
    output_rootdir = os.path.abspath(output_rootdir)
    proc_anat_rootdir = os.path.abspath(proc_anat_rootdir)
    session_dir = os.path.join(output_rootdir, participant_id, session_id)
    # a session subfolder of a given proc_dir, so that cleanup deletes only it
    if proc_dir:
        proc_dir = os.path.join(os.path.abspath(proc_dir),
                                '{0}_{1}_proc_dwi'.format(participant_id, session_id))
    else:
        proc_dir = os.path.join(session_dir, PROC_DIR)

    # a session that cannot be processed gets a failed status, so the rounds
    #   count it (and stop after max_attempts) rather than resubmit it
    have, sel = input_check(raw_rootdir, proc_anat_rootdir, participant_id, session_id)
    if not processable(have, sel):
        errmsg = not_processed_msg(have, sel)
        print('cmig_proc_bids_dwi: {0} {1}: {2}'.format(participant_id, session_id, errmsg))
        _write_status(session_dir, participant_id, session_id, 1, errmsg, _now())
        return 1, errmsg

    if not forceflag and _status_ok(session_dir, participant_id, session_id) \
            and not check_session(session_dir, max_rotation, proc_dir):
        print('cmig_proc_bids_dwi: {0} {1}: outputs complete; skipped'.format(
            participant_id, session_id))
        return 0, None

    if not fs_license:
        fs_license = os.environ.get('FS_LICENSE', '')
    if not fs_license or not os.path.isfile(fs_license):
        errmsg = 'FreeSurfer licence not found: set FS_LICENSE or pass --fs_license'
        print('cmig_proc_bids_dwi: ERROR: {0}'.format(errmsg))
        return 1, errmsg

    log_dir = os.path.join(proc_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    with open(os.path.join(log_dir, LOG_NAME), 'w') as log:
        lock = threading.Lock()
        out, err = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = _Tee(out, log, lock), _Tee(err, log, lock)
        try:
            errcode, errmsg = _run_pipeline(
                participant_id, session_id, raw_rootdir, proc_anat_rootdir,
                output_rootdir, session_dir, proc_dir, os.path.abspath(fs_license),
                n_cpus, output_res, output_size, output_orientation, max_rotation)
        finally:
            sys.stdout, sys.stderr = out, err
    # after the log is closed, so no open file is left in the folder deleted
    if cleanup_flag and not errcode:
        _print_stages(os.path.join(proc_dir, 'benchmark', 'stages.tsv'))
        shutil.rmtree(proc_dir)
        print('cmig_proc_bids_dwi: cleanup: deleted {0}'.format(proc_dir))
    return errcode, errmsg


def _run_pipeline(participant_id, session_id, raw_rootdir, proc_anat_rootdir,
                  output_rootdir, session_dir, proc_dir, fs_license, n_cpus, output_res,
                  output_size, output_orientation, max_rotation):
    startdate = _now()
    bench_dir = os.path.join(proc_dir, 'benchmark')
    os.makedirs(bench_dir, exist_ok=True)
    bench_tsv = os.path.join(bench_dir, 'stages.tsv')
    with open(bench_tsv, 'w') as f:
        f.write(STAGES_HEADER)

    if not n_cpus:
        n_cpus = len(os.sched_getaffinity(0))
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS=str(n_cpus),
               ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=str(n_cpus),
               GW_DIR=cmig_sing_utils.cmig_getenv('GRADIENT_DISTORTION_FIELDS_DIR'),
               BENCH_STAGE_TSV=bench_tsv)
    cmd = ['bash', os.path.abspath(PIPELINE), '--study=ABCD',
           '--prf={0}'.format(raw_rootdir), '--smri_dir={0}'.format(proc_anat_rootdir),
           '--sub_id={0}'.format(participant_id[len('sub-'):]),
           '--ses_id={0}'.format(session_id[len('ses-'):]),
           '-o={0}'.format(output_rootdir), '--gpu=0',
           '--output_res={0}'.format(output_res), '--output_size={0}'.format(output_size),
           '--output_orientation={0}'.format(output_orientation),
           '--fs_license={0}'.format(fs_license), '--proc_dir={0}'.format(proc_dir)]
    print('cmig_proc_bids_dwi: {0}'.format(' '.join(cmd)))
    sys.stdout.flush()
    # in its own process group, so a stop reaches every tool it started; from
    #   proc_dir, so files a tool writes to its current folder stay there
    proc = subprocess.Popen(cmd, cwd=proc_dir, env=env, start_new_session=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, errors='replace')
    copiers = [threading.Thread(target=_copy_lines, args=(proc.stdout, sys.stdout)),
               threading.Thread(target=_copy_lines, args=(proc.stderr, sys.stderr))]
    for t in copiers:
        t.start()

    stopped = []

    def stop(signum, frame):
        stopped.append(signum)
        os.killpg(proc.pid, signal.SIGTERM)

    old = {s: signal.signal(s, stop) for s in (signal.SIGTERM, signal.SIGUSR1)}
    try:
        rc = proc.wait()
    finally:
        for s, h in old.items():
            signal.signal(s, h)
    for t in copiers:
        t.join()
    print('cmig_proc_bids_dwi: pipeline exit status: {0}'.format(rc))

    problems = check_session(session_dir, max_rotation, proc_dir)
    if stopped:
        problems = ['stopped by signal {0}'.format(stopped[0])] + problems
    for q in problems:
        print('cmig_proc_bids_dwi: OUTPUT CHECK: {0}'.format(q))
    errcode = int(bool(problems))
    errmsg = '; '.join(problems) if problems else None
    _write_status(session_dir, participant_id, session_id, errcode, errmsg, startdate)
    return errcode, errmsg


def main():
    errcode, errmsg = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
