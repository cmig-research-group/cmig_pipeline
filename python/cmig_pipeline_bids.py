#!/usr/bin/env python3
"""
Pipeline for one session, anat, func and dwi together:
  cmig_analyze_bids_behav -> cmig_proc_bids_anat -> cmig_proc_bids_func
  -> cmig_proc_bids_dwi -> cmig_fsurf_bids -> cmig_analyze_bids_anat
  -> [for each task] cmig_analyze_bids_func -> cmig_analyze_bids_dwi
  -> cmig_postqc_bids_anat -> cmig_postqc_bids_fsurf -> cmig_postqc_bids_func
  -> cmig_postqc_bids_fmap (acq-func) -> cmig_postqc_bids_dwi
  -> cmig_postqc_bids_fmap (acq-dwi)

The MATLAB modules run through cmig_run_bids_module.py. The func steps (the
behavioral step, proc_func, analyze_func, postqc_func and postqc_fmap
acq-func) run when --func_flag is 1 (the default) and the session has task
fMRI data in rawdata; the dwi steps (proc_dwi, analyze_dwi, postqc_dwi and
postqc_fmap acq-dwi) the same way, with --dwi_flag and dwi data in rawdata.
cmig_proc_bids_dwi needs the dwi tools (the CMIG BIDS app image) and 2 or more
CPUs (Synb0-DISCO's inference); it runs its tools with --dwi_n_cpus threads
(default 8), so the thread count does not depend on the host; its
intermediate folder is {tmpdir}/{pid}_{sid}_proc_dwi, tmpdir the job's
(SCRATCHROOT, else TMPDIR or /tmp).
--analyze_anat_flag 0 skips cmig_analyze_bids_anat. --forceflag goes to every
step, --verbose to those that take it (cmig_analyze_bids_behav,
cmig_analyze_bids_dwi and the postqc modules).

A failed step skips only the steps that need its output:
  cmig_proc_bids_anat: cmig_fsurf_bids and everything needing it
    (analyze_anat, postqc_anat, postqc_fsurf, analyze_func), and the dwi
    steps (proc_dwi registers to the processed T1w); proc_func still runs,
    without registration to T1, and postqc_func and postqc_fmap
  cmig_fsurf_bids: analyze_anat, postqc_fsurf and analyze_func (which samples
    the surfaces); postqc_anat, postqc_func and postqc_fmap still run
  cmig_proc_bids_func: analyze_func, postqc_func and postqc_fmap acq-func
  cmig_proc_bids_dwi: analyze_dwi, postqc_dwi and postqc_fmap acq-dwi;
    cmig_fsurf_bids failing skips analyze_dwi only
postqc_fmap acq-func is also skipped when proc_func processed no acq-func
field map and rawdata has no usable one (proc_func then succeeds without B0
unwarping); postqc_fmap acq-dwi likewise for acq-dwi.
Any other failure is recorded and the rest runs; a failed analyze_func task
is not counted while another task succeeded. errcode is 1 if anything counted failed, with every
message joined in errmsg.

Created:  10/05/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/10/2026 by Don Hagler
Prev mod: 10/12/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/02/2026 by Don Hagler
Prev mod: 11/03/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/13/2026 by Don Hagler
"""

import argparse
import os
import sys
from argparse import Namespace

import cmig_get_bids_data
import cmig_analyze_bids_behav
import cmig_select_bids_func_runs
import cmig_fsurf_bids
import cmig_proc_bids_dwi
import cmig_run_bids_module

DEFAULT_APPVER = 'mmps_dev'
DEFAULT_APPDIR = '/usr/pubsw/packages/MMPS/apps'
DEFAULT_PUBSDIR = '/usr/pubsw'
DEFAULT_MRTDIR = '/usr/pubsw/packages/matlab/R2023a'
DEFAULT_MRTDIR_TARGET = '/opt/matlabruntime/R2023a'


def skip_postqc_fmap(pid, sid, proc_rootdir, raw_rootdir, acq_type='func'):
    """True when the proc step processed no field map of acq_type (func or
    dwi) and rawdata has no usable one, so cmig_postqc_bids_fmap would have
    nothing to check; the one test of it, for this pipeline and
    cmig_mkjobs_bids_pipeline.py."""
    for rootdir in (proc_rootdir, raw_rootdir):
        fnames, _ = cmig_get_bids_data.get_bids_data(
            pid, sid, rootdir, modality='fmap', series_type='_epi',
            acq_type=acq_type)
        if fnames:
            return False
    return True


def has_dwi(pid, sid, raw_rootdir):
    """True when rawdata has a usable dwi series for the session; the one test
    of it, for this pipeline and cmig_mkjobs_bids_pipeline.py."""
    fnames, _ = cmig_get_bids_data.get_bids_data(
        pid, sid, raw_rootdir, modality='dwi', series_type='_dwi')
    return bool(fnames)


def _demarcate(module):
    print('\n' + '=' * 60)
    print(module)
    print('=' * 60)


def build_parser():
    parser = argparse.ArgumentParser(
        description='Run the BIDS pipeline for one session: anat (proc, '
                    'FreeSurfer, analysis, postqc), func (behavioral, proc, '
                    'analysis per task, postqc) and dwi (proc, analysis, postqc).'
    )
    parser.add_argument('--participant_id', required=True,
                        help='Participant ID (e.g., sub-G9MENY7V)')
    parser.add_argument('--session_id', required=True,
                        help='Session ID (e.g., ses-00A)')
    parser.add_argument('--raw_rootdir', required=True,
                        help='Root directory of raw BIDS data (read-only input)')
    parser.add_argument('--proc_rootdir', required=True,
                        help='Root directory for processed output (created if absent)')
    parser.add_argument('--fsurf_rootdir', required=True,
                        help='Root directory for FreeSurfer output (created if absent)')
    parser.add_argument('--analysis_rootdir', required=True,
                        help='Root directory for analysis output (created if absent)')
    parser.add_argument('--postqc_rootdir', required=True,
                        help='Root directory for post-processing QC output (created if absent)')
    parser.add_argument('--func_flag', type=int, default=1, choices=[0, 1],
                        help='Run the func steps, for a session with task fMRI data '
                             '(default: 1)')
    parser.add_argument('--dwi_flag', type=int, default=1, choices=[0, 1],
                        help='Run the dwi steps, for a session with dwi data '
                             '(default: 1)')
    parser.add_argument('--dwi_n_cpus', type=int, default=8,
                        help='Threads for the dwi processing tools (default: 8)')
    parser.add_argument('--analyze_anat_flag', type=int, default=1, choices=[0, 1],
                        help='Run cmig_analyze_bids_anat (default: 1)')
    parser.add_argument('--task', nargs='+', default=None,
                        help='Task labels to analyze (default: whatever the session has '
                             'usable data for in rawdata)')
    parser.add_argument('--ignore_delay_tasks', nargs='+', default=[],
                        help='Tasks whose REDCap eprime_mismatch_ok flag is set, '
                             'exempt from the eprime-scan delay check (default: none)')
    parser.add_argument('--exclude_noncompliant_flag', type=int, default=0, choices=[0, 1],
                        help='The run selection treats runs with compliant 0 in the '
                             'rawdata _scans.tsv as invalid (default: 0)')
    parser.add_argument('--pydir', default=None,
                        help='Directory of the python scripts cmig_analyze_bids_func '
                             'calls (default: $CMIG_PIPELINE_DIR/python)')
    parser.add_argument('--gradunwarp_flag', type=int, default=2, choices=[0, 1, 2],
                        help='Gradient unwarp: 0=skip, 1=apply if gradwarp type known, '
                             '2=fail on unknown gradwarp type (default: 2)')
    parser.add_argument('--cleanup_flag', type=int, default=1, choices=[0, 1],
                        help='Cleanup intermediate files (default: 1)')
    parser.add_argument('--std_orient', default=None,
                        help='cmig_proc_bids_func output orientation, e.g. LAS or '
                             'LPS; none keeps the input\'s (default: the module\'s, LAS)')
    parser.add_argument('--T2_flag', type=int, default=0, choices=[0, 1, 2],
                        help='T2w handling for FreeSurfer: 0=ignore, 1=use if available, '
                             '2=require (default: 0)')
    parser.add_argument('--fsurf_options', default='-all',
                        help='recon-all options (default: -all)')
    parser.add_argument('--FS_version', type=int,
                        default=cmig_fsurf_bids.DEFAULT_FS_VERSION,
                        help='FreeSurfer version number '
                             '(default: %d)' % cmig_fsurf_bids.DEFAULT_FS_VERSION)
    parser.add_argument('--setup_script', default=None,
                        help='Path to SetUpFreeSurfer.csh (default: $PUBSW/bin/SetUpFreeSurfer.csh)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='Display status messages, in the steps that take it '
                             '(default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing output, in every step (default: 0)')
    parser.add_argument('--compiled_flag', type=int, default=1, choices=[0, 1],
                        help='0 = MATLAB code, 1 = compiled app (default: 1)')
    parser.add_argument('--reproducible_flag', type=int, default=1, choices=[0, 1],
                        help=cmig_run_bids_module.REPRODUCIBLE_HELP)
    parser.add_argument('--appver', default=DEFAULT_APPVER,
                        help='Container image version tag (default: %(default)s)')
    parser.add_argument('--appdir', default=DEFAULT_APPDIR,
                        help='Directory containing the compiled app (default: %(default)s)')
    parser.add_argument('--pubsdir', default=DEFAULT_PUBSDIR,
                        help='Public software directory bound to /usr/pubsw (default: %(default)s)')
    parser.add_argument('--mrtdir', default=DEFAULT_MRTDIR,
                        help='MATLAB runtime directory on host (default: %(default)s)')
    parser.add_argument('--mrtdir_target', default=DEFAULT_MRTDIR_TARGET,
                        help='MATLAB runtime mount point inside container (default: %(default)s)')
    return parser


def parse_args():
    return build_parser().parse_args()


def _join_errmsgs(errmsgs):
    """Join accumulated messages into one line, or None if there are none."""
    errmsgs = [m for m in errmsgs if m]
    return '; '.join(errmsgs) if errmsgs else None


def run(args):
    infra = dict(
        compiled_flag=args.compiled_flag,
        reproducible_flag=args.reproducible_flag,
        appver=args.appver,
        appdir=args.appdir,
        pubsdir=args.pubsdir,
        mrtdir=args.mrtdir,
        mrtdir_target=args.mrtdir_target,
    )
    pid, sid = args.participant_id, args.session_id
    errmsgs = []     # failures that count against errcode
    notes = []       # messages that travel with errcode 0

    # the func steps run for a session with task fMRI data, unless turned off
    tasks = []
    if args.func_flag:
        tasks = args.task
        if tasks is None:
            tasks = cmig_get_bids_data.get_bids_tasks(pid, sid, args.raw_rootdir)
            print('tasks with usable data in rawdata: {0}'.format(
                ' '.join(tasks) if tasks else 'none'))
    run_func = bool(args.func_flag)

    # behavioral metrics first: they need only rawdata; their tasks come from
    #   the events files (rec, in beh/, included), not from the bold files
    if run_func:
        _demarcate('cmig_analyze_bids_behav')
        behav_tasks = sorted(cmig_select_bids_func_runs.events_tasks(
            args.raw_rootdir, pid, sid))
        if behav_tasks:
            errcode, errmsg = cmig_analyze_bids_behav.analyze_bids_behav(
                pid, sid, args.raw_rootdir, args.analysis_rootdir,
                tasks=behav_tasks,
                ignore_delay_flags=[int(t in args.ignore_delay_tasks)
                                    for t in behav_tasks],
                exclude_noncompliant_flag=bool(args.exclude_noncompliant_flag),
                forceflag=bool(args.forceflag), verbose=bool(args.verbose))
            if errcode != 0:
                print(f'cmig_pipeline_bids: ERROR: cmig_analyze_bids_behav '
                      f'failed: {errmsg}')
                errmsgs.append(f'cmig_analyze_bids_behav: {errmsg}')
        else:
            print('no task events in rawdata func/ or beh/')

    _demarcate('cmig_proc_bids_anat')
    errcode, errmsg = cmig_run_bids_module.run_module(
        'cmig_proc_bids_anat',
        participant_id=pid, session_id=sid,
        raw_rootdir=args.raw_rootdir, output_rootdir=args.proc_rootdir,
        gradunwarp_flag=args.gradunwarp_flag, cleanup_flag=args.cleanup_flag,
        forceflag=args.forceflag,
        **infra,
    )
    anat_ok = errcode == 0
    if not anat_ok:
        errmsgs.append(f'cmig_proc_bids_anat: {errmsg}')

    # without a processed T1w, proc_func processes without registration to T1
    func_ok = run_func and bool(tasks)
    if func_ok:
        _demarcate('cmig_proc_bids_func')
        errcode, errmsg = cmig_run_bids_module.run_module(
            'cmig_proc_bids_func',
            participant_id=pid, session_id=sid,
            raw_rootdir=args.raw_rootdir, output_rootdir=args.proc_rootdir,
            gradunwarp_flag=args.gradunwarp_flag, cleanup_flag=args.cleanup_flag,
            std_orient=args.std_orient, forceflag=args.forceflag,
            **infra,
        )
        if errcode != 0:
            errmsgs.append(f'cmig_proc_bids_func: {errmsg}')
            func_ok = False

    # the dwi steps run for a session with dwi data, unless turned off; the
    #   dwi processing registers to the processed T1w
    dwi_ok = bool(args.dwi_flag) and anat_ok and has_dwi(pid, sid, args.raw_rootdir)
    if args.dwi_flag and not dwi_ok:
        print('cmig_pipeline_bids: {0}; dwi steps skipped'.format(
            'no processed T1w' if not anat_ok else 'no usable dwi data in rawdata'))
    if dwi_ok:
        _demarcate('cmig_proc_bids_dwi')
        errcode, errmsg = cmig_proc_bids_dwi.proc_dwi(
            pid, sid, args.raw_rootdir, args.proc_rootdir,
            n_cpus=args.dwi_n_cpus, cleanup_flag=args.cleanup_flag,
            forceflag=args.forceflag)
        if errcode != 0:
            errmsgs.append(f'cmig_proc_bids_dwi: {errmsg}')
            dwi_ok = False

    fsurf_ok = False
    if anat_ok:
        _demarcate('cmig_fsurf_bids')
        errcode, errmsg = cmig_fsurf_bids.run(Namespace(
            participant_id=pid, session_id=sid,
            proc_rootdir=args.proc_rootdir, output_rootdir=args.fsurf_rootdir,
            T2_flag=args.T2_flag, options=args.fsurf_options,
            forceflag=args.forceflag,
            FS_version=args.FS_version, setup_script=args.setup_script,
        ))
        fsurf_ok = errcode == 0
        if not fsurf_ok:
            errmsgs.append(f'cmig_fsurf_bids: {errmsg}')

    # anat analysis and QC; each is independent of the others
    if fsurf_ok and args.analyze_anat_flag:
        _demarcate('cmig_analyze_bids_anat')
        errcode, errmsg = cmig_run_bids_module.run_module(
            'cmig_analyze_bids_anat',
            participant_id=pid, session_id=sid,
            proc_rootdir=args.proc_rootdir, fsurf_rootdir=args.fsurf_rootdir,
            output_rootdir=args.analysis_rootdir, cleanup_flag=args.cleanup_flag,
            forceflag=args.forceflag,
            **infra,
        )
        if errcode != 0:
            errmsgs.append(f'cmig_analyze_bids_anat: {errmsg}')

    postqc = []
    if anat_ok:
        postqc.append(('cmig_postqc_bids_anat', {}))
    if fsurf_ok:
        postqc.append(('cmig_postqc_bids_fsurf', {}))

    if func_ok and fsurf_ok:
        # tasks are independent of each other: a failed task does not count
        #   against errcode while another task succeeded
        task_errmsgs = []
        n_task_ok = 0
        for task in tasks:
            _demarcate('cmig_analyze_bids_func (task-{0})'.format(task))
            # options not given take the module's own defaults
            errcode, errmsg = cmig_run_bids_module.run_module(
                'cmig_analyze_bids_func',
                participant_id=pid, session_id=sid,
                proc_rootdir=args.proc_rootdir, fsurf_rootdir=args.fsurf_rootdir,
                output_rootdir=args.analysis_rootdir,
                task=task,
                cleanup_flag=args.cleanup_flag,
                prep_tseries_flag=1,
                trim_minutes=[5],
                rest_trim_minutes=[10],
                ignore_delay_flag=int(task in args.ignore_delay_tasks),
                exclude_noncompliant_flag=args.exclude_noncompliant_flag,
                pydir=args.pydir, forceflag=args.forceflag,
                **infra,
            )
            if errcode != 0:
                task_errmsgs.append(f'cmig_analyze_bids_func (task-{task}): {errmsg}')
            else:
                n_task_ok += 1
        if n_task_ok == 0:
            errmsgs.extend(task_errmsgs or ['no tasks analyzed'])
        else:
            notes.extend(task_errmsgs)

    # analyze_dwi samples the surfaces, so it needs a complete recon
    if dwi_ok and fsurf_ok:
        _demarcate('cmig_analyze_bids_dwi')
        errcode, errmsg = cmig_run_bids_module.run_module(
            'cmig_analyze_bids_dwi',
            participant_id=pid, session_id=sid,
            proc_rootdir=args.proc_rootdir, fsurf_rootdir=args.fsurf_rootdir,
            output_rootdir=args.analysis_rootdir, cleanup_flag=args.cleanup_flag,
            verbose=args.verbose, forceflag=args.forceflag,
            **infra,
        )
        if errcode != 0:
            errmsgs.append(f'cmig_analyze_bids_dwi: {errmsg}')

    # analyze is not an input to postqc; fsurf is optional to postqc_func,
    #   postqc_dwi and postqc_fmap, which need proc_func or proc_dwi
    if func_ok:
        postqc.append(('cmig_postqc_bids_func', {}))
        # it still runs, and fails, when usable field maps were not processed
        if skip_postqc_fmap(pid, sid, args.proc_rootdir, args.raw_rootdir):
            print('cmig_pipeline_bids: no processed and no usable acq-func '
                  'field map; cmig_postqc_bids_fmap skipped')
        else:
            postqc.append(('cmig_postqc_bids_fmap', {'acq_types': ['func']}))
    if dwi_ok:
        postqc.append(('cmig_postqc_bids_dwi', {}))
        if skip_postqc_fmap(pid, sid, args.proc_rootdir, args.raw_rootdir,
                            acq_type='dwi'):
            print('cmig_pipeline_bids: no processed and no usable acq-dwi '
                  'field map; cmig_postqc_bids_fmap (acq-dwi) skipped')
        else:
            postqc.append(('cmig_postqc_bids_fmap', {'acq_types': ['dwi']}))
    for name, extra in postqc:
        _demarcate(name)
        errcode, errmsg = cmig_run_bids_module.run_module(
            name,
            participant_id=pid, session_id=sid,
            proc_rootdir=args.proc_rootdir, fsurf_rootdir=args.fsurf_rootdir,
            output_rootdir=args.postqc_rootdir, cleanup_flag=args.cleanup_flag,
            verbose=args.verbose, forceflag=args.forceflag,
            **dict(extra, **infra),
        )
        if errcode != 0:
            errmsgs.append(f'{name}: {errmsg}')

    return (1 if errmsgs else 0), _join_errmsgs(errmsgs + notes)


def main():
    errcode, errmsg = run(parse_args())
    if errmsg:
        print(errmsg, file=sys.stderr)
    # the closing lines cmig_check_batch_errors.py judges a pipeline job by
    print(f'cmig_pipeline_bids: errcode = {errcode}')
    print('cmig_pipeline_bids: finished')
    sys.exit(errcode)


if __name__ == '__main__':
    main()
