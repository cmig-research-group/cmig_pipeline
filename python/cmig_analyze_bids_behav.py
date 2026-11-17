#!/usr/bin/env python3
"""
cmig_analyze_bids_behav.py

Behavioral metrics of one session's E-Prime tasks, from BIDS rawdata events.

Unlike cmig_analyze_bids_{anat,dwi,func}.py, this wraps no MATLAB app: it
selects runs with select_bids_func_runs() (cmig_select_bids_func_runs.py),
the same selection cmig_analyze_bids_func uses, and computes the metrics with
events_to_behavioral() (cmig_abcd_events_to_behavioral.py).

For each task, separately for the runs in rawdata func/ (scanned) and beh/
(no scan), never combined, writes
  {output_rootdir}/{pid}/{sid}/{func,beh}/{pid}_{sid}_task-{task}_behavioral.tsv
  {output_rootdir}/{pid}/{sid}/{func,beh}/{pid}_{sid}_task-{task}_desc-behav_status.tsv
in each directory where the task has events, and in func/ the run selection,
  {output_rootdir}/{pid}/{sid}/func/{pid}_{sid}_task-{task}_desc-select_runs.tsv
which cmig_analyze_bids_func reads. REC has no runs; its single events file is
in beh/. Old-style (release 7.0) events have none of the
response columns the metrics need, so a directory whose selected runs are
old-style is skipped with a WARNING.

Usage
-----
  cmig_analyze_bids_behav.py --participant_id sub-XXXXXXXX --session_id ses-XXA \\
      --raw_rootdir /path/to/rawdata --output_rootdir /path/to/output \\
      [--task mid sst nback rec] [--ignore_delay_flag 0|1 [0|1 ...]] \\
      [--overshoot_thresh 5] [--max_delay 12.5]
      [--exclude_noncompliant_flag 0|1] [--forceflag 0|1] [--verbose 0|1]

--ignore_delay_flag takes one value for every task, or one per --task.
--exclude_noncompliant_flag 1 makes the run selection count runs with
compliant 0 in the rawdata _scans.tsv as invalid.
Exit status 0 when no task failed.

Created:  10/02/2026 by Don Hagler
Prev mod: 10/02/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Last mod: 11/11/2026 by Don Hagler
"""
from __future__ import print_function
import argparse
import datetime
import os
import sys

import cmig_sing_utils
import cmig_write_status as cmig_status
from cmig_abcd_events_to_behavioral import events_to_behavioral
from cmig_select_bids_func_runs import (select_bids_func_runs, is_old_style,
                                        read_events, events_tasks,
                                        behav_status_stale, EVENTS_TASKS)

APPNAME = 'cmig_analyze_bids_behav'
TASKS = list(EVENTS_TASKS)


def build_parser():
    parser = argparse.ArgumentParser(
        description='Behavioral metrics of one session\'s E-Prime tasks.')
    parser.add_argument('--participant_id', required=True)
    parser.add_argument('--session_id', required=True)
    parser.add_argument('--raw_rootdir', required=True,
                        help='BIDS rawdata root holding the events and bold')
    parser.add_argument('--output_rootdir', required=True,
                        help='root of the output tree (created if absent)')
    parser.add_argument('--task', nargs='+', default=TASKS,
                        help='BIDS task labels (default: %(default)s)')
    parser.add_argument('--ignore_delay_flag', nargs='+', type=int,
                        default=[0], choices=[0, 1],
                        help='REDCap eprime_mismatch_ok flag: one value, or '
                             'one per --task')
    parser.add_argument('--overshoot_thresh', type=float, default=5.0)
    parser.add_argument('--max_delay', type=float, default=12.5)
    parser.add_argument('--exclude_noncompliant_flag', type=int, default=0,
                        choices=[0, 1],
                        help='runs with compliant 0 in the rawdata '
                             '_scans.tsv are invalid')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1])
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1])
    return parser


def parse_args():
    return build_parser().parse_args()


def _now():
    now = datetime.datetime.now()
    return now.strftime('%Y-%m-%dT%H:%M:%S.') + '{0:03d}'.format(
        now.microsecond // 1000)


def status_module(task):
    return 'task-{}_desc-behav'.format(task)


def prior_success(outdir, pid, sid, task):
    """True if a status file for this task in outdir records errcode 0."""
    status = cmig_status.read_status(
        cmig_status.status_fpath(outdir, pid, sid, status_module(task)))
    try:
        return status is not None and int(status.get('errcode', 1)) == 0
    except ValueError:
        return False


def write_status(outdir, pid, sid, task, errcode, errmsg, startdate,
                 ignore_delay_flag, exclude_noncompliant_flag=0):
    # the flags used are recorded, so a changed flag makes it stale
    fpath = cmig_status.cmig_write_status(
        pid, sid, outdir, status_module(task), errcode, errmsg,
        extra={'ignore_delay_flag': int(ignore_delay_flag),
               'exclude_noncompliant_flag': int(exclude_noncompliant_flag)},
        startdate=startdate)
    if errcode == 0:
        print('{0}: finished'.format(APPNAME))
        print('{0}: completed successfully'.format(APPNAME))
    # pure python with no compiled app; software_version identifies the code
    cmig_sing_utils.add_status_versions(fpath, app_name=APPNAME)


def remove_output(fname, verbose):
    if os.path.exists(fname):
        os.remove(fname)
        if verbose:
            print('{}: removed the earlier {}'.format(APPNAME, fname))


def behav_dir(pid, sid, task, modality, raw_rootdir, outdir, ignore_delay_flag,
              overshoot_thresh, max_delay, forceflag, verbose,
              exclude_noncompliant_flag=False):
    """Metrics for one task's runs in one rawdata directory: (errcode, errmsg).

    Called once the status gate has let the directory through, so the output
    is always rewritten; a path that writes none removes an earlier one, so a
    stale _behavioral.tsv is never left beside a new status.
    """
    stem = '{}_{}_task-{}'.format(pid, sid, task)
    fname_out = os.path.join(outdir, stem + '_behavioral.tsv')
    if task == 'rec':
        fnames = [os.path.join(raw_rootdir, pid, sid, 'beh', stem + '_events.tsv')]
        run_slots = None
    else:
        # func/: the selection is written beside the behavioral output as
        #   _desc-select_runs.tsv, which cmig_analyze_bids_func then reads
        runs, errcode, errmsg = select_bids_func_runs(
            pid, sid, task, raw_rootdir, modality=modality,
            outdir=outdir if modality == 'func' else None,
            overshoot_thresh=overshoot_thresh, max_delay=max_delay,
            ignore_delay_flag=ignore_delay_flag,
            exclude_noncompliant_flag=exclude_noncompliant_flag,
            verbose=verbose)
        if errcode:
            return errcode, errmsg
        selected = sorted((r for r in runs if r['selected']),
                          key=lambda r: r['slot'])
        if not selected:
            print('{}: WARNING: {} {} {} {}/: no valid runs, no behavioral '
                  'output'.format(APPNAME, pid, sid, task, modality))
            remove_output(fname_out, verbose)
            return 0, ''
        fnames = [r['fname_events'] for r in selected]
        run_slots = [r['run'] for r in selected]

    old = [f for f in fnames if is_old_style(read_events(f)[0])]
    if old:
        print('{}: WARNING: {} {} {} {}/: old-style events, which lack the '
              'response columns the metrics need; no behavioral output: {}'.format(
                  APPNAME, pid, sid, task, modality,
                  ', '.join(os.path.basename(f) for f in old)))
        remove_output(fname_out, verbose)
        return 0, ''

    behav, errcode, errmsg = events_to_behavioral(
        fnames, outdir=outdir, outstem=stem, run_slots=run_slots,
        forceflag=True, verbose=verbose)
    if errcode:
        remove_output(fname_out, verbose)
    return errcode, errmsg


def analyze_bids_behav(participant_id, session_id, raw_rootdir, output_rootdir,
                       tasks=None, ignore_delay_flags=None, overshoot_thresh=5.0,
                       max_delay=12.5, exclude_noncompliant_flag=False,
                       forceflag=False, verbose=True):
    """Behavioral metrics for each task; returns (errcode, errmsg)."""
    pid, sid = participant_id, session_id
    tasks = [t.lower() for t in (tasks or TASKS)]
    flags = list(ignore_delay_flags or [0])
    if len(flags) == 1:
        flags = flags * len(tasks)
    if len(flags) != len(tasks):
        return 1, 'ignore_delay_flag: give one value or one per task'
    sesdir = os.path.join(raw_rootdir, pid, sid)
    if not os.path.isdir(sesdir):
        return 1, 'rawdata session not found: {}'.format(sesdir)

    errmsgs = []
    found_mods = events_tasks(raw_rootdir, pid, sid)
    for task, flag in zip(tasks, flags):
        stem = '{}_{}_task-{}'.format(pid, sid, task)
        if task not in TASKS:
            print('{}: WARNING: {} {} {}: not an E-Prime task with events '
                  '({}), skipped'.format(APPNAME, pid, sid, task,
                                         ' '.join(TASKS)))
            continue
        found = task in found_mods
        for modality in found_mods.get(task, []):
            outdir = os.path.join(output_rootdir, pid, sid, modality)
            if not forceflag and prior_success(outdir, pid, sid, task):
                why = behav_status_stale(
                    cmig_status.status_fpath(outdir, pid, sid, status_module(task)),
                    raw_rootdir, pid, sid, task, modality, flag,
                    int(exclude_noncompliant_flag))
                if not why:
                    if verbose:
                        print('{}: {} {} {} {}/: done before, skipping '
                              '(forceflag=0)'.format(APPNAME, pid, sid, task, modality))
                    continue
                print('{}: {} {} {} {}/: rerunning: {}'.format(
                    APPNAME, pid, sid, task, modality, why))
            startdate = _now()
            os.makedirs(outdir, exist_ok=True)
            try:
                errcode, errmsg = behav_dir(
                    pid, sid, task, modality, raw_rootdir, outdir, bool(flag),
                    overshoot_thresh, max_delay, forceflag, verbose,
                    exclude_noncompliant_flag=bool(exclude_noncompliant_flag))
            except Exception as exc:
                errcode, errmsg = 1, '{}: {}'.format(type(exc).__name__, exc)
                remove_output(os.path.join(outdir, stem + '_behavioral.tsv'),
                              verbose)
            if errcode:
                errmsgs.append('{} {}/: {}'.format(task, modality, errmsg))
                print('{}: ERROR: {} {} {} {}/: {}'.format(
                    APPNAME, pid, sid, task, modality, errmsg), file=sys.stderr)
            write_status(outdir, pid, sid, task, errcode, errmsg, startdate,
                         flag, exclude_noncompliant_flag)
        if not found and verbose:
            print('{}: NOTE: {} {} {}: no events in func/ or beh/'.format(
                APPNAME, pid, sid, task))
    if errmsgs:
        return 1, '; '.join(errmsgs)
    return 0, ''


def run(args):
    return analyze_bids_behav(
        args.participant_id, args.session_id, args.raw_rootdir,
        args.output_rootdir, tasks=args.task,
        ignore_delay_flags=args.ignore_delay_flag,
        overshoot_thresh=args.overshoot_thresh, max_delay=args.max_delay,
        exclude_noncompliant_flag=bool(args.exclude_noncompliant_flag),
        forceflag=bool(args.forceflag), verbose=bool(args.verbose))


def main():
    errcode, errmsg = run(parse_args())
    if errmsg:
        print('{}: ERROR: {}'.format(APPNAME, errmsg), file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
