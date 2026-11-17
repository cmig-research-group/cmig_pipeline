#!/usr/bin/env python
"""
cmig_select_bids_func_runs.py

Select the runs of one session's task-fMRI task to analyze, from BIDS rawdata.

Public function
---------------
select_bids_func_runs(participant_id, session_id, task, raw_rootdir,
                      outdir=None, modality='func', overshoot_thresh=5.0,
                      max_delay=12.5, ignore_delay_flag=False,
                      exclude_noncompliant_flag=False, forceflag=False,
                      verbose=True)
Returns (runs, errcode, errmsg).

  runs is a list of dicts, one per run in BIDS run order, with keys run (int),
  filename (the bold file, relative to the session directory, or ''),
  fname_events (full path, or ''), selected (0/1), slot (1, 2 or None) and
  notes (why the run was excluded, or the checks skipped for it).

  The first 1 or 2 valid runs are selected, to slots 1 and 2 (the release
  run slots r01 and r02). A run is invalid when its bold is not usable (usable
  0 in the rawdata _scans.tsv, or not listed there, as cmig_get_bids_data
  takes it; func only), it has no events file, its
  events file is desc-invalid (a trial count other than the task's, or, from
  the converter of 10/05/2026, no responses), it has no response recorded (no response_time in any task event), its
  timing_ok is 0, its events end more than overshoot_thresh seconds past the
  end of its scan, or its E-Prime start_time is more than max_delay minutes
  from its bold's AcquisitionDateTime (the bold json; the first stored
  volume, so the same on every vendor, unlike SeriesTime, whose offset from
  it runs to ~2.5 min either way). The scan end is volumes x TR,
  volumes from the bold NIfTI header, onsets being from the first stored
  volume as BIDS requires. ignore_delay_flag (the REDCap eprime_mismatch_ok
  flag) skips the start_time check. With exclude_noncompliant_flag, a run
  whose bold has compliant 0 in the rawdata _scans.tsv is invalid too.

  When no run is valid, the first 1 or 2 runs that pass the events checks
  (every check above but usable and compliant) are selected instead, for the behavioral
  metrics only, with the note "behavioral metrics only: no valid run".

  modality 'func' takes the runs of the session's bold files; 'beh' the runs
  of its beh/ events files, which have no scan, so only the events-file
  criteria apply.

  raw_rootdir may also be the cmig_proc_bids_func tree, whose bold files are
  _desc-preproc_bold and whose events files are copies with the rawdata
  names; filename then still gives the rawdata name of the bold.

  Old-style events files, with none of the per-run timing columns, are
  checked without timing_ok and start_time, with a WARNING.

  With outdir, also writes {pid}_{sid}_task-{task}_desc-select_runs.tsv
  there (func only), on every call (forceflag does not apply to it): one row per bold run, with filename, modality, task and
  run, so cmig_compile_bids_tabdat.py compiles it into mr_y_qc__post__series.
  An empty cell (the slot of a run not selected, no notes) is written n/a.

Usage
-----
  cmig_select_bids_func_runs.py --participant_id sub-XXXXXXXX \\
      --session_id ses-XXA --task mid --raw_rootdir /path/to/rawdata \\
      [--outdir /path/to/output] [--modality func|beh] \\
      [--overshoot_thresh 5] [--max_delay 12.5] [--ignore_delay_flag 0|1] \\
      [--exclude_noncompliant_flag 0|1] [--forceflag 0|1] [--verbose 0|1]

Created:  10/02/2026 by Don Hagler
Prev mod: 10/02/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 11/10/2026 by Don Hagler
Last mod: 11/11/2026 by Don Hagler
"""
from __future__ import print_function
import argparse
import csv
import datetime
import glob
import gzip
import json
import os
import re
import struct
import sys

SCRIPTNAME = 'cmig_select_bids_func_runs'

# per-run timing columns written by cmig_abcd_eprime_to_events.py; an events
# file with none of them is old-style (release 7.0)
TIMING_COLS = ('start_delay', 'start_delay_diff', 'init_delay_diff',
               'timing_ok', 'start_time')
SPECIAL_TYPES = ('dummy', 'scanner_trigger', 'task_init')
# SST trial types that record no response: a go trial unanswered, a stop trial
#   correctly withheld. SST writes response_time 0 for them, but 0 is also the
#   response_time of a late response (correctlate_go, incorrectlate_go), so
#   for SST a response is read from the trial type
SST_NORESP_TYPES = ('noresp_go', 'correct_stop')
NMAX_RUNS = 2
# the notes of a run that passes the events checks but is not used for imaging
UNUSABLE_NOTES = ('not usable', 'not compliant')
BEHAV_ONLY_NOTE = 'behavioral metrics only: no valid run'
# the proc bold (cmig_proc_bids_func) carries this before its suffix
PROC_DESC = '_desc-preproc'
# the E-Prime tasks with events files (REC has one, in beh/, with no run)
EVENTS_TASKS = ('mid', 'sst', 'nback', 'rec')
COL_PREFIX = 'mr_y_qc__post__series__tfmri__'


def parse_args():
    parser = argparse.ArgumentParser(
        description='Select the task-fMRI runs of one session to analyze.')
    parser.add_argument('--participant_id', required=True)
    parser.add_argument('--session_id', required=True)
    parser.add_argument('--task', required=True,
                        help='BIDS task label, e.g. mid, sst, nback')
    parser.add_argument('--raw_rootdir', required=True,
                        help='BIDS rawdata root holding the bold and events')
    parser.add_argument('--outdir', default=None,
                        help='directory for the _desc-select_runs.tsv')
    parser.add_argument('--modality', default='func', choices=['func', 'beh'])
    parser.add_argument('--overshoot_thresh', type=float, default=5.0,
                        help='seconds events may end past the scan end')
    parser.add_argument('--max_delay', type=float, default=12.5,
                        help='minutes allowed between E-Prime start and scan')
    parser.add_argument('--ignore_delay_flag', type=int, default=0,
                        choices=[0, 1])
    parser.add_argument('--exclude_noncompliant_flag', type=int, default=0,
                        choices=[0, 1],
                        help='runs with compliant 0 in the rawdata _scans.tsv '
                             'are invalid')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1])
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1])
    return parser.parse_args()


def events_tasks(raw_rootdir, participant_id, session_id, tasks=None):
    """{task: [modality, ...]}: for each of EVENTS_TASKS (or of tasks among
    them) with events files, the rawdata directories (func, beh) holding them.
    The one definition of a session's tasks with events, for the behavioral
    step, cmig_pipeline_bids.py and mkjobs."""
    sesdir = os.path.join(raw_rootdir, participant_id, session_id)
    found = {}
    for task in EVENTS_TASKS:
        if tasks is not None and task not in tasks:
            continue
        stem = '{}_{}_task-{}'.format(participant_id, session_id, task)
        pattern = stem + ('_events.tsv' if task == 'rec' else '_run-*_events.tsv')
        mods = [mod for mod in ('func', 'beh')
                if glob.glob(os.path.join(sesdir, mod, pattern))]
        if mods:
            found[task] = mods
    return found


def events_mtime(raw_rootdir, participant_id, session_id, task=None,
                 modality=None):
    """The newest modification time (datetime) of a session's events files,
    of one task and rawdata directory when given; None if there are none."""
    sesdir = os.path.join(raw_rootdir, participant_id, session_id)
    pattern = '{}_{}_task-{}*_events.tsv'.format(participant_id, session_id,
                                                 task or '*')
    fnames = glob.glob(os.path.join(sesdir, modality or '*', pattern))
    if not fnames:
        return None
    return datetime.datetime.fromtimestamp(max(os.path.getmtime(f)
                                               for f in fnames))


def behav_status_stale(fname_status, raw_rootdir, participant_id, session_id,
                       task, modality, ignore_delay_flag=None,
                       exclude_noncompliant_flag=None):
    """Why a behavioral step's status is out of date, or '' when it is not:
    the events files changed after the run started, or (func/) the REDCap
    ignore_delay_flag or the exclude_noncompliant_flag differs from the one
    it recorded, or (func/, exclude_noncompliant_flag) a bold's compliant
    changed (compliant_changed). A status written before ignore_delay_flag
    was recorded is not judged on it; one written before
    exclude_noncompliant_flag was recorded is taken as 0."""
    import cmig_write_status
    status = cmig_write_status.read_status(fname_status)
    if status is None:
        return ''
    started = cmig_write_status.status_time(fname_status, 'startdate')
    newest = events_mtime(raw_rootdir, participant_id, session_id, task, modality)
    if started is not None and newest is not None and newest > started:
        return 'events changed after the behavioral step ran'
    recorded = (status.get('ignore_delay_flag') or '').strip()
    if (modality == 'func' and ignore_delay_flag is not None and
            recorded in ('0', '1') and int(recorded) != int(ignore_delay_flag)):
        return 'ignore_delay_flag changed from {} to {}'.format(
            recorded, int(ignore_delay_flag))
    recorded = (status.get('exclude_noncompliant_flag') or '0').strip()
    if (modality == 'func' and exclude_noncompliant_flag is not None and
            recorded in ('0', '1') and
            int(recorded) != int(exclude_noncompliant_flag)):
        return 'exclude_noncompliant_flag changed from {} to {}'.format(
            recorded, int(exclude_noncompliant_flag))
    if modality == 'func' and exclude_noncompliant_flag and started is not None:
        why = compliant_changed(os.path.dirname(fname_status), raw_rootdir,
                                participant_id, session_id, task, started)
        if why:
            return why
    return ''


def compliant_changed(outdir, raw_rootdir, participant_id, session_id, task,
                      started):
    """Why the compliant 0 bolds of a task in the rawdata _scans.tsv differ
    from the runs its _desc-select_runs.tsv in outdir notes as not compliant,
    or ''; compared only where the _scans.tsv changed after started."""
    sesdir = os.path.join(raw_rootdir, participant_id, session_id)
    fscans = os.path.join(sesdir, '{}_{}_scans.tsv'.format(participant_id,
                                                            session_id))
    stem = '{}_{}_task-{}'.format(participant_id, session_id, task)
    fsel = os.path.join(outdir, stem + '_desc-select_runs.tsv')
    if not os.path.isfile(fscans) or not os.path.isfile(fsel):
        return ''
    if datetime.datetime.fromtimestamp(os.path.getmtime(fscans)) <= started:
        return ''
    compliant = read_scans_column(sesdir, participant_id, session_id,
                                  'compliant') or {}
    pattern = re.compile(r'^func/' + re.escape(stem) +
                         r'_run-\d+_bold\.nii(\.gz)?$')
    now = set(f for f, v in compliant.items() if v == 0 and pattern.match(f))
    with open(fsel) as fid:
        rows = list(csv.DictReader(fid, delimiter='\t'))
    notes = COL_PREFIX + 'select_notes'
    was = set(r.get('filename', '') for r in rows
              if 'not compliant' in (r.get(notes) or '').split('; '))
    if now == was:
        return ''
    parts = [label + ' ' + ' '.join(sorted(os.path.basename(f) for f in names))
             for label, names in (('now', now - was), ('was', was - now)) if names]
    return 'compliant changed in scans.tsv: not compliant ' + '; '.join(parts)


def _warn(verbose, msg):
    if verbose:
        print('{}: WARNING: {}'.format(SCRIPTNAME, msg))


def is_old_style(header):
    """True for an events file with none of the per-run timing columns."""
    return not any(c in header for c in TIMING_COLS)


def read_events(fname):
    """Header and rows (dicts) of an events TSV, quotes stripped."""
    with open(fname) as fid:
        rows = list(csv.reader(fid, delimiter='\t'))
    if not rows:
        return [], []
    header = [h.strip().strip('"') for h in rows[0]]
    body = [dict(zip(header, [v.strip().strip('"') for v in r]))
            for r in rows[1:] if r]
    return header, body


def _float(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def _finite(val):
    """True for a value that parses as a finite number."""
    try:
        v = float(val)
    except (TypeError, ValueError):
        return False
    return v == v and v not in (float('inf'), float('-inf'))


def has_responses(rows, task):
    """True if any task event of one run's events (dicts) records a response:
    for SST, a trial type other than noresp_go and correct_stop, or a finite
    response_time_go_adjusted; for other tasks, a finite response_time."""
    events = []
    for row in rows:
        tt = str(row.get('trial_type') or '').strip().strip('"')
        if tt and tt.lower() not in SPECIAL_TYPES + ('n/a', 'nan'):
            events.append((tt, row))
    if str(task).lower() == 'sst':
        return (any(tt not in SST_NORESP_TYPES for tt, _ in events) or
                any(_finite(row.get('response_time_go_adjusted')) for _, row in events))
    return any(_finite(row.get('response_time')) for _, row in events)


def _first_value(body, col):
    for row in body:
        v = row.get(col, '')
        if v and v.lower() != 'n/a':
            return v
    return None


def nifti_nvols(fname):
    """Number of volumes, dim[4], from a .nii or .nii.gz header."""
    opener = gzip.open if fname.endswith('.gz') else open
    with opener(fname, 'rb') as fid:
        hdr = fid.read(348)
    endian = '<' if struct.unpack('<i', hdr[:4])[0] == 348 else '>'
    dim = struct.unpack(endian + '8h', hdr[40:56])
    return dim[4] if dim[0] >= 4 else 1


def bold_acq_datetime(fname_bold):
    """AcquisitionDateTime of a bold's json sidecar, as a datetime, or None.

    The first stored volume, the BIDS onset 0; in the rawdata and proc trees
    alike, unlike _scans.tsv, which only rawdata has.
    """
    fname_json = re.sub(r'\.nii(\.gz)?$', '.json', fname_bold)
    if not os.path.isfile(fname_json):
        return None
    try:
        with open(fname_json) as fid:
            val = json.load(fid).get('AcquisitionDateTime')
    except (IOError, ValueError):
        return None
    return _parse_time(val)


def _parse_time(val):
    for fmt in ('%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f'):
        try:
            return datetime.datetime.strptime(val, fmt)
        except (TypeError, ValueError):
            continue
    return None


def read_scans_column(sesdir, participant_id, session_id, column,
                      verbose=False):
    """{filename: value (float or None)} of one column of the session's
    rawdata _scans.tsv, the filename relative to the session directory; None
    when there is no _scans.tsv or no such column."""
    fname = os.path.join(sesdir, '{}_{}_scans.tsv'.format(participant_id,
                                                           session_id))
    if not os.path.isfile(fname):
        _warn(verbose, 'no {}: {} not checked'.format(fname, column))
        return None
    with open(fname) as fid:
        rows = list(csv.DictReader(fid, delimiter='\t'))
    if not rows or column not in rows[0]:
        _warn(verbose, 'no {} column in {}: {} not checked'.format(
            column, fname, column))
        return None
    return dict((row.get('filename', ''), _float(row.get(column)))
                for row in rows)


def read_scans_usable(sesdir, participant_id, session_id, verbose=False):
    """{filename: usable (0/1)} from the session's rawdata _scans.tsv; None
    when there is no _scans.tsv or no usable column, all runs then counting
    as usable, as cmig_get_bids_data takes them."""
    vals = read_scans_column(sesdir, participant_id, session_id, 'usable',
                             verbose)
    if vals is None:
        return None
    return dict((k, int(v == 1)) for k, v in vals.items())


def run_usable(usable, sesdir, fname_bold):
    """True unless the bold is marked unusable in usable (from
    read_scans_usable()) or, usable having been read, is not listed there, as
    cmig_get_bids_data keeps only the scans listed with usable 1."""
    if usable is None or not fname_bold:
        return True
    return bool(usable.get(os.path.relpath(fname_bold, sesdir), 0))


def run_compliant(compliant, sesdir, fname_bold):
    """False only for a bold with compliant 0 in compliant (from
    read_scans_column()); -1 (not applicable), blank and unlisted pass."""
    if compliant is None or not fname_bold:
        return True
    return compliant.get(os.path.relpath(fname_bold, sesdir)) != 0


def find_runs(sesdir, stem, modality):
    """BIDS run number -> (bold file or '', events file or '')."""
    runs = {}
    datadir = os.path.join(sesdir, modality)
    if modality == 'func':
        # the proc tree's bold files when there are any, else the rawdata
        #   ones, so a rawdata-named file left in the proc tree is not read
        fnames = sorted(glob.glob(os.path.join(
            datadir, stem + '_run-*' + PROC_DESC + '_bold.nii*')))
        desc = PROC_DESC
        if not fnames:
            fnames = sorted(glob.glob(os.path.join(datadir, stem + '_run-*_bold.nii*')))
            desc = ''
        for f in fnames:
            m = re.search(r'_run-(\d+)' + re.escape(desc) + r'_bold\.nii(\.gz)?$', f)
            if m:
                # with both, the .nii.gz
                entry = runs.setdefault(int(m.group(1)), ['', ''])
                if not entry[0] or f.endswith('.gz'):
                    entry[0] = f
    for f in glob.glob(os.path.join(datadir, stem + '_run-*_events.tsv')):
        m = re.search(r'_run-(\d+)(?:_desc-invalid)?_events\.tsv$', f)
        if m:
            entry = runs.setdefault(int(m.group(1)), ['', ''])
            # a valid events file wins over a stale desc-invalid one
            if not entry[1] or 'desc-invalid' in entry[1]:
                entry[1] = f
    return runs


def check_run(run, fname_bold, fname_events, sesdir, modality,
              overshoot_thresh, max_delay, ignore_delay_flag, verbose, task,
              usable=True, compliant=True):
    """(valid, notes, old_style, events_ok) for one run: events_ok from the
    events checks alone, valid also requiring the run to be usable and
    compliant."""
    events_ok, notes, old_style = check_events(
        run, fname_bold, fname_events, sesdir, modality, overshoot_thresh,
        max_delay, ignore_delay_flag, verbose, task)
    if not compliant:
        notes = ['not compliant'] + notes
    if not usable:
        notes = ['not usable'] + notes
    return events_ok and usable and compliant, notes, old_style, events_ok


def check_events(run, fname_bold, fname_events, sesdir, modality,
                 overshoot_thresh, max_delay, ignore_delay_flag, verbose, task):
    """(valid, notes, old_style) for one run, from its events checks."""
    if modality == 'func' and not fname_bold:
        return False, ['no bold file'], False
    if not fname_events:
        return False, ['no events file'], False
    header, body = read_events(fname_events)
    if 'desc-invalid' in os.path.basename(fname_events):
        # the converter writes desc-invalid for a wrong trial count or no
        #   responses; the rows say which
        return False, ['desc-invalid ({})'.format(
            'trial count' if has_responses(body, task) else 'no responses')], False
    old_style = is_old_style(header)
    reasons = []
    skipped = []

    # a run with no responses at all (a participant asleep or not engaged);
    #   its response-defined conditions (e.g. SST correct_go/noresp_go)
    #   would misrepresent it
    if task.lower() != 'sst' and 'response_time' not in header:
        skipped.append('no response_time column: response check skipped')
    elif not has_responses(body, task):
        reasons.append('no responses')

    if old_style:
        skipped.append('old-style events: timing_ok and eprime-scan checks '
                       'skipped')
    else:
        timing_ok = _first_value(body, 'timing_ok')
        if timing_ok == '0':
            reasons.append('timing_ok=0')
        elif timing_ok is None:
            # the behavioral version, or a Siemens-version file with no
            #   scanner info: start delays not measured
            skipped.append('timing_ok n/a: timing check skipped')

    if modality == 'func':
        # events past the scan end, in event time
        fname_json = re.sub(r'\.nii(\.gz)?$', '.json', fname_bold)
        TR = None
        if os.path.isfile(fname_json):
            try:
                with open(fname_json) as fid:
                    TR = _float(json.load(fid).get('RepetitionTime'))
            except (IOError, OSError, ValueError):
                return False, ['unreadable bold json'], old_style
        try:
            nvols = nifti_nvols(fname_bold)
        except (IOError, OSError, EOFError, struct.error):
            return False, ['unreadable bold header'], old_style
        events_end = None
        for row in body:
            onset = _float(row.get('onset'))
            if onset is None:
                continue
            if row.get('trial_type', '').lower() in SPECIAL_TYPES:
                continue
            end = onset + (_float(row.get('duration')) or 0.0)
            if events_end is None or end > events_end:
                events_end = end
        if TR is None:
            skipped.append('no RepetitionTime: events-past-scan-end check '
                           'skipped')
        elif events_end is not None:
            # onsets are from the first stored volume, dummy volumes included
            scan_end = nvols * TR
            overshoot = events_end - scan_end
            if overshoot > overshoot_thresh:
                reasons.append('events end {:.1f} s past the scan end'.format(
                    overshoot))

        # eprime-scan time difference
        if not old_style:
            if ignore_delay_flag:
                skipped.append('ignore_delay_flag: eprime-scan check skipped')
            else:
                start_time = _parse_time(_first_value(body, 'start_time'))
                acq_time = bold_acq_datetime(fname_bold)
                if start_time is None or acq_time is None:
                    skipped.append('no {}: eprime-scan check skipped'.format(
                        'start_time' if start_time is None
                        else 'AcquisitionDateTime'))
                else:
                    delay = (start_time - acq_time).total_seconds() / 60.0
                    if abs(delay) > max_delay:
                        reasons.append('eprime-scan delay {:.2f} min'.format(
                            delay))
                        _warn(verbose, 'run {}: eprime-scan delay {:.2f} min, '
                              'over {} min'.format(run, delay, max_delay))
    return not reasons, reasons + skipped, old_style


def write_select_tsv(fname, participant_id, session_id, task, runs):
    cols = ['participant_id', 'session_id', 'filename', 'modality', 'task',
            'run', COL_PREFIX + 'select_indicator', COL_PREFIX + 'slot_num',
            COL_PREFIX + 'select_notes']
    with open(fname, 'w') as fid:
        writer = csv.writer(fid, delimiter='\t', lineterminator='\n')
        writer.writerow(cols)
        for r in runs:
            if not r['filename']:
                continue
            # BIDS: n/a for an empty cell
            writer.writerow([participant_id, session_id, r['filename'], 'bold',
                             task, '{:02d}'.format(r['run']), r['selected'],
                             r['slot'] if r['slot'] else 'n/a',
                             '; '.join(r['notes']) or 'n/a'])


def select_bids_func_runs(participant_id, session_id, task, raw_rootdir,
                          outdir=None, modality='func', overshoot_thresh=5.0,
                          max_delay=12.5, ignore_delay_flag=False,
                          exclude_noncompliant_flag=False, forceflag=False,
                          verbose=True):
    """Select the runs of one session's task to analyze; see module doc."""
    task = task.lower()
    sesdir = os.path.join(raw_rootdir, participant_id, session_id)
    if not os.path.isdir(sesdir):
        return [], 1, 'rawdata session not found: {}'.format(sesdir)
    stem = '{}_{}_task-{}'.format(participant_id, session_id, task)
    found = find_runs(sesdir, stem, modality)

    usable = None
    compliant = None
    if modality == 'func':
        usable = read_scans_usable(sesdir, participant_id, session_id, verbose)
        if exclude_noncompliant_flag:
            compliant = read_scans_column(sesdir, participant_id, session_id,
                                          'compliant', verbose)

    runs = []
    old_runs = []
    for run in sorted(found):
        fname_bold, fname_events = found[run]
        valid, notes, old_style, events_ok = check_run(
            run, fname_bold, fname_events, sesdir, modality,
            overshoot_thresh, max_delay, ignore_delay_flag, verbose, task,
            usable=run_usable(usable, sesdir, fname_bold),
            compliant=run_compliant(compliant, sesdir, fname_bold))
        if old_style:
            old_runs.append(run)
        runs.append({
            'run': run,
            # the rawdata name, also for a proc bold: the series identity
            'filename': (os.path.relpath(fname_bold, sesdir).replace(
                PROC_DESC + '_bold.nii', '_bold.nii') if fname_bold else ''),
            'fname_events': fname_events,
            'selected': 0,
            'slot': None,
            'notes': notes,
            'valid': valid,
            'events_ok': events_ok,
        })

    # the first valid runs fill the slots; with no valid run, the runs that
    #   pass the events checks do, for the behavioral metrics only
    behav_only = not any(r['valid'] for r in runs)
    nselected = 0
    for r in runs:
        if not (r['events_ok'] if behav_only else r['valid']):
            continue
        if nselected < NMAX_RUNS:
            nselected += 1
            r['selected'] = 1
            r['slot'] = nselected
            if behav_only:
                n = len([x for x in r['notes'] if x in UNUSABLE_NOTES])
                r['notes'] = r['notes'][:n] + [BEHAV_ONLY_NOTE] + r['notes'][n:]
        else:
            r['notes'] = ['beyond the first {}'.format(NMAX_RUNS)] + r['notes']
    for r in runs:
        del r['valid'], r['events_ok']

    if old_runs:
        _warn(verbose, '{} {} {}: old-style events in run(s) {}: timing_ok and '
              'eprime-scan checks skipped'.format(
                  participant_id, session_id, task,
                  ', '.join(str(r) for r in old_runs)))
    if verbose:
        for r in runs:
            print('{}: {} {} {} run {}: {}{}'.format(
                SCRIPTNAME, participant_id, session_id, task, r['run'],
                'slot {}'.format(r['slot']) if r['selected'] else 'excluded',
                (' -- ' + '; '.join(r['notes'])) if r['notes'] else ''))

    if outdir and modality == 'func':
        os.makedirs(outdir, exist_ok=True)
        # written every time, so it always records the selection just made
        fname = os.path.join(outdir, stem + '_desc-select_runs.tsv')
        write_select_tsv(fname, participant_id, session_id, task, runs)
    return runs, 0, ''


def run(args):
    runs, errcode, errmsg = select_bids_func_runs(
        args.participant_id, args.session_id, args.task, args.raw_rootdir,
        outdir=args.outdir, modality=args.modality,
        overshoot_thresh=args.overshoot_thresh, max_delay=args.max_delay,
        ignore_delay_flag=bool(args.ignore_delay_flag),
        exclude_noncompliant_flag=bool(args.exclude_noncompliant_flag),
        forceflag=bool(args.forceflag), verbose=bool(args.verbose))
    return errcode, errmsg


def main():
    errcode, errmsg = run(parse_args())
    if errmsg:
        print('{}: ERROR: {}'.format(SCRIPTNAME, errmsg), file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
