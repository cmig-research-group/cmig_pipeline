#!/usr/bin/env python
"""
Convert one session's ABCD E-Prime file(s) to BIDS _events.tsv in rawdata.

Per-session driver for cmig_abcd_eprime_to_events.eprime_to_events() (F1).
Writes only into rawdata; the timing files (F2) and behavioral (F3) outputs,
timing QC included, go to derivatives and are produced elsewhere.

The session's rows are assembled here from the two standard ABCD metadata
files, given by --eprime_csv and --pcinfo_csv, so no prepared manifest is
needed:

  DAL_ABCD_QC_eprime.csv           (175 MB)
  DAL_ABCD_QC_combined_pcinfo.csv  (1.15 GB)

Both are cut down with grep before any CSV parsing -- a full parse of the pair
costs ~28 s against ~0.64 s for grep, which is what makes per-session assembly
practical. pcinfo.csv supplies Manufacturer and SoftwareVersion, which
eprime.csv lacks and without which skipTRs falls back to 0 and shifts every
onset. To run against a subset, cut both files down to the rows wanted and
point --eprime_csv / --pcinfo_csv at them; the column names are theirs.

A subject-session-task can span more than one raw E-Prime file; each is
processed in turn. An E-Prime run goes to func/ only when eprime.csv links it
to a scan by SeriesInstanceUID; its BIDS run number is then read from the
rawdata bold json with that SeriesInstanceUID, under --output_rootdir, so the
events file name matches the bold file name, and a linked scan with no bold
json fails the task. For testing into a directory with no bold files,
--require_bold 0 instead numbers the runs by ranking SeriesTime over all of the
task's series in pcinfo.csv. A run eprime.csv links to no scan in an otherwise
linked file is skipped. A file with no linked run, or any row
flagged eprime_behav_only, goes to beh/ with every run, numbered by its E-Prime
run position, continuing across a task's behavioral-only files; so does a file with no scanner info (Manufacturer and
SoftwareVersion from no bold json, task scan or pcinfo.csv row), whose skipTRs,
and so onset zero point, is unknown.

REC is not a taskname in eprime.csv. It is carried by the nBack rows, whose
own eprime_file_name is the WM file and whose eprime_extra_file_name is the
REC file; both are read under their own names rather than being renamed.

Usage:
  cmig_convert_abcd_eprime.py --participant_id sub-XXXXXXXX \\
      --session_id ses-XXA [--task MID SST nBack REC] \\
      --output_rootdir /path/to/rawdata [--require_bold 0]

Exit status is 0 when no requested task failed. A task with no rows is skipped
rather than failed -- about 5% of sessions lack one.

Created:  08/03/2026 by Don Hagler
Prev mod: 09/28/2026 by Don Hagler
Prev mod: 09/29/2026 by Don Hagler
Prev mod: 10/01/2026 by Don Hagler
Last mod: 10/05/2026 by Don Hagler
"""

from __future__ import print_function

import argparse
import csv
import glob
import json
import os
import re
import subprocess
import sys

import pandas as pd

from cmig_abcd_eprime_to_events import eprime_to_events

SCRIPTNAME = os.path.basename(__file__)

TASKS = ['MID', 'SST', 'nBack', 'REC']

# taskname in eprime.csv holding each task's rows. REC has none of its own.
SOURCE_TASK = {'MID': 'MID', 'SST': 'SST', 'nBack': 'nBack', 'REC': 'nBack'}

# EventName <-> BIDS ses- label. No 09A entry: ABCD does not acquire imaging at
# the 9-year event, so any E-Prime file labeled 9_year_follow_up_y_arm_1 is an
# upload labeling mistake, not real data.
EVENT_TO_SES = {
    'baseline_year_1_arm_1':    '00A',
    '2_year_follow_up_y_arm_1': '02A',
    '4_year_follow_up_y_arm_1': '04A',
    '6_year_follow_up_y_arm_1': '06A',
    '8_year_follow_up_y_arm_1': '08A',
}
SES_TO_EVENT = dict((v, k) for k, v in EVENT_TO_SES.items())


def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert one session\'s ABCD E-Prime file(s) to BIDS '
                    '_events.tsv in rawdata.')
    parser.add_argument('--participant_id', required=True,
                        help='Participant ID (e.g., sub-XXXXXXXX)')
    parser.add_argument('--session_id', required=True,
                        help='Session ID (e.g., ses-00A)')
    parser.add_argument('--task', nargs='+', default=TASKS, choices=TASKS,
                        help='Task label(s) (default: MID SST nBack REC)')
    parser.add_argument('--output_rootdir', required=True,
                        help='Root directory of BIDS rawdata (created if absent)')
    parser.add_argument('--eprime_csv', required=True,
                        help='ABCD E-Prime metadata CSV (DAL_ABCD_QC_eprime.csv)')
    parser.add_argument('--pcinfo_csv', required=True,
                        help='ABCD protocol compliance CSV, source of '
                             'Manufacturer and SoftwareVersion '
                             '(DAL_ABCD_QC_combined_pcinfo.csv)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1],
                        help='Display status messages (default: 1)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing output (default: 0)')
    parser.add_argument('--require_bold', type=int, default=1, choices=[0, 1],
                        help='Take run numbers from the bold jsons under '
                             '--output_rootdir and fail on a linked scan with '
                             'none; 0 (testing) ranks SeriesTime over the '
                             'task\'s series in pcinfo.csv instead (default: 1)')
    return parser.parse_args()


def bids_sub(pguid):
    return 'sub-' + pguid.replace('NDAR_INV', '')


def bids_ses(event_name):
    ses = EVENT_TO_SES.get(event_name)
    if ses is None:
        raise ValueError('unknown EventName for BIDS ses- mapping: '
                         + repr(event_name))
    return 'ses-' + ses


def pguid_from_bids(participant_id):
    return 'NDAR_INV' + participant_id.replace('sub-', '')


def event_from_bids(session_id):
    event = SES_TO_EVENT.get(session_id.replace('ses-', ''))
    if event is None:
        raise ValueError('unknown session_id for EventName mapping: '
                         + repr(session_id))
    return event


def bids_outstem(participant_id, session_id, task):
    return '{0}_{1}_task-{2}'.format(participant_id, session_id, task.lower())


def modality(task, no_scan_match=False):
    """REC is behavioral-only; so is any session with no matched scan series."""
    return 'beh' if task.upper() == 'REC' or no_scan_match else 'func'


def read_header(fname_csv):
    """Field names of a CSV, without reading the rest of the file."""
    with open(fname_csv) as f:
        return next(csv.reader(f))


def grep_rows(fname_csv, pattern):
    """Rows of fname_csv whose text contains pattern, as dicts.

    grep first, csv parse second: a full parse of the standard metadata files
    costs tens of seconds, against well under a second for grep. The header is
    read separately because grep never matches it -- it carries no pGUID.
    Matched lines go through the csv module rather than a split on commas,
    since fields such as SeriesDescription contain commas inside quotes.
    """
    fieldnames = read_header(fname_csv)
    proc = subprocess.Popen(['grep', '-F', pattern, fname_csv],
                            stdout=subprocess.PIPE)
    out = proc.communicate()[0].decode('utf-8', 'replace')
    if proc.returncode not in (0, 1):     # 1 is grep's "no lines matched"
        raise RuntimeError('grep failed on ' + fname_csv)
    if not out:
        return []
    return list(csv.DictReader(out.splitlines(), fieldnames=fieldnames))


def load_scanner_info(fname_csv, pguid, rows=None):
    """VisitID -> (Manufacturer, SoftwareVersion) for one participant.

    VisitID is per-series (the _seNNNN suffix) and carried by both metadata
    files, so each E-Prime row joins to its own scan's scanner info. Rows with
    no imaging have an empty VisitID, are absent from pcinfo.csv, and so do not
    join -- their scanner info stays unset, which is the no_scan_match case
    that routes output to beh/ anyway.
    """
    if rows is None:
        rows = grep_rows(fname_csv, pguid)
    info = {}
    for row in rows:
        visit_id = (row.get('VisitID') or '').strip()
        if visit_id and visit_id not in info:
            info[visit_id] = (row.get('Manufacturer') or '',
                              row.get('SoftwareVersion') or '')
    return info


def session_scanner_info(rows, scanner_info, verbose=True):
    """Single (Manufacturer, SoftwareVersion) for a session's rows."""
    manu = []
    sw = []
    for row in rows:
        pair = scanner_info.get((row.get('VisitID') or '').strip())
        if not pair:
            continue
        if pair[0] and pair[0].strip() and pair[0] not in manu:
            manu.append(pair[0])
        if pair[1] and pair[1].strip() and pair[1] not in sw:
            sw.append(pair[1])
    if len(manu) > 1 and verbose:
        print('WARNING: multiple Manufacturer values %s -- using first' % manu)
    if len(sw) > 1 and verbose:
        print('WARNING: multiple SoftwareVersion values %s -- using first' % sw)
    return (manu[0] if manu else None), (sw[0] if sw else None)


def session_rows(rows, pguid, event):
    """One session's rows out of a participant's grepped eprime.csv rows."""
    return [r for r in rows
            if r.get('pGUID') == pguid and r.get('EventName') == event]


def task_rows(rows, task):
    """One task's rows out of a session's rows."""
    return [r for r in rows if r.get('taskname') == SOURCE_TASK[task]]


def load_bold_runs(output_rootdir, participant_id, session_id, task,
                   bold_scanner=None):
    """SeriesInstanceUID -> BIDS run number, from the rawdata bold jsons.

    If bold_scanner is a dict, it is filled with SeriesInstanceUID ->
    (Manufacturer, SoftwareVersions) from the same jsons, the scan's own
    scanner info for skipTRs.
    """
    pattern = os.path.join(output_rootdir, participant_id, session_id, 'func',
                           '{0}_run-*_bold.json'.format(
                               bids_outstem(participant_id, session_id, task)))
    bold_runs = {}
    for fname_json in glob.glob(pattern):
        m = re.search(r'_run-(\d+)_bold\.json$', fname_json)
        if not m:
            continue
        try:
            with open(fname_json) as f:
                info = json.load(f)
            uid = (info.get('SeriesInstanceUID') or '').strip()
        except (IOError, ValueError):
            continue
        if uid:
            bold_runs[uid] = int(m.group(1))
            if bold_scanner is not None:
                bold_scanner[uid] = (info.get('Manufacturer') or '',
                                     info.get('SoftwareVersions') or '')
    return bold_runs


def pcinfo_scan_runs(pcinfo_rows, event, task):
    """SeriesInstanceUID -> rank of SeriesTime over the task's series.

    For testing with no rawdata bold files (--require_bold 0): every series of
    the task's SeriesType in pcinfo.csv, complete or not, compliant or not.
    """
    series_type = 'fMRI_{0}_task'.format(SOURCE_TASK[task])
    times = {}
    for row in pcinfo_rows:
        if row.get('EventName') != event or row.get('SeriesType') != series_type:
            continue
        uid = (row.get('SeriesInstanceUID') or '').strip()
        try:
            series_time = float(row.get('SeriesTime'))
        except (TypeError, ValueError):
            continue
        if uid and uid not in times:
            times[uid] = series_time
    ranked = sorted(times, key=lambda u: times[u])
    return dict((uid, i + 1) for i, uid in enumerate(ranked))


def build_candidates(task, rows, scanner_info, bold_runs=None,
                     bold_scanner=None, verbose=True):
    """Return one candidate dict per distinct E-Prime file for this task.

    Each candidate carries what eprime_to_events() needs: the file name, the
    WM file (REC), the run_numbers/run_times mapping, the scanner info, and
    whether output routes to beh/ instead of func/. An empty list means the
    task has no usable rows, which the caller treats as a skip, not an error.

    bold_runs maps SeriesInstanceUID to the BIDS run number: from the rawdata
    bold jsons (load_bold_runs), or for testing a pcinfo.csv SeriesTime rank
    (pcinfo_scan_runs). A file with no row linked to a scan (no
    SeriesInstanceUID, or eprime_behav_only=1) routes to beh/. A linked
    scan absent from bold_runs gets no events file and is listed in the
    candidate's missing_scans, which convert_task reports.
    """
    if not rows:
        return []

    if task.upper() == 'REC':
        # The nBack row's own file is the WM file and its extra file is REC.
        # Both keep their eprime.csv names; only their role here differs. REC
        # has no run entity, so one candidate per distinct file.
        cands = []
        seen = set()
        for row in rows:
            fname_rec = (row.get('eprime_extra_file_name') or '').strip()
            fname_wm = (row.get('eprime_file_name') or '').strip()
            if not fname_rec or fname_rec.upper() == 'NA' or fname_rec in seen:
                continue
            seen.add(fname_rec)
            cands.append({
                'fname': fname_rec,
                'fname_wm': fname_wm or None,
                'run_numbers': None,
                'run_times': None,
                'no_scan_match': False,   # REC always routes to beh/
                'Manufacturer': None,
                'SoftwareVersion': None,
            })
        return cands

    if bold_runs is None:
        bold_runs = {}

    df = pd.DataFrame(rows)

    # a row is linked to a scan only by its SeriesInstanceUID; without one,
    # eprime.csv copies the session's last SeriesTime, which is not its scan
    df['_behav_only'] = pd.to_numeric(df['eprime_behav_only'],
                                      errors='coerce') == 1
    df['_uid'] = df['SeriesInstanceUID'].fillna('').astype(str).str.strip()
    df['_linked'] = (df['_uid'] != '') & ~df['_behav_only']

    # scanner info for skipTRs: from the rawdata bold json of a file's linked
    # scan (bold_scanner), else of any of the task's scans, else pcinfo.csv,
    # whose SoftwareVersion is empty for GE DV25 (skipTRs 15 there, not 4)
    if bold_scanner is None:
        bold_scanner = {}
    manufacturer, software_version = session_scanner_info(rows, scanner_info,
                                                          verbose=verbose)
    task_scanner = [v for v in bold_scanner.values() if v[0]]
    if task_scanner:
        manufacturer, software_version = task_scanner[0]

    cands = []
    fnames = [f for f in df['eprime_file_name'].unique()
              if str(f).strip().upper() not in ('', 'NA')]
    # beh/ runs are numbered on from one behavioral-only file to the next, so
    #   two such files for one task do not both write run-01
    beh_next = 1
    for fname in fnames:
        rows_f = df[df['eprime_file_name'] == fname]
        nruns_vals = pd.to_numeric(rows_f['eprime_nruns'],
                                   errors='coerce').dropna()
        if not len(nruns_vals):
            if verbose:
                print('WARNING: no eprime_nruns for %s, skipping' % fname)
            continue
        eprime_nruns = int(nruns_vals.iloc[0])
        missing_scans = []

        behav_only = (bool(rows_f['_behav_only'].any())
                      or not bool(rows_f['_linked'].any()))
        if behav_only:
            # no scan: every run, numbered by its E-Prime run position
            if verbose:
                print('%s: NOTE: %s has no scan linked in eprime.csv -- '
                      'writing all %d run(s) to beh/ by eprime_run position'
                      % (SCRIPTNAME, os.path.basename(fname), eprime_nruns))
            run_numbers = list(range(beh_next, beh_next + eprime_nruns))
            run_times = None
        else:
            run_numbers = [0] * eprime_nruns
            run_times = [None] * eprime_nruns
            for _, row in rows_f[rows_f['_linked']].iterrows():
                try:
                    eprime_run = int(float(str(row['eprime_run']).strip()))
                except (ValueError, TypeError):
                    continue
                if not 1 <= eprime_run <= eprime_nruns:
                    continue
                bids_run = bold_runs.get(row['_uid'])
                if bids_run is None:
                    missing_scans.append(row['_uid'])
                    continue
                if run_numbers[eprime_run - 1] not in (0, bids_run):
                    print('WARNING: %s run %d: linked to BIDS runs %d and %d '
                          '-- keeping %d' % (os.path.basename(fname),
                                             eprime_run,
                                             run_numbers[eprime_run - 1],
                                             bids_run,
                                             run_numbers[eprime_run - 1]))
                    continue
                if bids_run in run_numbers[:eprime_run - 1] + \
                        run_numbers[eprime_run:]:
                    print('WARNING: %s run %d: BIDS run %d already taken by '
                          'another run of this file -- skipped'
                          % (os.path.basename(fname), eprime_run, bids_run))
                    continue
                run_numbers[eprime_run - 1] = bids_run
                run_times[eprime_run - 1] = row.get('eprime_time')
        file_scanner = [bold_scanner[u] for u in rows_f.loc[rows_f['_linked'], '_uid']
                        if u in bold_scanner and bold_scanner[u][0]]
        file_manufacturer, file_software_version = (
            file_scanner[0] if file_scanner else (manufacturer, software_version))
        # with no scanner info, skipTRs is unknown and a GE-version file's
        #   onsets would stay timed from the end of the dummy volumes rather
        #   than the first stored one, so its runs go to beh/ (behavioral
        #   metrics only), numbered by E-Prime run position
        if not behav_only and not (file_manufacturer or file_software_version):
            print('WARNING: %s: no scanner info (bold json, task scan or '
                  'pcinfo.csv) -- writing all %d run(s) to beh/ by eprime_run '
                  'position' % (os.path.basename(fname), eprime_nruns))
            behav_only = True
            run_numbers = list(range(beh_next, beh_next + eprime_nruns))
            run_times = None
        if behav_only:
            beh_next += eprime_nruns
        cands.append({
            'fname': fname,
            'fname_wm': None,
            'run_numbers': run_numbers,
            'run_times': None if behav_only else run_times,
            'no_scan_match': behav_only,
            'Manufacturer': file_manufacturer,
            'SoftwareVersion': file_software_version,
            'missing_scans': missing_scans,
        })
    return cands


def convert_task(participant_id, session_id, task, cands, output_rootdir,
                 verbose=True, forceflag=False, require_bold=True):
    """Convert one task's candidates. Returns (errcode, errmsg, status).

    status is 'converted', 'skipped' (no rows for this task) or 'failed'. A
    file named in the metadata but absent from disk is a failure, unlike a task
    that was simply never administered; so is a scan eprime.csv links that has
    no rawdata bold json, unless require_bold is False, when it is a WARNING.
    A file's runs with run numbers are written either way.
    """
    if not cands:
        if verbose:
            print('%s: NOTE: no %s rows for %s %s -- skipping'
                  % (SCRIPTNAME, task, participant_id, session_id))
        return 0, None, 'skipped'

    outstem = bids_outstem(participant_id, session_id, task)

    nfail = 0
    nmissing = 0
    nmissing_scans = 0
    outdirs = []
    for cand in cands:
        fname = cand['fname']
        for uid in cand.get('missing_scans') or []:
            if require_bold:
                print('ERROR: %s: linked scan %s has no rawdata bold json'
                      % (os.path.basename(fname), uid), file=sys.stderr)
                nmissing_scans += 1
            else:
                print('WARNING: %s: linked scan %s is not in pcinfo.csv -- no '
                      'events file written for it'
                      % (os.path.basename(fname), uid))
        if cand['run_numbers'] is not None and not any(cand['run_numbers']):
            continue
        outdir = os.path.join(output_rootdir, participant_id, session_id,
                              modality(task, cand['no_scan_match']))
        if not os.path.isfile(fname):
            print('ERROR: E-Prime file listed but not on disk: %s' % fname,
                  file=sys.stderr)
            nmissing += 1
            continue
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        if outdir not in outdirs:
            outdirs.append(outdir)
        n_runs, runs_ok, errcode, errmsg = eprime_to_events(
            fname, outdir=outdir, outstem=outstem, task=task,
            fname_wm=cand['fname_wm'],
            run_numbers=cand['run_numbers'],
            run_times=cand['run_times'],
            Manufacturer=cand['Manufacturer'],
            SoftwareVersion=cand['SoftwareVersion'],
            forceflag=forceflag, verbose=verbose)
        if errcode:
            print('ERROR: %s: %s' % (os.path.basename(fname), errmsg),
                  file=sys.stderr)
            nfail += 1
        elif verbose:
            print('%s: %s: %d run(s), ok=%s' % (task, os.path.basename(fname),
                                                n_runs, runs_ok))

    if nmissing or nfail or nmissing_scans:
        return 1, '{0}: {1} of {2} E-Prime file(s) missing on disk, {3} ' \
                  'failed, {4} linked scan(s) not in rawdata'.format(
                      task, nmissing, len(cands), nfail,
                      nmissing_scans), 'failed'
    if verbose:
        print('%s: output written to %s' % (task, ', '.join(outdirs)))
    return 0, None, 'converted'


def run(args):
    for fname_csv in (args.eprime_csv, args.pcinfo_csv):
        if not os.path.isfile(fname_csv):
            return 1, 'metadata CSV not found: ' + fname_csv

    pguid = pguid_from_bids(args.participant_id)
    try:
        event = event_from_bids(args.session_id)
    except ValueError as e:
        return 1, str(e)

    tasks = []
    for task in args.task:
        if task not in tasks:
            tasks.append(task)

    try:
        rows = session_rows(grep_rows(args.eprime_csv, pguid), pguid, event)
    except RuntimeError as e:
        return 1, str(e)
    if not rows:
        return 1, 'no rows in {0} for {1} {2}'.format(
            os.path.basename(args.eprime_csv), pguid, event)

    # REC takes no scanner info, so skip the larger grep when REC is all we want
    scanner_info = {}
    pcinfo_rows = []
    if [t for t in tasks if t.upper() != 'REC']:
        try:
            pcinfo_rows = grep_rows(args.pcinfo_csv, pguid)
        except RuntimeError as e:
            return 1, str(e)
        scanner_info = load_scanner_info(args.pcinfo_csv, pguid,
                                         rows=pcinfo_rows)

    converted = []
    skipped = []
    errmsgs = []
    for task in tasks:
        bold_runs = {}
        bold_scanner = {}
        if task.upper() != 'REC' and args.require_bold:
            bold_runs = load_bold_runs(args.output_rootdir,
                                       args.participant_id, args.session_id,
                                       task, bold_scanner=bold_scanner)
        elif task.upper() != 'REC':
            bold_runs = pcinfo_scan_runs(pcinfo_rows, event, task)
        cands = build_candidates(task, task_rows(rows, task), scanner_info,
                                 bold_runs=bold_runs,
                                 bold_scanner=bold_scanner,
                                 verbose=bool(args.verbose))
        errcode, errmsg, status = convert_task(
            args.participant_id, args.session_id, task, cands,
            args.output_rootdir, verbose=bool(args.verbose),
            forceflag=bool(args.forceflag),
            require_bold=bool(args.require_bold))
        if status == 'converted':
            converted.append(task)
        elif status == 'skipped':
            skipped.append(task)
        else:
            errmsgs.append(errmsg)

    if args.verbose:
        print('%s: %s %s: converted %s; skipped %s'
              % (SCRIPTNAME, args.participant_id, args.session_id,
                 ', '.join(converted) or 'none', ', '.join(skipped) or 'none'))
    if errmsgs:
        return 1, '; '.join(errmsgs)
    return 0, None


def main():
    errcode, errmsg = run(parse_args())
    if errmsg:
        print(errmsg, file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
