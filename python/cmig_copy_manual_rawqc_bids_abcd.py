#!/usr/bin/env python3
"""
Copy ABCD manual raw QC scores from DAIC merged QC csv to BIDS derivatives.

Created:  08/13/2026 by Don Hagler
Prev mod: 08/13/2026 by Don Hagler
Prev mod: 08/21/2026 by Don Hagler
Prev mod: 09/22/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Last mod: 11/12/2026 by Don Hagler
"""

import argparse
import csv
import datetime
import json
import os
import re
import sys

import cmig_sing_utils

APPNAME = 'cmig_copy_manual_rawqc_bids_abcd'

# csv key columns kept for each series row
KEY_COLUMNS = [
    'pGUID', 'EventName', 'SeriesType', 'SeriesTime', 'AcquisitionTime',
    'StudyDate', 'StudyInstanceUID', 'SeriesInstanceUID',
]

# manual QC columns kept for each series row (non-empty values only)
#   consensus/summary: QC, QU_motion, nrev, revdisp
#   per-reviewer: QC_XX, QU_XX_motion / QU_motion_XX, notes_XX, started_XX,
#     completed_XX, duration_XX, revdate_XX
QC_COLUMN_PATTERN = \
    r'^(QC$|QC_|QU_|notes_|started_|completed_|duration_|revdate_|nrev$|revdisp$)'

# scans.tsv type -> csv SeriesType
SERIES_TYPE_MAP = {
    'ABCD-T1':              'T1',
    'ABCD-T1-NORM':         'T1_NORM',
    'ABCD-T2':              'T2',
    'ABCD-T2-NORM':         'T2_NORM',
    'ABCD-rsfMRI':          'rsfMRI',
    'ABCD-MID-fMRI':        'fMRI_MID_task',
    'ABCD-SST-fMRI':        'fMRI_SST_task',
    'ABCD-nBack-fMRI':      'fMRI_nBack_task',
    'ABCD-fMRI-FM':         'fMRI_FM',
    'ABCD-fMRI-FM-AP':      'fMRI_FM_AP',
    'ABCD-fMRI-FM-PA':      'fMRI_FM_PA',
    'ABCD-DTI':             'dMRI',
    'ABCD-Diffusion-FM':    'dMRI_FM',
    'ABCD-Diffusion-FM-AP': 'dMRI_FM_AP',
    'ABCD-Diffusion-FM-PA': 'dMRI_FM_PA',
}

# series type -> modality group, for notes-code parsing
SMRI_TYPES = set(['T1', 'T1_NORM', 'T2', 'T2_NORM'])
DMRI_TYPES = set(['dMRI', 'dMRI_FM', 'dMRI_FM_AP', 'dMRI_FM_PA'])
FMRI_TYPES = set(['rsfMRI', 'fMRI_MID_task', 'fMRI_SST_task',
                  'fMRI_nBack_task', 'fMRI_FM', 'fMRI_FM_AP', 'fMRI_FM_PA'])

# notes code -> canonical issue, following parse_notes in
#   autoqc_load_raw_qcinfo.m: a horizontal band, a bad slice, motion and a
#   multiband recon artifact all present as line artifacts on a diffusion or
#   functional scan, and distortion and cutoff each have several spellings.
#   gh and moire are the two departures from that function, which routes
#   neither anywhere; both are written often enough in the notes to parse
NOTES_CODE_ISSUE = {
    'dis': 'sus', 'sd': 'sus', 'si': 'sus', 'sus': 'sus',
    'dco': 'fov', 'vco': 'fov', 'wr': 'fov', 'fov': 'fov',
    'hb': 'line', 'SLICE': 'line', 'line': 'line', 'mo': 'line',
    'mbr': 'line',
    'gh': 'gh', 'fa': 'fa',
    'moire': 'moire',
}

# spellings reviewers use for a code, beyond the code itself: the artifact
#   written out, a typo seen repeatedly, and a transposition. bh and recon
#   resolve to line, which is where hb and mbr go. moir cannot be reached by
#   the moire searches, which need the final e, and cannot double count: the
#   bare form requires a non-letter after the code
NOTES_CODE_ALIAS = {
    'ghost': 'gh', 'ghpst': 'gh', 'bh': 'line', 'recon': 'line',
    'moir': 'moire',
}

# the artifacts each modality scores. sMRI carries the reviewer's recorded
#   motion rating and parses nothing: every sMRI review has one, and where a
#   note also carries a code the two agree in 89 per cent of cases, so a
#   derived score would only overwrite an explicit one
QU_ISSUES = {
    'sMRI': ['mo'],
    'dMRI': ['fov', 'sus', 'line', 'fa', 'gh'],
    'fMRI': ['fov', 'sus', 'line', 'gh', 'moire'],
}

# canonical issue -> output QU column
QU_COLUMN_MAP = {
    'mo': 'QU_motion', 'sus': 'QU_sus', 'fov': 'QU_cutoff',
    'line': 'QU_line', 'fa': 'QU_fa', 'gh': 'QU_ghost',
    'moire': 'QU_moire',
}

# appended rather than ordered with the rest, so the existing columns keep
#   their positions
QU_COLUMNS = ['QU_motion', 'QU_sus', 'QU_cutoff', 'QU_line', 'QU_fa',
              'QU_ghost', 'QU_moire']

# a bare code asserts the artifact without grading it. On a failed scan where
#   nothing else on the note is numbered it has to be what failed the scan;
#   otherwise something else accounts for the outcome
BARE_SCORE_FAILED = 3
BARE_SCORE_OTHER = 1.5

# the severity at which an artifact fails the scan
SEVERE_SCORE = 3

# how a reviewer's QU scores for a scan were arrived at
QU_RECORDED = 'recorded'
QU_PARSED = 'parsed'

# output detail file columns (long format: one row per scan x reviewer)
DETAIL_COLUMNS = [
    'participant_id', 'session_id', 'filename', 'series_type',
    'reviewer_id', 'source', 'QC',
] + QU_COLUMNS + [
    'qu_source', 'notes', 'started', 'completed', 'duration', 'revdate',
    'series_uid',
]

# output summary file columns (one row per scan)
SUMMARY_COLUMNS = [
    'participant_id', 'session_id', 'filename', 'series_type',
    'QC',
] + QU_COLUMNS + ['qu_source', 'nrev', 'revdisp', 'notes', 'series_uid']

# the summary file is compiled across sessions into the release table
#   mr_y_qc__raw__man__series, so its header carries release names; the rows
#   keep the short keys, and the ids and filename are unchanged
SUMMARY_TABLE = 'mr_y_qc__raw__man__series'
SUMMARY_NAMES = {
    'series_type': SUMMARY_TABLE + '_type',
    'QC':          SUMMARY_TABLE + '_score',
    'QU_motion':   SUMMARY_TABLE + '__mot_score',
    'QU_sus':      SUMMARY_TABLE + '__susc_score',
    'QU_cutoff':   SUMMARY_TABLE + '__cutoff_score',
    'QU_line':     SUMMARY_TABLE + '__line_score',
    'QU_fa':       SUMMARY_TABLE + '__famap_score',
    'QU_ghost':    SUMMARY_TABLE + '__ghost_score',
    'QU_moire':    SUMMARY_TABLE + '__moire_score',
    'qu_source':   SUMMARY_TABLE + '__score_source',
    'nrev':        SUMMARY_TABLE + '__rev_count',
    'revdisp':     SUMMARY_TABLE + '__rev__disp_indicator',
    'notes':       SUMMARY_TABLE + '_notes',
    'series_uid':  SUMMARY_TABLE + '_uid',
}

# the detail file, one row per scan and reviewer, is compiled into
#   mr_y_qc__raw__man__series__review likewise; rows keep the short keys
DETAIL_TABLE = 'mr_y_qc__raw__man__series__review'
DETAIL_NAMES = {
    'series_type': DETAIL_TABLE + '_type',
    'source':      DETAIL_TABLE + '_source',
    'QC':          DETAIL_TABLE + '_score',
    'QU_motion':   DETAIL_TABLE + '__mot_score',
    'QU_sus':      DETAIL_TABLE + '__susc_score',
    'QU_cutoff':   DETAIL_TABLE + '__cutoff_score',
    'QU_line':     DETAIL_TABLE + '__line_score',
    'QU_fa':       DETAIL_TABLE + '__famap_score',
    'QU_ghost':    DETAIL_TABLE + '__ghost_score',
    'QU_moire':    DETAIL_TABLE + '__moire_score',
    'qu_source':   DETAIL_TABLE + '__score_source',
    'notes':       DETAIL_TABLE + '_notes',
    'started':     DETAIL_TABLE + '__start_dtt',
    'completed':   DETAIL_TABLE + '__end_dtt',
    'duration':    DETAIL_TABLE + '_dur',
    'revdate':     DETAIL_TABLE + '__save_dtt',
    'series_uid':  DETAIL_TABLE + '_uid',
}
# release name -> short key, for reading; short names (files written before
#   the rename) are kept as they are
DETAIL_KEYS = dict((v, k) for k, v in DETAIL_NAMES.items())

SOURCE_METADATA = 'metadata'
SOURCE_REVIEWER = 'reviewer'
NA = 'n/a'

#############################################################################
# ID and session mapping
#############################################################################

def pid_to_pguid(participant_id):
    # 'sub-XXXXXXXX' -> 'NDAR_INVXXXXXXXX'
    return 'NDAR_INV' + re.sub(r'^sub-', '', participant_id)


def pguid_to_pid(pguid):
    # 'NDAR_INVXXXXXXXX' -> 'sub-XXXXXXXX'
    return 'sub-' + re.sub(r'^NDAR_INV', '', pguid)


def sid_to_eventname(session_id):
    # 'ses-00A' -> 'baseline_year_1_arm_1'
    # 'ses-02A' -> '2_year_follow_up_y_arm_1', etc.
    m = re.match(r'^ses-(\d+)[A-Z]$', session_id)
    if not m:
        return None
    year = int(m.group(1))
    if year == 0:
        return 'baseline_year_1_arm_1'
    return '{0}_year_follow_up_y_arm_1'.format(year)


def eventname_to_sid(eventname):
    # 'baseline_year_1_arm_1' -> 'ses-00A'
    # '2_year_follow_up_y_arm_1' -> 'ses-02A', etc.
    if eventname == 'baseline_year_1_arm_1':
        return 'ses-00A'
    m = re.match(r'^(\d+)_year_follow_up_y_arm_1$', eventname)
    if not m:
        return None
    return 'ses-{0:02d}A'.format(int(m.group(1)))


def norm_csv_datetime(study_date, series_time):
    # StudyDate '20180625' + SeriesTime '143326.81400000' -> '20180625T143326'
    # SeriesTime is truncated (not rounded) to whole seconds
    if not study_date or not series_time:
        return None
    hhmmss = series_time.split('.')[0].zfill(6)
    return '{0}T{1}'.format(study_date, hhmmss)


def norm_acq_time(acq_time):
    # '2018-06-25T14:33:26' -> '20180625T143326'
    if not acq_time:
        return None
    return re.sub(r'[-:]', '', acq_time.strip())

#############################################################################
# value normalization
#############################################################################

def norm_int(val):
    # '1 ' -> '1'; empty or unparseable -> None
    if val is None or not str(val).strip():
        return None
    try:
        return int(round(float(str(val).strip())))
    except ValueError:
        return None


def norm_float(val):
    # '36.78648200 ' -> 36.786482; empty or unparseable -> None
    if val is None or not str(val).strip():
        return None
    try:
        return float(str(val).strip())
    except ValueError:
        return None


def fmt_float(val):
    # minimal formatting: 0.5 not 0.50000000; 2.0 -> 2
    if val is None:
        return NA
    txt = '{0:.8f}'.format(val).rstrip('0').rstrip('.')
    return txt if txt else '0'


def norm_timestamp(val):
    # 'yyyymmdd.HHMMSSFF' -> ISO 8601 'YYYY-MM-DDTHH:MM:SS'
    # (centiseconds dropped, matching scans.tsv acq_time style)
    if val is None or not str(val).strip():
        return None
    txt = str(val).strip()
    m = re.match(r'^(\d{4})(\d{2})(\d{2})\.(\d{2})(\d{2})(\d{2})', txt)
    if not m:
        return None
    return '{0}-{1}-{2}T{3}:{4}:{5}'.format(*m.groups())


def norm_notes(val):
    # strip outer whitespace and trailing semicolons; content preserved
    if val is None:
        return None
    txt = str(val).strip().rstrip(';').strip()
    return None if is_empty(txt) else txt


def modality_group(series_type):
    if series_type in SMRI_TYPES:
        return 'sMRI'
    if series_type in DMRI_TYPES:
        return 'dMRI'
    if series_type in FMRI_TYPES:
        return 'fMRI'
    return None


def notes_code_map(group):
    # the codes this modality reads, each on the issue it scores; a code whose
    # issue the modality does not score is left out
    issues = QU_ISSUES.get(group)
    if not issues:
        return {}
    code_map = {}
    for src in (NOTES_CODE_ISSUE, NOTES_CODE_ALIAS):
        for code, issue in src.items():
            if issue in issues:
                code_map[code] = issue
    return code_map


def parse_notes_codes(notes, series_type, failed=False):
    # parse <code><severity> codes into canonical issues, taking max severity
    # per issue, allowing the plural a reviewer sometimes writes. A code written without a digit still asserts the artifact:
    # it scores 3 on a failed scan where no code is numbered, since it is then
    # what failed the scan, and 1.5 otherwise
    scores = {}
    code_map = notes_code_map(modality_group(series_type))
    if not notes or not code_map:
        return scores
    bare = {}
    numbered = False
    for code, issue in code_map.items():
        found = False
        for m in re.finditer(re.escape(code) + r's?\s*([1-3])', notes,
                             re.IGNORECASE):
            found = True
            numbered = True
            val = int(m.group(1))
            if val > scores.get(issue, 0):
                scores[issue] = val
        # a bare code has to stand as a word: mo is inside MOIRE and fa is
        # inside artifact, and either would otherwise score. The numbered
        # search needs no such guard, a digit having to follow
        if not found and re.search(
                r'(?<![A-Za-z])' + re.escape(code) + r's?(?![A-Za-z])',
                notes, re.IGNORECASE):
            bare[issue] = True
    bare_score = BARE_SCORE_FAILED if (failed and not numbered) \
        else BARE_SCORE_OTHER
    for issue in bare:
        if issue not in scores:
            scores[issue] = bare_score
    return scores

#############################################################################
# load QC csv (TODO 1a)
#############################################################################

def load_qcinfo(fname, pguid_filter=None, verbose=0):
    # read csv once; keep key columns plus non-empty manual QC columns
    # index by acq datetime (primary) and SeriesInstanceUID (fallback)
    # also collect the set of sessions that have any manual QC, which
    # drives session enumeration (no rawdata globbing needed)
    uid_index = {}
    time_index = {}
    qc_sessions = set()
    nrows = 0
    nkept = 0
    print('{0}: loading QC info from {1}...'.format(APPNAME, fname))
    with open(fname, newline='') as f:
        reader = csv.DictReader(f)
        qc_cols = [c for c in reader.fieldnames
                   if re.match(QC_COLUMN_PATTERN, c)]
        # reviewer ids from the per-reviewer QC columns (QC_XX)
        reviewer_ids = sorted(set(
            c[3:] for c in reader.fieldnames if re.match(r'^QC_\w+$', c)))
        rev_qc_cols = ['QC_{0}'.format(r) for r in reviewer_ids]
        for row in reader:
            nrows += 1
            pguid = (row.get('pGUID') or '').strip()
            if pguid_filter is not None and pguid not in pguid_filter:
                continue
            keep = {}
            for c in KEY_COLUMNS:
                keep[c] = (row.get(c) or '').strip()
            for c in qc_cols:
                val = (row.get(c) or '').strip()
                if val:
                    keep[c] = val
            uid = keep['SeriesInstanceUID']
            # a row with no pGUID is a stale QC container, and the uid
            # fallback would otherwise match a scan to one
            if uid and pguid:
                if uid in uid_index:
                    print('{0}: WARNING: duplicate SeriesInstanceUID '
                          'in csv: {1}'.format(APPNAME, uid))
                else:
                    uid_index[uid] = keep
            dt = norm_csv_datetime(keep['StudyDate'], keep['SeriesTime'])
            if pguid and dt:
                time_index.setdefault((pguid, dt), []).append(keep)
            # record sessions with at least one manual review
            if pguid and any(c in keep for c in rev_qc_cols):
                qc_sessions.add((pguid, keep['EventName']))
            nkept += 1
    print('{0}: read {1} csv rows, kept {2}; {3} sessions with manual '
          'QC'.format(APPNAME, nrows, nkept, len(qc_sessions)))
    return uid_index, time_index, reviewer_ids, qc_sessions

#############################################################################
# enumerate scans from rawdata scans.tsv (TODO 1b)
#############################################################################

def get_sessions(raw_rootdir, qc_sessions, participant_id=None,
                 session_id=None, verbose=0):
    # build the session list from the csv sessions that have manual QC,
    # constructing each scans.tsv path directly
    # globbing the rawdata tree (~12k subject dirs) is far slower and
    # would also enumerate sessions with nothing to copy
    pid_filter = set(participant_id) if participant_id else None
    sessions = []
    nmissing = 0
    nbadname = 0
    for pguid, eventname in qc_sessions:
        pid = pguid_to_pid(pguid)
        sid = eventname_to_sid(eventname)
        if sid is None:
            nbadname += 1
            continue
        if pid_filter is not None and pid not in pid_filter:
            continue
        if session_id is not None and sid != session_id:
            continue
        fpath = os.path.join(raw_rootdir, pid, sid,
                             '{0}_{1}_scans.tsv'.format(pid, sid))
        if not os.path.isfile(fpath):
            # in the csv but not in rawdata (e.g. 8-year events not yet
            # converted, or quarantined sessions)
            nmissing += 1
            continue
        sessions.append((pid, sid, fpath))
    sessions.sort()
    if nmissing:
        print('{0}: {1} sessions with manual QC have no scans.tsv in '
              'rawdata (skipped)'.format(APPNAME, nmissing))
    if nbadname:
        print('{0}: WARNING: {1} sessions with unrecognized EventName '
              '(skipped)'.format(APPNAME, nbadname))
    return sessions


def read_scans(scans_path):
    # read scans.tsv; return list of row dicts
    with open(scans_path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return [row for row in reader]


def get_sidecar_uids(sess_dir, filename):
    # read JSON sidecar for a scan; return (SeriesInstanceUID,
    # StudyInstanceUID)
    json_path = os.path.join(
        sess_dir, re.sub(r'\.nii(\.gz)?$', '.json', filename))
    if not os.path.isfile(json_path):
        return None, None
    try:
        with open(json_path) as f:
            meta = json.load(f)
    except ValueError:
        print('{0}: WARNING: failed to parse JSON sidecar {1}'.format(
            APPNAME, json_path))
        return None, None
    return meta.get('SeriesInstanceUID'), meta.get('StudyInstanceUID')

#############################################################################
# match scans to csv rows (TODO 1c)
#############################################################################

def match_session(pid, sid, scans_path, uid_index, time_index, verbose=0):
    # match each scans.tsv row to a csv row
    # returns list of (scan, csvrow, method) and per-session counts
    # acq_time is tried first because it comes from scans.tsv, already in
    # hand; the SeriesInstanceUID fallback costs a JSON sidecar read per
    # scan, which dominates runtime over the whole dataset
    sess_dir = os.path.dirname(scans_path)
    pguid = pid_to_pguid(pid)
    results = []
    counts = {'uid': 0, 'time': 0, 'unmatched': 0}
    for scan in read_scans(scans_path):
        fname = (scan.get('filename') or '').strip()
        stype = SERIES_TYPE_MAP.get((scan.get('type') or '').strip())
        csvrow = None
        method = None
        # primary match: acq_time == StudyDate + SeriesTime (truncated)
        dt = norm_acq_time(scan.get('acq_time'))
        rows = time_index.get((pguid, dt), [])
        if len(rows) == 1:
            csvrow = rows[0]
            method = 'time'
        elif len(rows) > 1 and stype:
            cand = [r for r in rows if r['SeriesType'] == stype]
            if len(cand) == 1:
                csvrow = cand[0]
                method = 'time+type'
        if csvrow is None:
            # fallback match: sidecar SeriesInstanceUID
            uid, study_uid = get_sidecar_uids(sess_dir, fname)
            # a scan with no csv row has no uid from the csv, so keep this one
            if uid:
                scan['series_uid'] = uid
            if uid and uid in uid_index:
                csvrow = uid_index[uid]
                method = 'uid'
        if csvrow is None:
            counts['unmatched'] += 1
            print('{0}: WARNING: no csv match for {1} {2} {3}'.format(
                APPNAME, pid, sid, fname))
        else:
            counts['uid' if method == 'uid' else 'time'] += 1
            # sanity checks
            if stype and csvrow['SeriesType'] != stype:
                print('{0}: WARNING: series type mismatch for {1} {2} {3}: '
                      'scans.tsv {4} vs csv {5}'.format(
                          APPNAME, pid, sid, fname,
                          stype, csvrow['SeriesType']))
            if csvrow['pGUID'] != pguid:
                print('{0}: WARNING: pGUID mismatch for {1} {2} {3}: '
                      '{4}'.format(APPNAME, pid, sid, fname, csvrow['pGUID']))
            eventname = sid_to_eventname(sid)
            if eventname and csvrow['EventName'] and \
               csvrow['EventName'] != eventname:
                print('{0}: WARNING: EventName mismatch for {1} {2} {3}: '
                      'expected {4}, csv has {5}'.format(
                          APPNAME, pid, sid, fname,
                          eventname, csvrow['EventName']))
        results.append((scan, csvrow, method))
    if verbose:
        print('{0}: {1} {2}: {3} scans, {4} matched by uid, '
              '{5} matched by time, {6} unmatched'.format(
                  APPNAME, pid, sid, len(results),
                  counts['uid'], counts['time'], counts['unmatched']))
    return results, counts

#############################################################################
# build output rows (TODO 1d)
#############################################################################

def get_reviewer_value(csvrow, prefix, revid):
    # per-reviewer column, e.g. notes_CP
    return csvrow.get('{0}_{1}'.format(prefix, revid))


def get_reviewer_qu_motion(csvrow, revid):
    # reviewer block naming is inconsistent: QU_XX_motion and
    # QU_motion_XX both occur; treat as the same field
    for col in ('QU_{0}_motion'.format(revid), 'QU_motion_{0}'.format(revid)):
        val = csvrow.get(col)
        if val is not None and str(val).strip():
            return val
    return None


def build_detail_rows(pid, sid, scan, csvrow, reviewer_ids):
    # one row per reviewer with any data for this scan
    rows = []
    fname = (scan.get('filename') or '').strip()
    stype = csvrow.get('SeriesType') or ''
    uid = csvrow.get('SeriesInstanceUID') or ''
    for revid in reviewer_ids:
        qc = norm_int(csvrow.get('QC_{0}'.format(revid)))
        qu_motion = norm_int(get_reviewer_qu_motion(csvrow, revid))
        notes = norm_notes(get_reviewer_value(csvrow, 'notes', revid))
        started = norm_timestamp(get_reviewer_value(csvrow, 'started', revid))
        completed = norm_timestamp(
            get_reviewer_value(csvrow, 'completed', revid))
        duration = norm_float(get_reviewer_value(csvrow, 'duration', revid))
        revdate = norm_timestamp(get_reviewer_value(csvrow, 'revdate', revid))
        # skip reviewers with no data at all for this scan
        if qc is None and qu_motion is None and notes is None and \
           started is None and completed is None and duration is None and \
           revdate is None:
            continue
        # if completed is empty but revdate is not, use revdate
        if completed is None and revdate is not None:
            completed = revdate
        # QU columns: sMRI motion is the reviewer's recorded rating, every
        # other score is parsed from the notes. A column the modality does not
        # score stays empty, which is not the same as scoring it zero
        qu_vals = dict((col, None) for col in QU_COLUMNS)
        if qu_motion is not None and \
           'mo' in QU_ISSUES.get(modality_group(stype), []):
            qu_vals['QU_motion'] = qu_motion
        for issue, val in parse_notes_codes(notes, stype, qc == 0).items():
            qu_vals[QU_COLUMN_MAP[issue]] = val
        row = {
            'participant_id': pid,
            'session_id':     sid,
            'filename':       fname,
            'series_type':    stype,
            'reviewer_id':    revid,
            'source':         SOURCE_METADATA,
            'QC':             qc,
            'notes':          notes,
            'started':        started,
            'completed':      completed,
            'duration':       duration,
            'revdate':        revdate,
            'series_uid':     uid,
            # every score a legacy review can yield beyond sMRI motion comes
            # from the notes, and an empty column there means nothing was
            # written rather than nothing was seen
            'qu_source':      QU_PARSED,
        }
        row.update(qu_vals)
        rows.append(row)
    return rows


def consensus_qc(row):
    # a scan carrying an artifact the reviewer called severe does not pass:
    # the score and the verdict were entered independently and nothing
    # checked one against the other. The detail row keeps the verdict as
    # given, and the scores were read against it, so an imputed bare code,
    # which never reaches 3 on a passing scan, cannot trigger this
    qc = row['QC']
    if qc == 1 and any(row.get(col) is not None and
                       row[col] >= SEVERE_SCORE for col in QU_COLUMNS):
        return 0
    return qc


def compute_consensus(pid, sid, scan_key, rows, uid=None):
    # consensus across reviewers, following MMIL_Summarize_QC.m:
    #   QC = round(sum(QC)/nrev - 0.01) -> majority vote, ties fail
    #   QU = sum(QU, missing as 0)/nrev
    #   revdisp = 1 if any QC disagreement
    #   notes concatenated with '; ' when distinct
    #   a reviewer's QC 1 with a severe QU counts as 0 (consensus_qc)
    qc_vals = [consensus_qc(r) for r in rows if r['QC'] is not None]
    nrev = len(qc_vals)
    # the reviewers agree on how their scores were arrived at unless a scan
    # has been reviewed in both eras, where the parsed reading governs
    qu_sources = set(r.get('qu_source') for r in rows if r.get('qu_source'))
    summary = {
        'participant_id': pid,
        'session_id':     sid,
        'filename':       scan_key[0],
        'series_type':    scan_key[1],
        'nrev':           nrev,
        'series_uid':     uid,
        'qu_source':      QU_PARSED if QU_PARSED in qu_sources
                          else (QU_RECORDED if qu_sources else None),
    }
    scored = QU_ISSUES.get(modality_group(scan_key[1]), [])
    scored_cols = [QU_COLUMN_MAP[i] for i in scored]
    if nrev:
        summary['QC'] = int(round((float(sum(qc_vals)) / nrev) - 0.01))
        summary['revdisp'] = 1 if len(set(qc_vals)) > 1 else 0
        for col in QU_COLUMNS:
            vals = [r[col] for r in rows if r['QC'] is not None]
            if col in scored_cols and any(v is not None for v in vals):
                total = sum(v for v in vals if v is not None)
                summary[col] = float(total) / nrev
            else:
                summary[col] = None
    else:
        summary['QC'] = None
        summary['revdisp'] = 0
        for col in QU_COLUMNS:
            summary[col] = None
    notes = []
    for r in rows:
        if r['notes'] and r['notes'] not in notes:
            notes.append(r['notes'])
    summary['notes'] = '; '.join(notes) if notes else None
    return summary


def is_empty(val):
    # None, blank or n/a (files written before n/a have empty cells)
    return val is None or str(val).strip() in ('', NA)


def fmt_detail_value(col, val):
    if is_empty(val):
        return NA
    if col == 'duration' or col in QU_COLUMNS:
        return fmt_float(val)
    return str(val)


def fmt_summary_value(col, val):
    if is_empty(val):
        return NA
    if col in QU_COLUMNS:
        return fmt_float(val)
    return str(val)

#############################################################################
# write output files (TODO 1d) and status file (TODO 1e)
#############################################################################

def get_output_paths(output_rootdir, pid, sid):
    outdir = os.path.join(output_rootdir, pid, sid)
    detail = os.path.join(
        outdir, '{0}_{1}_desc-reviews_manualrawqc.tsv'.format(pid, sid))
    summary = os.path.join(
        outdir, '{0}_{1}_desc-series_manualrawqc.tsv'.format(pid, sid))
    status = os.path.join(
        outdir, '{0}_{1}_desc-manualrawqc_status.tsv'.format(pid, sid))
    return outdir, detail, summary, status


def get_old_output_paths(output_rootdir, pid, sid):
    # the detail and summary files under their older names
    stem = os.path.join(output_rootdir, pid, sid,
                        '{0}_{1}_desc-manualrawqc'.format(pid, sid))
    return stem + '_qc.tsv', stem + '_scans.tsv'


def get_existing_detail_path(output_rootdir, pid, sid):
    # the detail file to read: the older name where the current one is missing
    detail = get_output_paths(output_rootdir, pid, sid)[1]
    old_detail = get_old_output_paths(output_rootdir, pid, sid)[0]
    if not os.path.isfile(detail) and os.path.isfile(old_detail):
        return old_detail
    return detail


def remove_old_output(output_rootdir, pid, sid):
    # remove the files of the older names, once the current ones are written
    for fname in get_old_output_paths(output_rootdir, pid, sid):
        if os.path.isfile(fname):
            os.remove(fname)


def read_existing_detail(detail_path):
    # read existing detail rows, keeping source=reviewer rows so a rerun
    # never discards reviews entered with the manual rawQC tool
    # a file's release names (or the short names of an older file) are
    # returned as the short keys
    if not os.path.isfile(detail_path):
        return []
    kept = []
    try:
        with open(detail_path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                row = detail_short_keys(row)
                if (row.get('source') or '').strip() == SOURCE_REVIEWER:
                    kept.append(row)
    except (IOError, OSError):
        print('{0}: WARNING: failed to read existing {1}'.format(
            APPNAME, detail_path))
    return kept


def detail_short_keys(row):
    # a detail row with its release names renamed to the short keys
    return dict((DETAIL_KEYS.get(k, k), v) for k, v in row.items())


def normalize_existing_row(row):
    # existing rows are read as strings with 'n/a' for missing values
    row = detail_short_keys(row)
    out = {}
    for col in DETAIL_COLUMNS:
        val = row.get(col)
        if is_empty(val):
            out[col] = None
        else:
            out[col] = str(val).strip()
    out['QC'] = norm_int(out['QC'])
    for col in QU_COLUMNS:
        out[col] = norm_float(out[col])
    out['duration'] = norm_float(out['duration'])
    return out


def write_detail(detail_path, rows):
    with open(detail_path, 'w') as f:
        f.write('\t'.join(DETAIL_NAMES.get(c, c)
                          for c in DETAIL_COLUMNS) + '\n')
        for row in rows:
            f.write('\t'.join(
                fmt_detail_value(c, row.get(c)) for c in DETAIL_COLUMNS) + '\n')


def write_summary(summary_path, rows):
    with open(summary_path, 'w') as f:
        f.write('\t'.join(SUMMARY_NAMES.get(c, c)
                          for c in SUMMARY_COLUMNS) + '\n')
        for row in rows:
            f.write('\t'.join(
                fmt_summary_value(c, row.get(c))
                for c in SUMMARY_COLUMNS) + '\n')


def write_status(status_path, pid, sid, errcode, errmsg):
    now = datetime.datetime.now()
    procdate = now.strftime('%Y-%m-%dT%H:%M:%S.') + \
        '{0:03d}'.format(now.microsecond // 1000)
    header = ['participant_id', 'session_id', 'procdate', 'errcode', 'errmsg']
    values = [pid, sid, procdate, str(errcode), errmsg or '']
    with open(status_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        f.write('\t'.join(values) + '\n')
    # pure python with no compiled app, so there is no app_version to record;
    # software_version identifies the code, since it ships in MMPS
    cmig_sing_utils.add_status_versions(status_path, app_name=APPNAME)


def write_session_output(pid, sid, results, reviewer_ids, output_rootdir,
                         verbose=0):
    # build detail rows from matched csv rows, merge with any existing
    # reviewer rows, recompute the summary, and write all three files
    #   (a session with output is here only with forceflag, or with its files
    #   under the older names, which are read and then removed)
    outdir, detail_path, summary_path, status_path = get_output_paths(
        output_rootdir, pid, sid)

    existing = [normalize_existing_row(r)
                for r in read_existing_detail(
                    get_existing_detail_path(output_rootdir, pid, sid))]

    detail_rows = []
    nscans_with_qc = 0
    for scan, csvrow, method in results:
        if csvrow is None:
            continue
        rows = build_detail_rows(pid, sid, scan, csvrow, reviewer_ids)
        if rows:
            nscans_with_qc += 1
        detail_rows.extend(rows)
    detail_rows.extend(existing)

    # group by scan for the summary, preserving scans.tsv order
    # scans with no csv row are included with nrev = 0, so the summary
    # lists every scan in the session
    order = []
    by_scan = {}
    uid_by_scan = {}
    for scan, csvrow, method in results:
        if csvrow is not None:
            stype = csvrow.get('SeriesType') or ''
        else:
            stype = SERIES_TYPE_MAP.get(
                (scan.get('type') or '').strip(),
                (scan.get('type') or '').strip())
        key = ((scan.get('filename') or '').strip(), stype)
        if key not in by_scan:
            by_scan[key] = []
            order.append(key)
        # the csv row is the same one build_detail_rows takes the uid from,
        # and it is in hand here even for a scan no reviewer opened
        if csvrow is not None and csvrow.get('SeriesInstanceUID'):
            uid_by_scan[key] = csvrow['SeriesInstanceUID']
        elif scan.get('series_uid'):
            uid_by_scan[key] = scan['series_uid']
    for row in detail_rows:
        key = (row.get('filename'), row.get('series_type'))
        if key not in by_scan:
            by_scan[key] = []
            order.append(key)
        by_scan[key].append(row)
        if key not in uid_by_scan and row.get('series_uid'):
            uid_by_scan[key] = row['series_uid']

    summary_rows = [compute_consensus(pid, sid, key, by_scan[key],
                                      uid_by_scan.get(key))
                    for key in order]

    os.makedirs(outdir, exist_ok=True)
    write_detail(detail_path, detail_rows)
    write_summary(summary_path, summary_rows)
    write_status(status_path, pid, sid, 0, '')
    remove_old_output(output_rootdir, pid, sid)
    if verbose:
        print('{0}: {1} {2}: wrote {3} detail rows ({4} preserved reviewer '
              'rows), {5} summary rows'.format(
                  APPNAME, pid, sid, len(detail_rows), len(existing),
                  len(summary_rows)))
    return len(detail_rows), nscans_with_qc


def check_prior_output(output_rootdir, pid, sid):
    # with forceflag=0, skip sessions that already have output
    outdir, detail_path, summary_path, status_path = get_output_paths(
        output_rootdir, pid, sid)
    return os.path.isfile(detail_path)

#############################################################################
# main
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Copy ABCD manual raw QC scores from DAIC merged QC '
                    'csv to BIDS derivatives.')
    parser.add_argument('--raw_rootdir', required=True,
                        help='root directory of BIDS rawdata')
    parser.add_argument('--output_rootdir', required=True,
                        help='root directory of BIDS derivatives output')
    parser.add_argument('--qcinfo_file', required=True,
                        help='csv file with merged manual QC info '
                             '(DAL_ABCD_QC_merged_pcqcinfo.csv)')
    parser.add_argument('--participant_id', nargs='+', default=None,
                        help='optional participant filter; one or more '
                             'ids (e.g. sub-XXXXXXXX)')
    parser.add_argument('--session_id', default=None,
                        help='optional session filter (e.g. ses-00A)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='display status messages')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='overwrite existing manual QC output files')
    return parser.parse_args()


def run(args):
    if not os.path.isdir(args.raw_rootdir):
        print('{0}: ERROR: raw_rootdir {1} not found'.format(
            APPNAME, args.raw_rootdir))
        return 1
    if not os.path.isfile(args.qcinfo_file):
        print('{0}: ERROR: qcinfo_file {1} not found'.format(
            APPNAME, args.qcinfo_file))
        return 1

    # load csv once; filter rows if participant filter given
    pguid_filter = None
    if args.participant_id:
        pguid_filter = set([pid_to_pguid(p) for p in args.participant_id])
    uid_index, time_index, reviewer_ids, qc_sessions = load_qcinfo(
        args.qcinfo_file, pguid_filter, args.verbose)
    if args.verbose:
        print('{0}: reviewer ids: {1}'.format(APPNAME, ' '.join(reviewer_ids)))

    # enumerate sessions from the csv, checking rawdata for scans.tsv
    sessions = get_sessions(args.raw_rootdir, qc_sessions,
                            args.participant_id, args.session_id,
                            args.verbose)
    if not sessions:
        print('{0}: ERROR: no sessions with manual QC and a scans.tsv '
              'found in {1}'.format(APPNAME, args.raw_rootdir))
        return 1
    print('{0}: {1} sessions to process'.format(APPNAME, len(sessions)))

    # skip sessions with existing output unless forceflag
    if not args.forceflag:
        nfound = len(sessions)
        sessions = [s for s in sessions
                    if not check_prior_output(args.output_rootdir,
                                              s[0], s[1])]
        nskip = nfound - len(sessions)
        if nskip:
            print('{0}: skipping {1} sessions with existing output '
                  '(use --forceflag 1 to overwrite)'.format(APPNAME, nskip))
        if not sessions:
            return 0

    # loop over sessions; match scans to csv rows and write output
    totals = {'sessions': 0, 'scans': 0, 'uid': 0, 'time': 0, 'unmatched': 0,
              'rows': 0, 'reviewed': 0}
    for pid, sid, scans_path in sessions:
        results, counts = match_session(
            pid, sid, scans_path, uid_index, time_index, args.verbose)
        totals['sessions'] += 1
        totals['scans'] += len(results)
        for key in ['uid', 'time', 'unmatched']:
            totals[key] += counts[key]
        nrows, nreviewed = write_session_output(
            pid, sid, results, reviewer_ids, args.output_rootdir,
            args.verbose)
        totals['rows'] += nrows
        totals['reviewed'] += nreviewed

    print('{0}: processed {1} sessions, {2} scans: {3} matched by uid, '
          '{4} matched by time, {5} unmatched'.format(
              APPNAME, totals['sessions'], totals['scans'],
              totals['uid'], totals['time'], totals['unmatched']))
    print('{0}: wrote {1} manual QC rows for {2} reviewed scans'.format(
        APPNAME, totals['rows'], totals['reviewed']))
    return 0


def main():
    sys.exit(run(parse_args()))


if __name__ == '__main__':
    main()
