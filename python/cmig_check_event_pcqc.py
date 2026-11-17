#!/usr/bin/env python3
"""
Session-level protocol compliance plus QC for BIDS-formatted ABCD MRI data.

Merges protocol compliance, automated raw QC and manual raw QC for one
session, then applies the event-level checks that abcd_check_event_pcqc
applies to a cohort-wide table.

Created:  09/10/2026 by Don Hagler
Prev mod: 09/21/2026 by Don Hagler
Prev mod: 09/24/2026 by Don Hagler
Prev mod: 09/25/2026 by Don Hagler
Prev mod: 10/08/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/09/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Last mod: 11/15/2026 by Don Hagler
"""

import argparse
import datetime
import glob
import os
import re
import sys

import pandas as pd

import cmig_sing_utils
import cmig_write_status as cmig_status
from cmig_pc_bids_abcd import SERIES_FM_COLUMN, decode_series_value

APPNAME = 'cmig_check_event_pcqc'

#############################################################################
# check parameters, from abcd_check_event_pcqc
#############################################################################

STYPES      = ['t1', 't2', 'dmri', 'rsfmri', 'mid', 'sst', 'nback']
SMRI_STYPES = ['t1', 't2']
DMRI_STYPES = ['dmri']
FMRI_STYPES = ['rsfmri', 'mid', 'sst', 'nback']
FMAP_STYPES = ['dmri', 'rsfmri', 'mid', 'sst', 'nback']

SERIES_TYPES = ['T1', 'T2', 'dMRI', 'rsfMRI',
                'fMRI_MID_task', 'fMRI_SST_task', 'fMRI_nBack_task']
MIN_NSCANS = [1, 1, 1, 3, 2, 2, 2]
MAX_NSCANS = [3, 3, 6, 12, 8, 8, 8]

# the autoQC measures, named as the rawqc TSVs name them, below the stem
COMMON_AUTOQC_TLIST = ['vol__vox_mean', 'vol__brain_mean',
                       'intensity__brain_mean', 'intensity__brain_std']
SMRI_AUTOQC_TLIST = []
DMRI_AUTOQC_TLIST = [
    'rep_count', 'mot_mean', 'mot__nody_mean',
    'fa__brain_mean', 'fa__brain_med', 'fa__brain_std',
    'md__brain_mean', 'md__brain_med', 'md__brain_std',
    'intensity__b0_mean', 'intensity__b0_med', 'intensity__b0_std',
    'dterr_mean', 'dterr_med', 'dterr_std',
]
FMRI_AUTOQC_TLIST = [
    'rep_count', 'tr', 'mot_mean', 'mot__nody_mean',
    'mot__filt_mean', 'mot__filt__nody_mean',
    'dtr__xdir_fwhm', 'dtr__ydir_fwhm', 'dtr__zdir_fwhm',
    'snr_mean', 'snr_med', 'snr_std',
]
FMRI_SUBTHRESH_TLIST = ['fdsub__02_sum', 'fdsub__02__nody_sum',
                        'fdsub__02__filt_sum', 'fdsub__02__filt__nody_sum']
# the FOV cutoffs of the dMRI and fMRI scans, taken as a max over the ok scans
FOV_MAX_TLIST = ['fov__cutoff__dorsal_max', 'fov__cutoff__ventral_max']

# the artifacts each modality scores. Motion is rated for the 3D structural
# scans and estimated from the data elsewhere, where it shows as a line
# artifact; an FA map exists only for dMRI
QU_ISSUES_SMRI = ['mo']
QU_ISSUES_DMRI = ['fov', 'sus', 'line', 'fa', 'gh']
# moire is scored for fMRI only: dMRI carries 18 reviews of it, too few to
#   rate, so the two dMRI moire counts stay out of the release
QU_ISSUES_FMRI = ['fov', 'sus', 'line', 'gh', 'moire']
QU_ISSUES_ALL = ['mo', 'fov', 'sus', 'line', 'fa', 'gh', 'moire']

QU_COLUMN = {'mo': 'QU_motion', 'sus': 'QU_sus', 'fov': 'QU_cutoff',
             'line': 'QU_line', 'fa': 'QU_fa', 'gh': 'QU_ghost',
             'moire': 'QU_moire'}
QU_COLUMNS = ['QU_motion', 'QU_sus', 'QU_cutoff', 'QU_line', 'QU_fa',
              'QU_ghost', 'QU_moire']

# the manual raw QC scans TSV's release names (mr_y_qc__raw__man__series),
#   by the short keys used here; a file written before the rename has the
#   short names, which are read where the release name is absent
MANUAL_TABLE = 'mr_y_qc__raw__man__series'
MANUAL_NAMES = {
    'QC': MANUAL_TABLE + '_score',
    'QU_motion': MANUAL_TABLE + '__mot_score',
    'QU_sus': MANUAL_TABLE + '__susc_score',
    'QU_cutoff': MANUAL_TABLE + '__cutoff_score',
    'QU_line': MANUAL_TABLE + '__line_score',
    'QU_fa': MANUAL_TABLE + '__famap_score',
    'QU_ghost': MANUAL_TABLE + '__ghost_score',
    'QU_moire': MANUAL_TABLE + '__moire_score',
    'qu_source': MANUAL_TABLE + '__score_source',
    'nrev': MANUAL_TABLE + '__rev_count',
    'revdisp': MANUAL_TABLE + '__rev__disp_indicator',
    'series_uid': MANUAL_TABLE + '_uid',
}

# the per-reviewer QC TSV's release names (mr_y_qc__raw__man__series__review),
#   read like MANUAL_NAMES
REVIEW_TABLE = 'mr_y_qc__raw__man__series__review'
REVIEW_NAMES = {
    'notes': REVIEW_TABLE + '_notes',
    'series_uid': REVIEW_TABLE + '_uid',
}

# a scan joins an artifact's session total when the reviewers judged it
# severe. A reviewer who passed a scan cannot have scored the artifact 3, so
# the threshold is never reached on a scan whose consensus QC is 1
QU_SEVERE = 2.5

SERIES_TAGS = ['qc_score', 'qc_notes', 'pc_score', 'pc_bad', 'pc_notes', 'complete']
FMAP_TAGS   = ['fm_missing', 'fm_wronglocus', 'fm_invalid',
               'fm_qc_score', 'fm_qc_notes',
               'fm_pc_score', 'fm_pc_bad', 'fm_pc_notes', 'fm_complete']
INFO_TAGS        = ['StudyInstanceUID', 'SeriesInstanceUID', 'StudyDate', 'SeriesTime']
INFO_TAGS_RENAME = ['suid', 'seuid', 'studydate', 'seriestime']

#############################################################################
# input assembly
#############################################################################

# the rawqc column stem each series type is measured under
STEM_BY_STYPE = {
    't1': 'smri__t1', 't2': 'smri__t2', 'dmri': 'dmri', 'rsfmri': 'rsfmri',
    'mid': 'tfmri__mid', 'sst': 'tfmri__sst', 'nback': 'tfmri__nback',
}

RAWQC_GLOBS = [
    os.path.join('anat', '*_desc-T1w_rawqc.tsv'),
    os.path.join('anat', '*_desc-T2w_rawqc.tsv'),
    os.path.join('dwi',  '*_desc-dwi_rawqc.tsv'),
    os.path.join('func', '*_desc-bold_rawqc.tsv'),
]

# the series type each rawqc file holds; the per-series files carry one
# modality each and name their measures under a neutral stem, so the type
# comes from the file, and for fMRI from the task entity of the row
STYPE_BY_GLOB = {
    RAWQC_GLOBS[0]: 't1',
    RAWQC_GLOBS[1]: 't2',
    RAWQC_GLOBS[2]: 'dmri',
    RAWQC_GLOBS[3]: None,
}

STYPE_BY_TASK = {
    'rest': 'rsfmri', 'mid': 'mid', 'sst': 'sst', 'nback': 'nback',
}

NA_STRINGS = ('', 'n/a', 'N/A', 'NaN', 'nan')


def _is_empty(val):
    """True for the values the merged inputs use to mean absent."""
    if val is None:
        return True
    if isinstance(val, float):
        return val != val
    if isinstance(val, str):
        return val.strip() in NA_STRINGS
    if isinstance(val, list):
        return len(val) == 0
    return False


def _num(val):
    """Return val as a number, or None where it is absent or unparsable."""
    if _is_empty(val):
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def read_tsv(fpath):
    """Read a TSV as strings, leaving absent values as written."""
    return pd.read_csv(fpath, sep='\t', dtype=str, keep_default_na=False)


def read_manual_qc(manualqc_rootdir, pid, sid):
    """Return manual QC keyed by series uid, and by filename where absent.

    The consensus comes from the scans TSV, which cmig_copy_manual_rawqc_bids_abcd
    computes; the notes come from the per-reviewer QC TSV, one row per scan and
    reviewer, which the merged table carried as one notes column per reviewer.

    Both keys are returned for every row. series_uid identifies the scan
    itself and survives a rename or a reconversion, where filename does not,
    so the caller tries it first; manual QC written before series_uid was
    added carries none, and those rows are reachable only by filename.
    """
    sessdir = os.path.join(manualqc_rootdir, pid, sid)
    consensus = {}
    notes = {}
    reviewers = []

    fpath = manualqc_fpath(sessdir, pid, sid, 'desc-series_manualrawqc',
                           'desc-manualrawqc_scans')
    if os.path.isfile(fpath):
        df = read_tsv(fpath)
        df = df.rename(columns=dict((v, k) for k, v in MANUAL_NAMES.items()
                                    if v in df.columns))
        for _, row in df.iterrows():
            fname = row.get('filename')
            if _is_empty(fname):
                continue
            con = {'QC': _num(row.get('QC')),
                   'nrev': _num(row.get('nrev')),
                   'revdisp': _num(row.get('revdisp')),
                   'qu_source': row.get('qu_source')}
            # a column the modality does not score stays absent rather than 0
            for col in QU_COLUMNS:
                con[col] = _num(row.get(col))
            consensus[fname] = con
            uid = row.get('series_uid')
            if not _is_empty(uid):
                consensus[uid] = con

    fpath = manualqc_fpath(sessdir, pid, sid, 'desc-reviews_manualrawqc',
                           'desc-manualrawqc_qc')
    if os.path.isfile(fpath):
        df = read_tsv(fpath)
        df = df.rename(columns=dict((v, k) for k, v in REVIEW_NAMES.items()
                                    if v in df.columns))
        for _, row in df.iterrows():
            fname = row.get('filename')
            revid = row.get('reviewer_id')
            if _is_empty(fname) or _is_empty(revid):
                continue
            if revid not in reviewers:
                reviewers.append(revid)
            if not _is_empty(row.get('notes')):
                notes.setdefault(fname, {})[revid] = row['notes']
                uid = row.get('series_uid')
                if not _is_empty(uid):
                    notes.setdefault(uid, {})[revid] = row['notes']

    return consensus, notes, sorted(reviewers)


def manualqc_fpath(sessdir, pid, sid, name, old_name):
    """A manual raw QC file, under its older name where the current one is
    missing (the older names are converted when the session is next saved)."""
    fpath = os.path.join(sessdir, '{0}_{1}_{2}.tsv'.format(pid, sid, name))
    old_fpath = os.path.join(sessdir, '{0}_{1}_{2}.tsv'.format(pid, sid, old_name))
    if not os.path.isfile(fpath) and os.path.isfile(old_fpath):
        return old_fpath
    return fpath


def stype_from_task(fname):
    """The fMRI series type of a rawqc row, from the task entity."""
    for part in os.path.basename(fname).split('_'):
        if part.startswith('task-'):
            return STYPE_BY_TASK.get(part[5:])
    return None


def read_rawqc(rawqc_rootdir, pid, sid):
    """Return automated QC measures keyed by file name, under check_event's names."""
    sessdir = os.path.join(rawqc_rootdir, pid, sid)
    measures = {}
    for pattern in RAWQC_GLOBS:
        matches = sorted(glob.glob(os.path.join(sessdir, pattern)))
        for fpath in matches:
            df = read_tsv(fpath)
            for _, row in df.iterrows():
                fname = row.get('filename')
                if _is_empty(fname):
                    continue
                stype = STYPE_BY_GLOB.get(pattern)
                if stype is None:
                    stype = stype_from_task(fname)
                if stype is None:
                    continue
                tlist = autoqc_tlist(stype) + fov_max_tlist(stype)
                if stype in FMRI_STYPES:
                    tlist = tlist + FMRI_SUBTHRESH_TLIST
                vals = {}
                for tag in tlist:
                    col = 'mr_y_qc__raw__series__{0}'.format(tag)
                    if col in df.columns and not _is_empty(row[col]):
                        vals[tag] = _num(row[col])
                if vals:
                    measures[fname] = vals
    return measures


# the field map columns of pc's series TSV, flags and text
PC_FMAP_FLAGS = ['fm_missing', 'fm_wronglocus', 'fm_invalid', 'fm_pc_score',
                 'fm_pc_bad', 'fm_complete']
PC_FMAP_TEXT = ['fm_pa_filename', 'fm_ap_filename', 'fm_pc_notes']


def build_series_table(pc_df, consensus, notes, reviewers, measures):
    """Assemble the one-row-per-series table check_event reads.

    The merged table cmig replaces carried protocol compliance, manual QC and
    automated QC in one row per series; this rebuilds that row from the three
    cmig sources, joined on the series uid, under the names check_event uses.
    """
    rows = []
    for _, pc in pc_df.iterrows():
        fname = pc.get('filename', '')
        row = {
            'SeriesType':        decode_series_value(
                'mr_y_qc__raw__series_type',
                pc.get('mr_y_qc__raw__series_type', '')),
            'Compliant':         decode_series_value(
                'mr_y_qc__raw__series__protcom_indicator',
                pc.get('mr_y_qc__raw__series__protcom_indicator', '')),
            'AdditionalInfo':    pc.get('mr_y_qc__raw__series_notes', ''),
            'Manufacturer':      pc.get('mr_y_qc__raw__series__dev_manufact', ''),
            'SeriesNumber':      pc.get('mr_y_qc__raw__series_num', ''),
            'ImageType':         pc.get('mr_y_qc__raw__series__img_type', ''),
            'PulseSequenceName': pc.get('mr_y_qc__raw__series_psn', ''),
            'StudyDate':         pc.get('mr_y_qc__raw__series__study_dt', ''),
            'SeriesTime':        pc.get('mr_y_qc__raw__series__series_t', ''),
            'StudyInstanceUID':  pc.get('mr_y_qc__raw__series__study_uid', ''),
            'SeriesInstanceUID': pc.get('mr_y_qc__raw__series__series_uid', ''),
            'filename':          fname,
        }
        # check_event tests Completed for truth, so the pc codes become 1 and 0
        completed = decode_series_value(
            'mr_y_qc__raw__series__complete_indicator',
            pc.get('mr_y_qc__raw__series__complete_indicator', ''))
        row['Completed'] = 1 if completed == 'Yes' else 0 if completed == 'No' else None

        # series_uid first, the scan's own identifier; filename is what the
        #   manual QC carried before it and is kept as the fallback
        suid = pc.get('mr_y_qc__raw__series__series_uid', '')
        con = consensus.get(suid) if not _is_empty(suid) else None
        if con is None:
            con = consensus.get(fname, {})
        row['QC'] = con.get('QC')
        row['nrev'] = con.get('nrev')
        row['revdisp'] = con.get('revdisp')
        row['qu_source'] = con.get('qu_source')
        for col in QU_COLUMNS:
            row[col] = con.get(col)

        rev_notes = notes.get(suid) if not _is_empty(suid) else None
        if rev_notes is None:
            rev_notes = notes.get(fname, {})
        for revid in reviewers:
            row['notes_{0}'.format(revid)] = rev_notes.get(revid, '')

        # the repetitions and TR come from pc, which reads them from the
        #   sidecar; rawqc writes neither, so its measures cannot supply them
        row['rep_count'] = _num(pc.get('mr_y_qc__raw__series__frame_count', ''))
        row['tr'] = _num(pc.get('mr_y_qc__raw__series_tr', ''))

        # the field map pair and its checks come from pc, which matches it
        for tag in PC_FMAP_FLAGS:
            val = _num(pc.get(SERIES_FM_COLUMN[tag], ''))
            row[tag] = None if val is None else int(val)
        for tag in PC_FMAP_TEXT:
            row[tag] = pc.get(SERIES_FM_COLUMN[tag], '')

        row.update(measures.get(fname, {}))
        rows.append(row)
    return pd.DataFrame(rows)

#############################################################################
# the checks, ported from abcd_check_event_pcqc
#############################################################################

def get_relevant_series(df):
    """Select ABCD-relevant series rows."""
    keep = df['Compliant'].apply(lambda v: str(v).strip() in ('Yes', 'No', 'Bad'))
    keep = keep & ~df['SeriesType'].astype(str).str.contains('Undefined', na=False)
    return df[keep].copy()


def convert_seriestype(df, seriestype_old, seriestype_new):
    """Rename SeriesType values (list-to-list mapping)."""
    for old, new in zip(seriestype_old, seriestype_new):
        df.loc[df['SeriesType'] == old, 'SeriesType'] = new
    return df


def find_qc_notes(df):
    """Return the per-reviewer notes columns, in column order."""
    return [c for c in df.columns if re.match(r'notes_\w+', c)]


def get_row_val(row, col, default=None):
    """Value from a row, or default where the column or the value is absent."""
    if col not in row.index:
        return default
    val = row[col]
    if _is_empty(val):
        return default
    return val


def get_qc_notes(row, qc_notes_tags, infix=None):
    """Assemble the QC notes string from the per-reviewer notes columns.

    Each note is labelled with the reviewer's initials. The cohort table
    labelled them R1, R2 and so on by position in a fixed column per reviewer,
    which one session cannot know.
    """
    notes = ''
    for tag in qc_notes_tags:
        note = get_row_val(row, tag, '')
        if note and str(note).strip():
            revid = tag[len('notes_'):]
            infix_str = str(infix) if infix is not None else ''
            notes += '{0}{1}: {2} '.format(revid, infix_str, note)
    return notes.strip()


def init_scan_info():
    """Initialize the scan-level counters."""
    scan_info = {
        'total_ser': 0, 'total_passpc': 0, 'total_passqc': 0,
        'good_ser': 0, 'ok_ser': 0, 'comp_ser': 0,
        'ser_failpc': 0, 'ser_failqc': 0, 'ser_incomp': 0, 'qcd': 0,
        'total_passpc_wofm': 0, 'total_passqc_wofm': 0, 'good_ser_wofm': 0,
        'fm_passpc': 0, 'fm_passqc': 0, 'fm_missing': 0, 'fm_wronglocus': 0,
        'fm_invalid': 0, 'fm_failpc': 0, 'fm_failqc': 0, 'fm_incomp': 0,
    }
    for issue in QU_ISSUES_ALL:
        scan_info[issue] = 0
        scan_info['fm_' + issue] = 0
    all_autoqc = list(dict.fromkeys(
        COMMON_AUTOQC_TLIST + SMRI_AUTOQC_TLIST + DMRI_AUTOQC_TLIST + FMRI_AUTOQC_TLIST))
    for t in all_autoqc:
        m = t.lower()
        scan_info[m] = 0
        scan_info['ok_' + m] = 0
        scan_info['comp_' + m] = 0
    for t in FMRI_SUBTHRESH_TLIST:
        m = t.lower()
        for prefix in ('', 'wofm_', 'good_', 'ok_', 'comp_'):
            scan_info[prefix + m] = 0
    return scan_info


def qu_issues(stype):
    """The artifacts a series type's reviewers score."""
    if stype in SMRI_STYPES:
        return QU_ISSUES_SMRI
    if stype in DMRI_STYPES:
        return QU_ISSUES_DMRI
    if stype in FMRI_STYPES:
        return QU_ISSUES_FMRI
    return []


def autoqc_tlist(stype):
    """The autoQC measures check_event reads for a series type."""
    if stype in SMRI_STYPES:
        return COMMON_AUTOQC_TLIST + SMRI_AUTOQC_TLIST
    if stype in DMRI_STYPES:
        return COMMON_AUTOQC_TLIST + DMRI_AUTOQC_TLIST
    if stype in FMRI_STYPES:
        return COMMON_AUTOQC_TLIST + FMRI_AUTOQC_TLIST
    return []


def fov_max_tlist(stype):
    """The FOV cutoff measures check_event takes a max of for a series type."""
    if stype in DMRI_STYPES or stype in FMRI_STYPES:
        return FOV_MAX_TLIST
    return []


def init_event_info(pid, sid, series_rows):
    """Create the event row with its column order fixed."""
    event_info = {'participant_id': pid, 'session_id': sid}

    for i, stype in enumerate(STYPES):
        for j in range(1, MAX_NSCANS[i] + 1):
            for k, tag in enumerate(INFO_TAGS):
                event_info['iqc_{0}_{1}_{2}'.format(stype, j, INFO_TAGS_RENAME[k])] = None
            if stype in FMAP_STYPES:
                for tag in FMAP_TAGS:
                    event_info['iqc_{0}_{1}_{2}'.format(stype, j, tag)] = None
            for tag in SERIES_TAGS:
                event_info['iqc_{0}_{1}_{2}'.format(stype, j, tag)] = None
            tlist = autoqc_tlist(stype)
            if stype in FMRI_STYPES:
                tlist = tlist + FMRI_SUBTHRESH_TLIST
            for f in tlist:
                event_info['iqc_{0}_{1}_{2}'.format(stype, j, f.lower())] = None

    init_event_aggregates(event_info)
    return event_info


def init_event_aggregates(event_info):
    """Create the across-scan columns, in the order the checks write them.

    A series type with no scans is never counted, so its aggregates would
    otherwise not appear at all and the column set would vary from session to
    session. One session per file makes that the reader's problem rather than
    a concatenation's, so every column exists here and stays empty.
    """
    for stype in STYPES:
        for tag in ['total_ser', 'total_passpc', 'ser_incomp', 'ok_ser',
                    'comp_ser', 'ser_pc_issues', 'total_passqc']:
            event_info['iqc_{0}_{1}'.format(stype, tag)] = None
        if stype in FMAP_STYPES:
            event_info['iqc_{0}_good_ser_wofm'.format(stype)] = None
        event_info['iqc_{0}_good_ser'.format(stype)] = None
        for issue in qu_issues(stype):
            event_info['iqc_{0}_{1}_qc'.format(stype, issue)] = None
        event_info['iqc_{0}_ser_qcs'.format(stype)] = None
        event_info['iqc_{0}_bad_ser'.format(stype)] = None

        if stype in FMAP_STYPES:
            for tag in ['fm_missing', 'fm_wronglocus', 'fm_invalid',
                        'fm_pc_issues', 'fm_incomp', 'fm_qcs']:
                event_info['iqc_{0}_{1}'.format(stype, tag)] = None
            for issue in qu_issues(stype):
                event_info['iqc_{0}_fm_{1}_qc'.format(stype, issue)] = None

        for t in autoqc_tlist(stype):
            m_tag = t.lower()
            for prefix in ('ok', 'comp'):
                event_info['iqc_{0}_{1}_{2}'.format(stype, prefix, m_tag)] = None
        if stype in FMRI_STYPES:
            for t in FMRI_SUBTHRESH_TLIST:
                m_tag = t.lower()
                if stype in FMAP_STYPES:
                    event_info['iqc_{0}_ok_{1}'.format(stype, m_tag)] = None
                    event_info['iqc_{0}_comp_{1}'.format(stype, m_tag)] = None
        for t in fov_max_tlist(stype):
            event_info['iqc_{0}_ok_{1}'.format(stype, t)] = None

    for stem, _ in FMAP_PTYPES:
        event_info['iqc_{0}_fm_ok'.format(stem)] = None
        event_info['iqc_{0}_fm_passpc_comp'.format(stem)] = None
    for tag in ['eventcomplete', 'eventcomplete_ok']:
        event_info['iqc_{0}'.format(tag)] = None
    return event_info


def get_scan_info_from_event(event_info, stype, j):
    return (event_info.get('iqc_{0}_{1}_qc_score'.format(stype, j)),
            event_info.get('iqc_{0}_{1}_pc_score'.format(stype, j)),
            event_info.get('iqc_{0}_{1}_pc_bad'.format(stype, j)),
            event_info.get('iqc_{0}_{1}_complete'.format(stype, j)))


def get_fmap_info(event_info, stype, j):
    if stype in FMAP_STYPES:
        return (1,
                event_info.get('iqc_{0}_{1}_fm_missing'.format(stype, j)),
                event_info.get('iqc_{0}_{1}_fm_wronglocus'.format(stype, j)),
                event_info.get('iqc_{0}_{1}_fm_invalid'.format(stype, j)),
                event_info.get('iqc_{0}_{1}_fm_qc_score'.format(stype, j)),
                event_info.get('iqc_{0}_{1}_fm_pc_score'.format(stype, j)),
                event_info.get('iqc_{0}_{1}_fm_pc_bad'.format(stype, j)),
                event_info.get('iqc_{0}_{1}_fm_complete'.format(stype, j)))
    return 0, 0, 0, 0, None, None, None, None


def count_qu_issues(row, counter, stype, prefix=''):
    """Count the artifacts the reviewers judged severe on one scan.

    The severity comes from the consensus QU scores, which
    cmig_copy_manual_rawqc_bids_abcd works out per scan; the notes are not
    read. A column the modality does not score is absent, which is not the
    same as a score of zero.
    """
    for issue in qu_issues(stype):
        val = _num(get_row_val(row, QU_COLUMN[issue], None))
        if val is not None and val >= QU_SEVERE:
            counter[prefix + issue] += 1


def count_qu_issues_pair(fwd_row, rev_row, counter, stype, prefix=''):
    """The same over a field map pair, either member carrying the artifact."""
    for issue in qu_issues(stype):
        vals = [_num(get_row_val(r, QU_COLUMN[issue], None))
                for r in (fwd_row, rev_row)]
        vals = [v for v in vals if v is not None]
        if vals and max(vals) >= QU_SEVERE:
            counter[prefix + issue] += 1


def check_series_pcqc(event_info, scan_info, series_rows, stype, j, qc_notes_tags):
    """Compliance and manual QC for scan j (1-indexed) of stype."""
    row = series_rows.iloc[j - 1]

    compliant = str(get_row_val(row, 'Compliant', '')).strip().lower()
    event_info['iqc_{0}_{1}_pc_score'.format(stype, j)] = 1 if compliant == 'yes' else 0
    event_info['iqc_{0}_{1}_pc_bad'.format(stype, j)] = 1 if compliant == 'bad' else 0
    event_info['iqc_{0}_{1}_pc_notes'.format(stype, j)] = get_row_val(row, 'AdditionalInfo', '')

    qc_score = get_row_val(row, 'QC', None)
    event_info['iqc_{0}_{1}_qc_score'.format(stype, j)] = qc_score

    qc_notes = get_qc_notes(row, qc_notes_tags)
    event_info['iqc_{0}_{1}_qc_notes'.format(stype, j)] = qc_notes

    count_qu_issues(row, scan_info, stype)

    event_info['iqc_{0}_{1}_complete'.format(stype, j)] = get_row_val(row, 'Completed', None)
    return event_info, scan_info


def acquisition_slots(all_rows):
    """Number each row by its acquisition, rows sharing a time being one.

    A pepolar sequence acquires both phase encoding directions in one series,
    which the converter splits into two files, so those two rows carry the same
    StudyDate and SeriesTime. Grouping by time puts them in one slot, where
    counting rows would make them two.
    """
    slots, last, n = [], None, -1
    for _, row in all_rows.iterrows():
        key = (row.get('StudyDate', ''), row.get('SeriesTime', ''))
        if key != last:
            n += 1
            last = key
        slots.append(n)
    return slots


def check_series_fmap(event_info, scan_info, all_rows, series_rows, stype, j, qc_notes_tags):
    """Field map variables for scan j of stype.

    cmig_pc_bids_abcd matches the field map pair to the scan and judges its
    placement and compliance; this reads that and adds the pair's manual QC.
    """
    if stype not in FMAP_STYPES:
        return event_info, scan_info

    row = series_rows.iloc[j - 1]
    fm_missing = get_row_val(row, 'fm_missing', 1)
    fm_wronglocus = get_row_val(row, 'fm_wronglocus', 1)
    fm_invalid = get_row_val(row, 'fm_invalid', 1)
    fm_pc_score = get_row_val(row, 'fm_pc_score', None)
    fm_pc_bad = get_row_val(row, 'fm_pc_bad', None)
    fm_pc_notes = get_row_val(row, 'fm_pc_notes', None)
    fm_completed = get_row_val(row, 'fm_complete', None)
    fm_qc_score = fm_qc_notes = None

    fnames = all_rows['filename'].tolist()
    pa = get_row_val(row, 'fm_pa_filename', '')
    ap = get_row_val(row, 'fm_ap_filename', '')
    if pa in fnames and ap in fnames:
        fwd_row = all_rows.iloc[fnames.index(pa)]
        rev_row = all_rows.iloc[fnames.index(ap)]

        fwd_qc = get_row_val(fwd_row, 'QC', None)
        rev_qc = get_row_val(rev_row, 'QC', None)
        if fwd_qc is not None and rev_qc is not None:
            fm_qc_score = fwd_qc * rev_qc

        notes1 = get_qc_notes(fwd_row, qc_notes_tags, 'FM1')
        notes2 = get_qc_notes(rev_row, qc_notes_tags, 'FM2')
        if notes1 and notes2:
            fm_qc_notes = notes1 + ' - ' + notes2
        else:
            fm_qc_notes = notes1 or notes2 or None

        count_qu_issues_pair(fwd_row, rev_row, scan_info, stype, 'fm_')

    event_info['iqc_{0}_{1}_fm_missing'.format(stype, j)]    = fm_missing
    event_info['iqc_{0}_{1}_fm_wronglocus'.format(stype, j)] = fm_wronglocus
    event_info['iqc_{0}_{1}_fm_invalid'.format(stype, j)]    = fm_invalid
    event_info['iqc_{0}_{1}_fm_qc_score'.format(stype, j)]   = fm_qc_score
    event_info['iqc_{0}_{1}_fm_qc_notes'.format(stype, j)]   = fm_qc_notes
    event_info['iqc_{0}_{1}_fm_pc_score'.format(stype, j)]   = fm_pc_score
    event_info['iqc_{0}_{1}_fm_pc_bad'.format(stype, j)]     = fm_pc_bad
    event_info['iqc_{0}_{1}_fm_pc_notes'.format(stype, j)]   = fm_pc_notes
    event_info['iqc_{0}_{1}_fm_complete'.format(stype, j)]   = fm_completed

    return event_info, scan_info


def check_series_info(event_info, series_rows, stype, j):
    """Copy the series identifiers into the event row."""
    row = series_rows.iloc[j - 1]
    for i, i_tag in enumerate(INFO_TAGS):
        event_info['iqc_{0}_{1}_{2}'.format(stype, j, INFO_TAGS_RENAME[i])] = \
            get_row_val(row, i_tag)
    return event_info


def check_series_autoqc_various(event_info, scan_info, series_rows, stype, j):
    """Accumulate the autoQC measures for scan j."""
    row = series_rows.iloc[j - 1]
    qc_score, pc_score, pc_bad, completed = get_scan_info_from_event(event_info, stype, j)
    if _is_empty(completed):
        completed = 0
    if _is_empty(pc_bad):
        pc_bad = 0

    for q_tag in autoqc_tlist(stype):
        m_tag = q_tag.lower()
        qc_val = get_row_val(row, q_tag, None)
        event_info['iqc_{0}_{1}_{2}'.format(stype, j, m_tag)] = qc_val
        if not _is_empty(qc_val):
            scan_info[m_tag] = scan_info.get(m_tag, 0) + qc_val
            scan_info['n_' + m_tag] = scan_info.get('n_' + m_tag, 0) + 1
            if (_is_empty(qc_score) or qc_score == 1) and completed and not pc_bad:
                scan_info['ok_' + m_tag] = scan_info.get('ok_' + m_tag, 0) + qc_val
                scan_info['n_ok_' + m_tag] = scan_info.get('n_ok_' + m_tag, 0) + 1
            if completed and not pc_bad and pc_score:
                scan_info['comp_' + m_tag] = scan_info.get('comp_' + m_tag, 0) + qc_val
                scan_info['n_comp_' + m_tag] = scan_info.get('n_comp_' + m_tag, 0) + 1
    return event_info, scan_info


def check_series_autoqc_subthresh(event_info, scan_info, series_rows, stype, j):
    """Accumulate the fMRI sub-threshold motion counts for scan j."""
    if stype not in FMRI_STYPES:
        return event_info, scan_info

    row = series_rows.iloc[j - 1]
    qc_score, pc_score, pc_bad, completed = get_scan_info_from_event(event_info, stype, j)
    if _is_empty(completed):
        completed = 0
    (fm_flag, fm_missing, fm_wronglocus, fm_invalid,
     fm_qc_score, fm_pc_score, fm_pc_bad, fm_completed) = get_fmap_info(event_info, stype, j)
    if _is_empty(fm_completed):
        fm_completed = 0

    for t in FMRI_SUBTHRESH_TLIST:
        m_tag = t.lower()
        qc_val = get_row_val(row, t, None)
        event_info['iqc_{0}_{1}_{2}'.format(stype, j, m_tag)] = qc_val
        if _is_empty(qc_val):
            continue
        scan_info[m_tag] = scan_info.get(m_tag, 0) + qc_val
        scan_info['n_' + m_tag] = scan_info.get('n_' + m_tag, 0) + 1
        if (_is_empty(qc_score) or qc_score == 1) and completed and not pc_bad:
            scan_info['ok_' + m_tag] = scan_info.get('ok_' + m_tag, 0) + qc_val
            scan_info['n_ok_' + m_tag] = scan_info.get('n_ok_' + m_tag, 0) + 1
        if completed and not pc_bad and pc_score:
            scan_info['comp_' + m_tag] = scan_info.get('comp_' + m_tag, 0) + qc_val
            scan_info['n_comp_' + m_tag] = scan_info.get('n_comp_' + m_tag, 0) + 1
        if fm_flag:
            if pc_score and (not _is_empty(qc_score) and qc_score == 1) and completed:
                scan_info['wofm_' + m_tag] = scan_info.get('wofm_' + m_tag, 0) + qc_val
                if (not fm_missing and not fm_invalid and fm_pc_score and
                        not _is_empty(fm_qc_score) and fm_qc_score == 1 and fm_completed):
                    scan_info['good_' + m_tag] = scan_info.get('good_' + m_tag, 0) + qc_val
    return event_info, scan_info


def check_series_autoqc_max(event_info, scan_info, series_rows, stype, j):
    """Fold the FOV cutoffs of scan j into their max over the ok scans."""
    row = series_rows.iloc[j - 1]
    qc_score, pc_score, pc_bad, completed = get_scan_info_from_event(event_info, stype, j)
    if _is_empty(completed):
        completed = 0
    if _is_empty(pc_bad):
        pc_bad = 0
    if not ((_is_empty(qc_score) or qc_score == 1) and completed and not pc_bad):
        return event_info, scan_info
    for t in fov_max_tlist(stype):
        qc_val = get_row_val(row, t, None)
        if _is_empty(qc_val):
            continue
        # absent until an ok scan carries the measure, so none leaves it empty
        prev = scan_info.get('ok_' + t)
        scan_info['ok_' + t] = qc_val if prev is None else max(prev, qc_val)
    return event_info, scan_info


def count_series(event_info, scan_info, stype, j):
    """Update the scan counters after processing scan j."""
    qc_score, pc_score, pc_bad, completed = get_scan_info_from_event(event_info, stype, j)
    scan_info['total_ser'] += 1
    if not pc_score:
        scan_info['ser_failpc'] += 1
    if _is_empty(pc_score):
        pc_score = 0
    # an unreviewed scan counts as neither passed nor failed, where the MATLAB
    # counted an absent score as a failure
    if not _is_empty(qc_score) and qc_score != 1:
        scan_info['ser_failqc'] += 1
    if not completed:
        scan_info['ser_incomp'] += 1
    if _is_empty(completed):
        completed = 0

    (fm_flag, fm_missing, fm_wronglocus, fm_invalid,
     fm_qc_score, fm_pc_score, fm_pc_bad, fm_completed) = get_fmap_info(event_info, stype, j)

    if not fm_missing:
        if fm_wronglocus:
            scan_info['fm_wronglocus'] += 1
        if fm_invalid:
            scan_info['fm_invalid'] += 1
        if not fm_pc_score:
            scan_info['fm_failpc'] += 1
        if _is_empty(fm_pc_score):
            fm_pc_score = 0
        if fm_pc_score:
            scan_info['fm_passpc'] += 1
        if not _is_empty(fm_qc_score) and fm_qc_score != 1:
            scan_info['fm_failqc'] += 1
        if fm_qc_score:
            scan_info['fm_passqc'] += 1
        if not fm_completed:
            scan_info['fm_incomp'] += 1
        if _is_empty(fm_completed):
            fm_completed = 0
    else:
        scan_info['fm_missing'] += 1

    if not fm_flag:
        fm_pc_score = fm_qc_score = fm_completed = 1
    elif fm_missing:
        fm_pc_score = fm_qc_score = fm_completed = 0

    if pc_score and fm_pc_score:
        scan_info['total_passpc'] += 1
    if (not _is_empty(qc_score) and qc_score == 1) and \
       (not _is_empty(fm_qc_score) and fm_qc_score == 1):
        scan_info['total_passqc'] += 1
    if (pc_score and (not _is_empty(qc_score) and qc_score == 1) and completed and
            fm_pc_score and (not _is_empty(fm_qc_score) and fm_qc_score == 1) and fm_completed):
        scan_info['good_ser'] += 1
    if (_is_empty(qc_score) or qc_score == 1) and completed and not pc_bad:
        scan_info['ok_ser'] += 1
    if completed and not pc_bad and pc_score:
        scan_info['comp_ser'] += 1

    if fm_flag:
        if pc_score:
            scan_info['total_passpc_wofm'] += 1
        if not _is_empty(qc_score) and qc_score == 1:
            scan_info['total_passqc_wofm'] += 1
        if pc_score and (not _is_empty(qc_score) and qc_score == 1) and completed:
            scan_info['good_ser_wofm'] += 1
    return scan_info


def check_series_scans(event_info, all_rows, stype, i_type, qc_notes_tags):
    """Process every scan of one series type."""
    type_rows = all_rows[all_rows['SeriesType'] == SERIES_TYPES[i_type]].copy()
    scan_info = init_scan_info()
    scan_files = {}

    for j in range(1, MAX_NSCANS[i_type] + 1):
        if j > len(type_rows):
            continue
        event_info, scan_info = check_series_pcqc(
            event_info, scan_info, type_rows, stype, j, qc_notes_tags)
        event_info, scan_info = check_series_fmap(
            event_info, scan_info, all_rows, type_rows, stype, j, qc_notes_tags)
        event_info = check_series_info(event_info, type_rows, stype, j)
        event_info, scan_info = check_series_autoqc_various(
            event_info, scan_info, type_rows, stype, j)
        event_info, scan_info = check_series_autoqc_subthresh(
            event_info, scan_info, type_rows, stype, j)
        event_info, scan_info = check_series_autoqc_max(
            event_info, scan_info, type_rows, stype, j)
        scan_info = count_series(event_info, scan_info, stype, j)
        scan_files[j] = get_row_val(type_rows.iloc[j - 1], 'filename', '')

    qc_vals = type_rows['QC'].tolist() if 'QC' in type_rows.columns else []
    scan_info['qcd'] = 1 if any(not _is_empty(v) for v in qc_vals) else 0
    return event_info, scan_info, scan_files


def check_event_fmap(event_info, scan_info, stype, qcd=1):
    if stype not in FMAP_STYPES:
        return event_info
    event_info['iqc_{0}_fm_missing'.format(stype)]    = scan_info['fm_missing']
    event_info['iqc_{0}_fm_wronglocus'.format(stype)] = scan_info['fm_wronglocus']
    event_info['iqc_{0}_fm_invalid'.format(stype)]    = scan_info['fm_invalid']
    event_info['iqc_{0}_fm_pc_issues'.format(stype)]  = scan_info['fm_failpc']
    event_info['iqc_{0}_fm_incomp'.format(stype)]     = scan_info['fm_incomp']
    if qcd:
        event_info['iqc_{0}_fm_qcs'.format(stype)] = scan_info['fm_failqc']
        for issue in qu_issues(stype):
            event_info['iqc_{0}_fm_{1}_qc'.format(stype, issue)] = \
                scan_info['fm_' + issue]
    return event_info


def check_event_autoqc(event_info, scan_info, stype, qcd=1):
    for t in autoqc_tlist(stype):
        m_tag = t.lower()
        if not scan_info['total_ser']:
            continue
        # the mean over every scan is cmig_rawqc_bids's, needing only the
        # automated QC; no scan carrying the measure leaves a mean empty
        # rather than 0, the accumulator having stayed at its initial zero
        event_info['iqc_{0}_ok_{1}'.format(stype, m_tag)] = \
            (scan_info['ok_' + m_tag] / scan_info['ok_ser']
             if scan_info['ok_ser'] and scan_info.get('n_ok_' + m_tag) else None)
        event_info['iqc_{0}_comp_{1}'.format(stype, m_tag)] = \
            (scan_info['comp_' + m_tag] / scan_info['comp_ser']
             if scan_info['comp_ser'] and scan_info.get('n_comp_' + m_tag) else None)

    for t in fov_max_tlist(stype):
        if scan_info['total_ser']:
            event_info['iqc_{0}_ok_{1}'.format(stype, t)] = scan_info.get('ok_' + t)

    if stype not in FMRI_STYPES:
        return event_info
    for t in FMRI_SUBTHRESH_TLIST:
        m_tag = t.lower()
        if not scan_info['total_ser']:
            continue
        if stype in FMAP_STYPES:
            event_info['iqc_{0}_ok_{1}'.format(stype, m_tag)] = \
                (scan_info['ok_' + m_tag] if scan_info.get('n_ok_' + m_tag) else None)
            event_info['iqc_{0}_comp_{1}'.format(stype, m_tag)] = \
                (scan_info['comp_' + m_tag] if scan_info.get('n_comp_' + m_tag) else None)
    return event_info


def check_series_event(event_info, scan_info, stype):
    """Write the scan counters into the event row."""
    event_info['iqc_{0}_total_ser'.format(stype)]    = scan_info['total_ser']
    event_info['iqc_{0}_total_passpc'.format(stype)] = scan_info['total_passpc']
    event_info['iqc_{0}_ser_incomp'.format(stype)]   = scan_info['ser_incomp']
    event_info['iqc_{0}_ok_ser'.format(stype)]       = scan_info['ok_ser']
    event_info['iqc_{0}_comp_ser'.format(stype)]     = scan_info['comp_ser']
    event_info['iqc_{0}_ser_pc_issues'.format(stype)] = scan_info['ser_failpc']

    if scan_info['qcd']:
        event_info['iqc_{0}_total_passqc'.format(stype)] = scan_info['total_passqc']
        if stype in FMAP_STYPES:
            event_info['iqc_{0}_good_ser_wofm'.format(stype)] = scan_info['good_ser_wofm']
        event_info['iqc_{0}_good_ser'.format(stype)] = scan_info['good_ser']
        for issue in qu_issues(stype):
            event_info['iqc_{0}_{1}_qc'.format(stype, issue)] = scan_info[issue]
        event_info['iqc_{0}_ser_qcs'.format(stype)] = scan_info['ser_failqc']
    else:
        event_info['iqc_{0}_total_passqc'.format(stype)] = None
        if stype in FMAP_STYPES:
            event_info['iqc_{0}_good_ser_wofm'.format(stype)] = None
        event_info['iqc_{0}_good_ser'.format(stype)] = None
        event_info['iqc_{0}_bad_ser'.format(stype)] = None
        for issue in qu_issues(stype):
            event_info['iqc_{0}_{1}_qc'.format(stype, issue)] = None
        event_info['iqc_{0}_ser_qcs'.format(stype)] = None

    event_info['iqc_{0}_bad_ser'.format(stype)] = scan_info['total_ser'] - scan_info['ok_ser']
    event_info = check_event_fmap(event_info, scan_info, stype, scan_info['qcd'])
    event_info = check_event_autoqc(event_info, scan_info, stype, scan_info['qcd'])
    return event_info


def check_event_pcqc(event_info, series_rows, qc_notes_tags):
    scan_index = {}
    for i, stype in enumerate(STYPES):
        event_info, scan_info, scan_files = check_series_scans(
            event_info, series_rows, stype, i, qc_notes_tags)
        event_info = check_series_event(event_info, scan_info, stype)
        for j, fname in scan_files.items():
            scan_index[(stype, j)] = fname
    return event_info, scan_index


FMAP_PTYPES = [('dmri', 'dmri'), ('fmri', 'fmri')]


def count_fmap_pairs(event_info, all_rows):
    """Count the field map pairs of a session, and how many are usable.

    Judged on the field maps' own rows rather than through the scans that
    select them, so a pair shared by several runs is counted once and a pair
    no scan selects is counted at all. A pair is one forward and one reverse
    sharing an acquisition slot ordering; both members have to be usable for
    the pair to be.
    """
    sertypes = [str(t).lower() for t in all_rows['SeriesType'].tolist()]
    slots = acquisition_slots(all_rows)
    for stem, ptype in FMAP_PTYPES:
        fwd = [i for i, t in enumerate(sertypes) if t == ptype + '_fm_pa']
        rev = [i for i, t in enumerate(sertypes) if t == ptype + '_fm_ap']
        n_ok = n_comp = 0
        # pair them in acquisition order, which is how they were collected
        for f_pos, r_pos in zip(sorted(fwd, key=lambda p: slots[p]),
                                sorted(rev, key=lambda p: slots[p])):
            rows = [all_rows.iloc[f_pos], all_rows.iloc[r_pos]]
            done = [get_row_val(r, 'Completed', 0) for r in rows]
            comp = [str(get_row_val(r, 'Compliant', '')).strip().lower()
                    for r in rows]
            qc = [get_row_val(r, 'QC', None) for r in rows]
            if all(done) and all(c != 'bad' for c in comp) and \
               all(_is_empty(q) or q == 1 for q in qc):
                n_ok += 1
            if all(done) and all(c == 'yes' for c in comp):
                n_comp += 1
        event_info['iqc_{0}_fm_ok'.format(stem)] = n_ok
        event_info['iqc_{0}_fm_passpc_comp'.format(stem)] = n_comp
    return event_info


def check_event_completeness(event_info):
    event_info['iqc_eventcomplete'] = 1
    for i, stype in enumerate(STYPES):
        if (event_info.get('iqc_{0}_total_ser'.format(stype)) or 0) < MIN_NSCANS[i]:
            event_info['iqc_eventcomplete'] = 0
            break
    # over ok_ser rather than good_ser: good_ser needs a QC score of 1, so a
    # session nobody reviewed falls short on every series type at once, which
    # reads as unreviewed rather than as unusable
    event_info['iqc_eventcomplete_ok'] = 1
    for i, stype in enumerate(STYPES):
        if (event_info.get('iqc_{0}_ok_ser'.format(stype)) or 0) < MIN_NSCANS[i]:
            event_info['iqc_eventcomplete_ok'] = 0
            break
    return event_info


def series_datetime(date, time):
    """The acquisition date-time of a series, for ordering.

    SeriesTime is HHMMSS with an optional fraction, possibly with colons or
    without a leading zero, so sorting it as a string can misorder series.
    Where either field does not parse, the series sorts after those that do,
    in string order.
    """
    try:
        hms, _, frac = str(time).replace(':', '').strip().partition('.')
        hms = hms.zfill(6)
        stamp = datetime.datetime.strptime(str(date).strip() + hms,
                                           '%Y%m%d%H%M%S')
        if frac:
            stamp += datetime.timedelta(microseconds=int(frac[:6].ljust(6, '0')))
        return (0, stamp, '', '')
    except ValueError:
        return (1, datetime.datetime.min, str(date), str(time))


def check_event(pid, sid, series_rows):
    """Roll one session's series rows up into its event row."""
    valid = series_rows[series_rows['StudyDate'].apply(lambda d: not _is_empty(d))]
    if len(valid) < len(series_rows):
        print('{0}: WARNING: excluding series from {1} {2} because of missing '
              'info'.format(APPNAME, pid, sid))
    series_rows = valid.copy()
    order = sorted(range(len(series_rows)), key=lambda i: series_datetime(
        series_rows['StudyDate'].iloc[i], series_rows['SeriesTime'].iloc[i]))
    series_rows = series_rows.iloc[order].reset_index(drop=True)

    qc_notes_tags = find_qc_notes(series_rows)
    event_info = init_event_info(pid, sid, series_rows)
    event_info, scan_index = check_event_pcqc(event_info, series_rows, qc_notes_tags)
    event_info = count_fmap_pairs(event_info, series_rows)
    event_info = check_event_completeness(event_info)
    return event_info, scan_index, series_rows

#############################################################################
# output
#############################################################################

PER_RUN = re.compile(r'^iqc_[a-z0-9]+_\d+_')

# the manual QC of the scan's field map pair, the one per-scan value that needs
# pc's match and the manual QC together; pc's series TSV carries the rest. The
# series data is one table, so the columns take its name as their prefix
SERIES_TABLE = 'mr_y_qc__raw__series'
SERIES_TAGS_OUT = ['fm_qc_score', 'fm_qc_notes']
# the series column of each of those
SERIES_COLUMN_OUT = {'fm_qc_score': SERIES_TABLE + '__fm__qc_score',
                     'fm_qc_notes': SERIES_TABLE + '__fm__qc_notes'}

BIDS_ENTITIES = ['task', 'acq', 'dir', 'run', 'rec']

# the released name of every across-scan column is its stem and one suffix,
# the suffix being the same on every stem that defines the column
COUNT_NAMES = {
    'total_ser':       'series_count',
    'total_passpc':    'pass__protcom_count',
    'total_passqc':    'pass__qc_count',
    'good_ser':        'pass__qc__protcom_count',
    'good_ser_fmri':   'pass__qc__protcom__comp_count',
    'good_ser_wofm':   'pass__qc__protcom__nofm_count',
    'ok_ser':          'ok_count',
    'comp_ser':        'pass__protcom__comp_count',
    'bad_ser':         'fail_count',
    'ser_incomp':      'incomp_count',
    'ser_pc_issues':   'fail__protcom_count',
    'ser_qcs':         'fail__qc_count',
    'fm_missing':      'fm__miss_count',
    'fm_wronglocus':   'fm__wrngloc_count',
    'fm_invalid':      'fm__invalid_count',
    'fm_pc_issues':    'fm__fail__protcom_count',
    'fm_incomp':       'fm__incomp_count',
    'fm_qcs':          'fm__fail__qc_count',
}
ISSUE_NAMES = {
    'mo':   'artif__mot_count',
    'fov':  'fov__cutoff_count',
    'sus':  'artif__susc_count',
    'line': 'line_count',
    'fa':   'artif__famap_count',
    'gh':   'artif__ghost_count',
    'moire': 'artif__moire_count',
}

# the aggregate over scans makes a mean of a series-level value that is not
# itself a statistic; the frame count is named frame, as in the series file
AGG_MEASURE = {'rep_count': 'frame_mean', 'tr': 'tr_mean'}

# the two columns that span every series type rather than sitting under one
FMAP_PAIR_NAMES = {
    'fm_ok':          'fm__ok_count',
    'fm_passpc_comp': 'fm__pass__protcom__comp_count',
}

EVENT_NAMES = {
    'eventcomplete':    'event__full__protocol_indicator',
    'eventcomplete_ok': 'event__full__protocol__ok_indicator',
}


def scoped_name(stype, scope, measure):
    """The released name of one across-scan aggregate.

    The scope goes in front of the measure word, which is what separates
    rsfmri__mot__ok_mean from the plain rsfmri__mot_mean the rawqc module
    writes.
    """
    m = AGG_MEASURE.get(measure, measure)
    i = m.rfind('_')
    return 'mr_y_qc__raw__{0}__{1}__{2}{3}'.format(
        STEM_BY_STYPE[stype], m[:i], scope, m[i:])


def released_name(col):
    """The released name of an event column, or None where it has none."""
    if col[4:] in EVENT_NAMES:
        return 'mr_y_qc__raw__{0}'.format(EVENT_NAMES[col[4:]])
    m = re.match(r'^iqc_([a-z0-9]+)_(.+)$', col)
    if not m:
        return None
    stype, tag = m.group(1), m.group(2)
    # the field map pair counts sit under the modality, not a series type
    if tag in FMAP_PAIR_NAMES and stype in [p[0] for p in FMAP_PTYPES]:
        return 'mr_y_qc__raw__{0}__{1}'.format(stype, FMAP_PAIR_NAMES[tag])
    if stype not in STEM_BY_STYPE:
        return None
    stem = STEM_BY_STYPE[stype]
    # the fMRI stems spell good_ser with the complete their description names,
    # where t1, t2 and dmri leave it out of both
    if tag == 'good_ser' and stype in FMRI_STYPES:
        tag = 'good_ser_fmri'
    if tag in COUNT_NAMES:
        return 'mr_y_qc__raw__{0}__{1}'.format(stem, COUNT_NAMES[tag])
    m = re.match(r'^(fm_)?([a-z]+)_qc$', tag)
    if m and m.group(2) in ISSUE_NAMES:
        name = ISSUE_NAMES[m.group(2)]
        if m.group(1):
            name = 'fm__' + name
        return 'mr_y_qc__raw__{0}__{1}'.format(stem, name)
    m = re.match(r'^(ok|comp)_(.+)$', tag)
    if m:
        return scoped_name(stype, m.group(1), m.group(2))
    return None


def bids_entities(fname):
    """The entities of a BIDS filename, keyed as the rawqc series TSVs key them."""
    stem = re.sub(r'\.nii(\.gz)?$', '', os.path.basename(fname))
    out = dict((k, '') for k in BIDS_ENTITIES)
    for part in stem.split('_')[:-1]:
        if '-' in part:
            key, val = part.split('-', 1)
            if key in out:
                out[key] = val
    return out


# the released variables that need protocol compliance alone, which
#   cmig_pc_bids_abcd writes from a check_event of its own; this module still works
#   them out, its QC counts being built from them, but leaves them out of the
#   session file so that each is written once
PC_ONLY_TAGS = {
    'total_ser', 'total_passpc', 'ser_pc_issues', 'ser_incomp', 'comp_ser',
    'fm_missing', 'fm_wronglocus', 'fm_invalid', 'fm_pc_issues', 'fm_incomp',
    'comp_rep_count', 'comp_tr', 'fm_passpc_comp',
}
PC_ONLY_EVENT = {'eventcomplete'}


def is_pc_only(col):
    """True for an event column cmig_pc_bids_abcd writes in place of this one."""
    if col[4:] in PC_ONLY_EVENT:
        return True
    m = re.match(r'^iqc_([a-z0-9]+)_(.+)$', col)
    return bool(m) and m.group(2) in PC_ONLY_TAGS


def na_strings(row):
    """A row with its empty strings as n/a (absent values are written as n/a
    by to_csv's na_rep), as BIDS requires of a TSV."""
    return dict((k, 'n/a' if isinstance(v, str) and v.strip() == '' else v)
                for k, v in row.items())


def write_pcqc_tsv(fpath, event_info):
    """Write the event row, one row per session, without the per-run columns.

    Those values go to the series file, where the row identifies the scan and
    no run index is needed.
    """
    row = {}
    for key, val in event_info.items():
        if PER_RUN.match(key) or is_pc_only(key):
            continue
        # a column with no released name keeps the working one, which is the
        # duplicate and multiples family moving to cmig_pc_bids_abcd and the
        # two completeness flags, neither of which the dictionary carries
        row[released_name(key) or key] = val
    df = pd.DataFrame([na_strings(row)], columns=list(row.keys()))
    df.to_csv(fpath, sep='\t', index=False, na_rep='n/a')


def write_series_pcqc_tsv(fpath, pid, sid, event_info, scan_index, series_rows):
    """Write one row per scan, carrying what the merge works out for it."""
    rows = []
    for (stype, j), fname in sorted(scan_index.items(), key=lambda kv: kv[1]):
        # the series type and the uid come from cmig_pc_bids_abcd, which the
        # compiler joins to this on filename
        row = {'participant_id': pid, 'session_id': sid, 'filename': fname}
        row.update(bids_entities(fname))
        for tag in SERIES_TAGS_OUT:
            row[SERIES_COLUMN_OUT[tag]] = \
                event_info.get('iqc_{0}_{1}_{2}'.format(stype, j, tag))
        rows.append(row)
    cols = (['participant_id', 'session_id', 'filename'] + BIDS_ENTITIES +
            [SERIES_COLUMN_OUT[t] for t in SERIES_TAGS_OUT])
    rows = [na_strings(r) for r in rows]
    pd.DataFrame(rows, columns=cols).to_csv(fpath, sep='\t', index=False, na_rep='n/a')


def write_status(outdir, pid, sid, errcode, errmsg, startdate=None):
    fpath = cmig_status.cmig_write_status(
        pid, sid, outdir, 'desc-pcqc', errcode, errmsg, startdate=startdate)
    if errcode == 0:
        print('{0}: finished'.format(APPNAME))
        print('{0}: completed successfully'.format(APPNAME))
    cmig_sing_utils.add_status_versions(fpath, app_name=APPNAME)


def input_status_paths(pc_rootdir, rawqc_rootdir, manualqc_rootdir, pid, sid):
    """The status files of the three modules this one reads."""
    paths = [
        cmig_status.status_fpath(
            os.path.join(pc_rootdir, pid, sid), pid, sid, 'desc-pc'),
        cmig_status.status_fpath(
            os.path.join(manualqc_rootdir, pid, sid), pid, sid,
            'desc-manualrawqc'),
    ]
    paths.extend(sorted(glob.glob(os.path.join(
        rawqc_rootdir, pid, sid, '*',
        '{0}_{1}_desc-rawqc_status.tsv'.format(pid, sid)))))
    return paths


def check_prior_status(outdir, pid, sid, forceflag, input_paths=None):
    """True where a prior successful run exists, forceflag is 0, and no input
    has run since.

    The output is a merge of three modules' output, so a prior run is only
    good while none of them has run again.
    """
    fpath = cmig_status.status_fpath(outdir, pid, sid, 'desc-pcqc')
    status = cmig_status.read_status(fpath)
    if status is None or forceflag:
        return False
    try:
        if int(status.get('errcode', 1)) != 0:
            return False
    except ValueError:
        return False
    return not cmig_status.inputs_newer(fpath, input_paths or [])

#############################################################################
# CLI
#############################################################################

def build_parser():
    parser = argparse.ArgumentParser(
        description='Merge protocol compliance with QC for one BIDS session.')
    parser.add_argument('--participant_id', required=True,
                        help='Participant ID (e.g., sub-XXXXXXXX)')
    parser.add_argument('--session_id', required=True,
                        help='Session ID (e.g., ses-00A)')
    parser.add_argument('--pc_rootdir', required=True,
                        help='Root directory of protocol compliance output')
    parser.add_argument('--rawqc_rootdir', required=True,
                        help='Root directory of automated raw QC output')
    parser.add_argument('--manualqc_rootdir', required=True,
                        help='Root directory of manual raw QC output')
    parser.add_argument('--output_rootdir', required=True,
                        help='Root directory for pcqc output (created if absent)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='Display status messages (default: 0)')
    parser.add_argument('--forceflag', type=int, default=1, choices=[0, 1],
                        help='Overwrite existing output (default: 1)')
    return parser


def parse_args():
    return build_parser().parse_args()


def run(args):
    pid = args.participant_id
    sid = args.session_id
    outdir = os.path.join(args.output_rootdir, pid, sid)
    startdate = cmig_status.module_start()

    input_paths = input_status_paths(args.pc_rootdir, args.rawqc_rootdir,
                                     args.manualqc_rootdir, pid, sid)
    if check_prior_status(outdir, pid, sid, args.forceflag, input_paths):
        if args.verbose:
            print('{0}: prior successful run found, skipping'.format(APPNAME))
        return 0, ''

    pc_path = os.path.join(args.pc_rootdir, pid, sid,
                           '{0}_{1}_desc-series_pc.tsv'.format(pid, sid))
    if not os.path.isfile(pc_path):
        errmsg = 'series pc file not found: {0}'.format(pc_path)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg

    pc_df = read_tsv(pc_path)
    if pc_df.empty:
        errmsg = 'no series in {0}'.format(pc_path)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg
    # a series file from before the field map pair columns cannot be read
    missing = [c for c in ('mr_y_qc__raw__series_type',
                           SERIES_FM_COLUMN['fm_missing'])
               if c not in pc_df.columns]
    if missing:
        errmsg = 'columns {0} not found in {1}'.format(', '.join(missing), pc_path)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg

    consensus, notes, reviewers = read_manual_qc(args.manualqc_rootdir, pid, sid)
    if not consensus and args.verbose:
        print('{0}: no manual QC for {1} {2}'.format(APPNAME, pid, sid))
    measures = read_rawqc(args.rawqc_rootdir, pid, sid)
    if not measures and args.verbose:
        print('{0}: no automated QC for {1} {2}'.format(APPNAME, pid, sid))

    series_rows = build_series_table(pc_df, consensus, notes, reviewers, measures)
    series_rows = get_relevant_series(series_rows)
    if series_rows.empty:
        errmsg = 'no ABCD-relevant series for {0} {1}'.format(pid, sid)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg
    series_rows = convert_seriestype(series_rows, ['T1_NORM', 'T2_NORM'], ['T1', 'T2'])

    event_info, scan_index, sorted_rows = check_event(pid, sid, series_rows)

    os.makedirs(outdir, exist_ok=True)
    tsv_path = os.path.join(outdir, '{0}_{1}_desc-session_pcqc.tsv'.format(pid, sid))
    write_pcqc_tsv(tsv_path, event_info)
    ser_path = os.path.join(outdir, '{0}_{1}_desc-series_pcqc.tsv'.format(pid, sid))
    write_series_pcqc_tsv(ser_path, pid, sid, event_info, scan_index, sorted_rows)

    if not os.path.isfile(tsv_path):
        errmsg = 'failed to write {0}'.format(tsv_path)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg

    if args.verbose:
        print('{0}: wrote {1}'.format(APPNAME, tsv_path))

    write_status(outdir, pid, sid, 0, '', startdate)
    return 0, ''


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
