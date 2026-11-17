#!/usr/bin/env python3
"""
ABCD protocol compliance checking using BIDS rawdata input.

Created:  05/13/2026 by Don Hagler
Prev mod: 05/28/2026 by Don Hagler
Prev mod: 06/17/2026 by Don Hagler
Prev mod: 08/10/2026 by Don Hagler
Prev mod: 09/24/2026 by Don Hagler
Prev mod: 09/25/2026 by Don Hagler
Prev mod: 10/02/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/08/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Last mod: 11/15/2026 by Don Hagler
"""

import argparse
import collections
import csv
import datetime
import glob
import json
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='dicom')
import nibabel as nib
import os
import re
import sys

import cmig_sing_utils
import cmig_write_status as cmig_status
from cmig_select_bids_func_runs import (read_events, is_old_style,
                                        bold_acq_datetime, read_scans_usable,
                                        run_usable, read_scans_column,
                                        run_compliant)

APPNAME  = 'cmig_pc_bids_abcd'
MODALITIES = ['anat', 'func', 'dwi', 'fmap']

# Series types that are inherently 3D (single volume); Completed='Yes' if file exists
THREE_D_SERIES = {
    'T1', 'T1_NORM', 'T2', 'T2_NORM',
    'fMRI_FM_AP', 'fMRI_FM_PA', 'dMRI_FM_AP', 'dMRI_FM_PA',
}

# Expected volume counts per (manufacturer_key, series_type) or
# (manufacturer_key, series_type, sw_tag) for GE
EXPECTED_NVOLS = {
    # Siemens
    ('Siemens', 'rsfMRI'):          383,
    ('Siemens', 'fMRI_MID_task'):   411,
    ('Siemens', 'fMRI_SST_task'):   445,
    ('Siemens', 'fMRI_nBack_task'): 370,
    ('Siemens', 'dMRI'):            103,
    # Philips
    ('Philips', 'rsfMRI'):          383,
    ('Philips', 'fMRI_MID_task'):   411,
    ('Philips', 'fMRI_SST_task'):   445,
    ('Philips', 'fMRI_nBack_task'): 370,
    ('Philips', 'dMRI'):            51,
    # GE DV25
    ('GE', 'rsfMRI',          'DV25'): 380,
    ('GE', 'fMRI_MID_task',   'DV25'): 408,
    ('GE', 'fMRI_SST_task',   'DV25'): 442,
    ('GE', 'fMRI_nBack_task', 'DV25'): 367,
    ('GE', 'dMRI',            'DV25'): 103,
    # GE other (incl. DV26)
    ('GE', 'rsfMRI',          ''):     391,
    ('GE', 'fMRI_MID_task',   ''):     419,
    ('GE', 'fMRI_SST_task',   ''):     453,
    ('GE', 'fMRI_nBack_task', ''):     378,
    ('GE', 'dMRI',            ''):     104,
}

# Output TSV columns in order
INFO_COLUMNS = [
    'participant_id', 'session_id',
    'mr_y_adm__info_dt',
    'mr_y_adm__info_sessid',
    'mr_y_adm__info__dev_manufact',
    'mr_y_adm__info__dev_model',
    'mr_y_adm__info__dev_serial',
    'mr_y_adm__info__dev_fldstr',
    'mr_y_adm__info__dev__sftw_ver',
]

# pc's session columns take the released mr_y_qc__raw__ prefix, for
#   consistency with earlier releases, each under its series type's table
RAW = 'mr_y_qc__raw__'

# session variables that need only the compliance: the session's duration,
# the per-type duplicate counts, and whether the session was assembled from
# more than one study
SESSION_EVENT_COLUMNS = [
    RAW + 'event_dur',
    RAW + 'event__study__multi_indicator',
    RAW + 'event__study__dt__multi_indicator',
    RAW + 'event__dup_max',
] + [RAW + '{0}__dup_count'.format(st) for st in
     ['smri__t1', 'smri__t2', 'dmri', 'rsfmri',
      'tfmri__mid', 'tfmri__sst', 'tfmri__nback']]

# minutes added to the last series start for the session's duration, the
#   series' own length not being recorded
EVENT_DUR_LAST_MIN = 5

# the stems with a field map, whose released pass count requires it to pass
NOFM_STEMS = ['rsfmri', 'tfmri__mid', 'tfmri__sst', 'tfmri__nback', 'dmri']

SESSION_PC_COLUMNS = (
    ['participant_id', 'session_id'] +
    [RAW + 'smri__t1__norm__series_count',
     RAW + 'smri__t2__norm__series_count'] +
    [RAW + '{0}__{1}'.format(st, tail) for st in NOFM_STEMS
     for tail in ('pass__protcom__nofm_count', 'fail__protcom__comp_count')] +
    [RAW + '{0}__fm__{1}__{2}__protcom_count'.format(mod, d, pf)
     for mod in ('fmri', 'dmri') for d in ('ap', 'pa')
     for pf in ('pass', 'fail')] +
    [RAW + '{0}__fm__{1}'.format(mod, tail) for mod in ('fmri', 'dmri')
     for tail in ('pair_count', 'pair__invalid_count')])

# the stem each series type is counted under, NORM with its plain type as the
#   released counts fold it
SERIES_TYPE_TO_STEM = {
    'T1':              'smri__t1',
    'T1_NORM':         'smri__t1',
    'T2':              'smri__t2',
    'T2_NORM':         'smri__t2',
    'rsfMRI':          'rsfmri',
    'fMRI_MID_task':   'tfmri__mid',
    'fMRI_SST_task':   'tfmri__sst',
    'fMRI_nBack_task': 'tfmri__nback',
    'fMRI_FM_AP':      'fmri__fm__ap',
    'fMRI_FM_PA':      'fmri__fm__pa',
    'dMRI':            'dmri',
    'dMRI_FM_AP':      'dmri__fm__ap',
    'dMRI_FM_PA':      'dmri__fm__pa',
}

SESSION_PC_COLUMNS = SESSION_PC_COLUMNS + SESSION_EVENT_COLUMNS

# the released variables that need protocol compliance alone, recreated here
#   with the release's definitions by check_event below. The series types and
#   their minimums are those of cmig_check_event_pcqc, from
#   abcd_check_event_pcqc; unlike it, every series of a type is counted, the
#   release's maximum having served only its per-run columns
CHECK_STYPES      = ['t1', 't2', 'dmri', 'rsfmri', 'mid', 'sst', 'nback']
CHECK_SERIES_TYPES = ['T1', 'T2', 'dMRI', 'rsfMRI',
                       'fMRI_MID_task', 'fMRI_SST_task', 'fMRI_nBack_task']
CHECK_MIN_NSCANS  = [1, 1, 1, 3, 2, 2, 2]
CHECK_FMAP_STYPES = ['dmri', 'rsfmri', 'mid', 'sst', 'nback']
CHECK_FMRI_STYPES = ['rsfmri', 'mid', 'sst', 'nback']
CHECK_STEM = {
    't1': 'smri__t1', 't2': 'smri__t2', 'dmri': 'dmri', 'rsfmri': 'rsfmri',
    'mid': 'tfmri__mid', 'sst': 'tfmri__sst', 'nback': 'tfmri__nback',
}
# the series measures averaged over complete and compliant scans
CHECK_COMP_MEASURES = {'dmri': ['frame'], 'rsfmri': ['frame', 'tr'],
                        'mid': ['frame', 'tr'], 'sst': ['frame', 'tr'],
                        'nback': ['frame', 'tr']}
CHECK_COUNT_NAMES = [
    ('total_ser',  'series_count'),
    ('total_passpc', 'pass__protcom_count'),
    ('ser_failpc', 'fail__protcom_count'),
    ('ser_incomp', 'incomp_count'),
    ('comp_ser',   'pass__protcom__comp_count'),
]
CHECK_FMAP_NAMES = [
    ('fm_missing',    'fm__miss_count'),
    ('fm_wronglocus', 'fm__wrngloc_count'),
    ('fm_invalid',    'fm__invalid_count'),
    ('fm_failpc',     'fm__fail__protcom_count'),
    ('fm_incomp',     'fm__incomp_count'),
]
CHECK_PAIR_MODALITIES = ['dmri', 'fmri']


def released_pc_columns():
    """The released names check_event writes, in output order."""
    cols = []
    for stype in CHECK_STYPES:
        stem = CHECK_STEM[stype]
        for _, name in CHECK_COUNT_NAMES:
            cols.append('mr_y_qc__raw__{0}__{1}'.format(stem, name))
        if stype in CHECK_FMAP_STYPES:
            for _, name in CHECK_FMAP_NAMES:
                cols.append('mr_y_qc__raw__{0}__{1}'.format(stem, name))
        for m in CHECK_COMP_MEASURES.get(stype, []):
            cols.append('mr_y_qc__raw__{0}__{1}__comp_mean'.format(stem, m))
    for mod in CHECK_PAIR_MODALITIES:
        for d in ('ap', 'pa'):
            cols.append('mr_y_qc__raw__{0}__fm__{1}__series_count'.format(mod, d))
        cols.append('mr_y_qc__raw__{0}__fm__pass__protcom__comp_count'.format(mod))
    cols.append('mr_y_qc__raw__event__full__protocol_indicator')
    return cols


SESSION_PC_COLUMNS = SESSION_PC_COLUMNS + released_pc_columns()

# E-Prime variables from the rawdata events files, per task: four run counts,
#   new, and the release's names, recreated from the events and bold jsons
EPRIME_TASKS = ['mid', 'sst', 'nback']
EPRIME_NAMES = [
    'scan__run_count', 'scan__run__valid_count', 'behav__run_count',
    'noscan__run_count', 'file_indicator', 'read_indicator',
    'content_indicator', 'corrtask_indicator', 'match_indicator',
    'start_t', 'dt', 'dt_indicator', 'tdiff_min', 'tdiff_max',
]
# the release's spelling: eprime_dt, not eprime__dt
EPRIME_COLUMN = dict(
    (n, ('{0}eprime_dt' if n == 'dt' else '{0}eprime__' + n)) for n in EPRIME_NAMES)
EPRIME_MAX_DELAY = 12.5   # minutes, the eprime-scan limit match_indicator tests


def eprime_columns():
    return [EPRIME_COLUMN[n].format(RAW + 'tfmri__{0}__'.format(task))
            for task in EPRIME_TASKS for n in EPRIME_NAMES]


SESSION_PC_COLUMNS = SESSION_PC_COLUMNS + eprime_columns()

# the series column of each field map pair value check_fmap returns
SERIES_FM_COLUMN = collections.OrderedDict(
    (tag, 'mr_y_qc__raw__series__' + col) for tag, col in (
        ('fm_missing', 'fm__missing_indicator'),
        ('fm_wronglocus', 'fm__wronglocus_indicator'),
        ('fm_invalid', 'fm__invalid_indicator'),
        ('fm_pa_filename', 'fm__pa_filename'),
        ('fm_ap_filename', 'fm__ap_filename'),
        ('fm_pc_score', 'fm__pc__pass_indicator'),
        ('fm_pc_bad', 'fm__pc__bad_indicator'),
        ('fm_pc_notes', 'fm__pc_notes'),
        ('fm_complete', 'fm__complete_indicator')))

TSV_COLUMNS = [
    'participant_id', 'session_id', 'filename',
    'mr_y_qc__raw__series_num',
    'mr_y_qc__raw__series_desc',
    'mr_y_qc__raw__series_type',
    'mr_y_qc__raw__series__protcom_indicator',
    'mr_y_qc__raw__series__psn__match_indicator',
    'mr_y_qc__raw__series__desc__match_indicator',
    'mr_y_qc__raw__series__frame_count',
    'mr_y_qc__raw__series__complete_indicator',
    'mr_y_qc__raw__series_notes',
    'mr_y_qc__raw__series__dev_manufact',
    'mr_y_qc__raw__series__dev_model',
    'mr_y_qc__raw__series__dev__sftw_ver',
    'mr_y_qc__raw__series__dev_serial',
    'mr_y_qc__raw__series_psn',
    'mr_y_qc__raw__series__img_type',
    'mr_y_qc__raw__series_position',
    'mr_y_qc__raw__series__acq__matrix_pe',
    'mr_y_qc__raw__series__pe_steps',
    'mr_y_qc__raw__series__pe__fov_pct',
    'mr_y_qc__raw__series_tr',
    'mr_y_qc__raw__series_te',
    'mr_y_qc__raw__series__pixel_bw',
    'mr_y_qc__raw__series__study_dt',
    'mr_y_qc__raw__series__study_uid',
    'mr_y_qc__raw__series__series_uid',
    'mr_y_qc__raw__series__series_t',
    'mr_y_qc__raw__series__acq_t',
    'mr_y_qc__raw__series__coil_str',
    'mr_y_qc__raw__series__coil_elements',
    'mr_y_qc__raw__series__coil_name',
] + list(SERIES_FM_COLUMN.values())

#############################################################################
# BIDS utilities
#############################################################################

def parse_bids_filename(fname):
    # BIDS filenames are underscore-delimited key-value pairs followed by a suffix,
    # e.g. sub-X_ses-Y_task-rest_run-01_bold.nii.gz -> suffix='bold', task='rest', run='01'
    stem = re.sub(r'\.nii(\.gz)?$', '', os.path.basename(fname))
    parts = stem.split('_')
    info = {'suffix': parts[-1]}  # suffix is always the last token
    for part in parts[:-1]:
        if '-' in part:
            key, val = part.split('-', 1)
            info[key] = val
    return info


def get_bids_scans(participant_id, session_id, raw_rootdir, verbose=False):
    # session directory: {raw_rootdir}/{participant_id}/{session_id}/
    sess_dir = os.path.join(raw_rootdir, participant_id, session_id)
    if not os.path.isdir(sess_dir):
        print('{0}: WARNING: session directory not found: {1}'.format(
            APPNAME, sess_dir))
        return []

    scans = []
    for modality in MODALITIES:
        mod_dir = os.path.join(sess_dir, modality)
        if not os.path.isdir(mod_dir):
            continue

        # collect .nii and .nii.gz files; sort for deterministic ordering
        nii_files = sorted(
            glob.glob(os.path.join(mod_dir, '*.nii')) +
            glob.glob(os.path.join(mod_dir, '*.nii.gz'))
        )

        for nii_path in nii_files:
            # JSON sidecar has the same path as the NIfTI with .json extension
            json_path = re.sub(r'\.nii(\.gz)?$', '.json', nii_path)
            if not os.path.isfile(json_path):
                print('{0}: WARNING: no JSON sidecar found for {1}'.format(
                    APPNAME, nii_path))
                json_path = None

            # parse BIDS filename entities for use in series identification
            info = parse_bids_filename(nii_path)
            scans.append({
                'nii_path':  nii_path,
                'json_path': json_path,
                'modality':  modality,
                'suffix':    info.get('suffix'),   # T1w, T2w, bold, dwi, epi
                'rec':       info.get('rec'),       # e.g. norm
                'task':      info.get('task'),      # e.g. rest, MID, nBack, SST
                'acq':       info.get('acq'),       # e.g. func, dwi
                'dir':       info.get('dir'),       # e.g. AP, PA
                'run':       info.get('run'),       # e.g. 01
                'series_type': None,               # filled in by run()
            })
            if verbose:
                print('{0}: found {1}'.format(APPNAME, nii_path))

    return scans

#############################################################################
# Metadata extraction
#############################################################################

def load_metadata(json_path):
    """Load a BIDS JSON sidecar. Returns empty dict if path is None or unreadable."""
    if not json_path:
        return {}
    try:
        with open(json_path, 'r') as fid:
            return json.load(fid)
    except (IOError, ValueError) as err:
        print('{0}: WARNING: could not read {1}: {2}'.format(
            APPNAME, json_path, err))
        return {}


def extract_info_fields(meta):
    """Extract session-level scanner metadata for the info TSV."""
    sw = meta.get('SoftwareVersions', '')
    if isinstance(sw, list):
        sw = '|'.join(str(v) for v in sw)
    return {
        'mr_y_adm__info_dt':            meta.get('StudyDate', ''),
        'mr_y_adm__info__dev_manufact': meta.get('Manufacturer', ''),
        'mr_y_adm__info__dev_model':    meta.get('ManufacturersModelName', ''),
        'mr_y_adm__info__dev_serial':   meta.get('DeviceSerialNumber', ''),
        'mr_y_adm__info__dev_fldstr':   meta.get('MagneticFieldStrength', ''),
        'mr_y_adm__info__dev__sftw_ver': sw,
    }


def extract_fields(meta):
    """Extract the BIDS JSON fields needed for compliance checks."""
    image_type = meta.get('ImageType', '')
    if isinstance(image_type, list):
        image_type = '\\'.join(image_type)
    return {
        'Manufacturer':              meta.get('Manufacturer', ''),
        'ManufacturersModelName':    meta.get('ManufacturersModelName', ''),
        'SoftwareVersions':          meta.get('SoftwareVersions', ''),
        'DeviceSerialNumber':        meta.get('DeviceSerialNumber', ''),
        'PulseSequenceName':         meta.get('PulseSequenceName', ''),
        'SeriesDescription':         meta.get('SeriesDescription', ''),
        'ScanOptions':               meta.get('ScanOptions', ''),
        'ImageType':                 image_type,
        'PatientPosition':           meta.get('PatientPosition', ''),
        'AcquisitionMatrixPE':       meta.get('AcquisitionMatrixPE', ''),
        'PhaseEncodingSteps':        meta.get('PhaseEncodingSteps', ''),
        'PercentPhaseFOV':           meta.get('PercentPhaseFOV', ''),
        'RepetitionTime':            meta.get('RepetitionTime', ''),
        'EchoTime':                  meta.get('EchoTime', ''),
        'PixelBandwidth':            meta.get('PixelBandwidth', ''),
        'PatientID':                 meta.get('PatientID', ''),
        'PatientName':               meta.get('PatientName', ''),
        'StudyDate':                 meta.get('StudyDate', ''),
        'StudyInstanceUID':          meta.get('StudyInstanceUID', ''),
        'SeriesInstanceUID':         meta.get('SeriesInstanceUID', ''),
        'SeriesNumber':              meta.get('SeriesNumber', ''),
        'SeriesTime':                meta.get('SeriesTime', ''),
        'AcquisitionTime':           meta.get('AcquisitionTime', ''),
        'NumberOfTemporalPositions':  meta.get('NumberOfTemporalPositions', ''),
        'ScanningSequence':           meta.get('ScanningSequence', ''),
        'SequenceVariant':            meta.get('SequenceVariant', ''),
        'Rows':                       meta.get('Rows', ''),
        'CoilString':                 meta.get('CoilString', ''),
        'ReceiveCoilActiveElements':  meta.get('ReceiveCoilActiveElements', ''),
        'ReceiveCoilName':            meta.get('ReceiveCoilName', ''),
    }

#############################################################################
# Series type identification from BIDS filename entities
#############################################################################

def get_series_type(scan):
    """Determine series type from BIDS filename entities."""
    modality = scan['modality']
    suffix   = scan.get('suffix') or ''
    rec      = (scan.get('rec') or '').lower()
    acq      = (scan.get('acq') or '').lower()

    # NORM is identified from the rec entity alone, for every manufacturer;
    #   rawdata writes it as rec-norm
    if modality == 'anat':
        if suffix == 'T1w':
            return 'T1_NORM' if rec == 'norm' else 'T1'
        if suffix == 'T2w':
            return 'T2_NORM' if rec == 'norm' else 'T2'
    elif modality == 'func':
        if suffix == 'bold':
            return 'fMRI'
    elif modality == 'dwi':
        if suffix == 'dwi':
            return 'dMRI'
        if suffix == 'epi':
            return 'dMRI_FM'
    elif modality == 'fmap':
        if suffix == 'epi':
            if 'dwi' in acq or 'dti' in acq or 'dmri' in acq:
                return 'dMRI_FM'
            return 'fMRI_FM'

    return 'Undefined'

#############################################################################
# Compliance check helpers
#############################################################################

def _fmri_task_type(scan):
    """Map BIDS task entity to fMRI series type label."""
    task = (scan.get('task') or '').lower()
    mapping = {
        'rest':  'rsfMRI',
        'mid':   'fMRI_MID_task',
        'sst':   'fMRI_SST_task',
        'nback': 'fMRI_nBack_task',
    }
    return mapping.get(task, 'fMRI')


def get_nifti_nvols(nii_path):
    """Return number of volumes in a NIfTI file (1 for 3D)."""
    try:
        shape = nib.load(nii_path).shape
        return shape[3] if len(shape) >= 4 else 1
    except Exception:
        return None


def get_expected_nvols(manufacturer, series_type, sw_versions):
    """Return expected volume count for this series, or None if not defined."""
    mfr = manufacturer or ''
    sw  = sw_versions or ''
    if 'GE' in mfr:
        sw_tag = 'DV25' if 'DV25' in sw else ''
        return EXPECTED_NVOLS.get(('GE', series_type, sw_tag))
    elif 'Siemens' in mfr:
        return EXPECTED_NVOLS.get(('Siemens', series_type))
    elif 'Philips' in mfr:
        return EXPECTED_NVOLS.get(('Philips', series_type))
    return None


def _to_int(val):
    """Convert value to int, returning None on failure."""
    try:
        return int(val)
    except (TypeError, ValueError):
        return None


def _check_patient_position(fields):
    """Return (bad, msg) where bad=True if position is not HFS."""
    pos = fields.get('PatientPosition') or ''
    if pos and pos != 'HFS':
        return True, 'PatientPosition is {0}; '.format(pos)
    return False, ''


def _check_matrix(fields, expected):
    """Check AcquisitionMatrixPE against expected value(s). Returns (ok, msg)."""
    val = _to_int(fields.get('AcquisitionMatrixPE'))
    if val is None:
        return False, 'AcquisitionMatrixPE unavailable; '
    if val in expected:
        return True, 'Valid matrix size; '
    return False, 'Invalid matrix size ({0}); '.format(val)


def _check_phase_steps(fields, expected):
    """Check PhaseEncodingSteps against expected value(s). Returns (ok, msg)."""
    val = _to_int(fields.get('PhaseEncodingSteps'))
    if val is None:
        return False, 'PhaseEncodingSteps unavailable; '
    if val in expected:
        return True, 'Valid number of phase encoding steps; '
    return False, 'Invalid number of phase encoding steps ({0}); '.format(val)


def _psn_match(psn, patterns):
    """Return True if psn contains any of the expected patterns."""
    psn_str = psn or ''
    return any(p in psn_str for p in patterns)


def _sd_match(sd, patterns):
    """Return True if SeriesDescription contains any of the expected patterns."""
    sd_str = sd or ''
    return any(p in sd_str for p in patterns)

#############################################################################
# GE compliance
#############################################################################

_GE_PSN_PATTERNS = {
    'T1':      ['research/mprage_promo', 'research/ABCD/mprage_promo',
                'research/ABCD/efgre3d'],
    'T2':      ['research/3dfse_promo', 'research/ABCD/3dfse_promo',
                'research/ABCD/3dfse'],
    'fMRI':    ['research/ABCD/muxepi', 'research/PCGC/muxepi'],
    'fMRI_FM': ['research/ABCD/epi_pepolar', 'research/PCGC/epi_pepolar'],
    'dMRI':    ['research/ABCD/muxepi2', 'research/PCGC/muxepi2'],
    'dMRI_FM': ['research/ABCD/muxepi2', 'research/PCGC/muxepi2'],
}

_GE_SD_PATTERNS = {
    'rsfMRI':          ['rsfMRI'],
    'fMRI_MID_task':   ['Monetary'],
    'fMRI_SST_task':   ['Stop'],
    'fMRI_nBack_task': ['Emotional'],
    'dMRI':            ['DTI', 'ORIG', 'REV'],
    'dMRI_FM':         ['Map'],
}


def check_ge(scan, fields, all_scans):
    """GE-specific compliance checks."""
    series_type = scan['series_type']
    psn = fields.get('PulseSequenceName') or ''
    sd  = fields.get('SeriesDescription') or ''

    # base type for PSN lookup (T1_NORM -> T1, T2_NORM -> T2)
    base_type = series_type.replace('_NORM', '') if '_NORM' in series_type else series_type

    # PulseSequenceNameMatch
    psn_patterns = _GE_PSN_PATTERNS.get(base_type, [])
    psn_match = 'Yes' if (psn_patterns and _psn_match(psn, psn_patterns)) else (
                'No'  if  psn_patterns else 'NA')

    sd_match  = 'NA'
    compliant = 'No'
    info      = 'GE {0}; '.format(series_type)

    if base_type in ('T1', 'T2'):
        ok, msg = _check_matrix(fields, [256])
        info += msg
        compliant = 'Yes' if ok else 'No'
        if series_type in ('T1_NORM', 'T2_NORM'):
            info += 'Normalized image; '

    elif base_type == 'fMRI':
        task_type = _fmri_task_type(scan)
        scan['series_type'] = task_type
        info = 'GE {0}; '.format(task_type)
        ok, msg = _check_matrix(fields, [90])
        info += msg
        compliant = 'Yes' if ok else 'No'
        sd_patterns = _GE_SD_PATTERNS.get(task_type, [])
        if sd_patterns:
            sd_match = 'Yes' if _sd_match(sd, sd_patterns) else 'No'

    elif base_type == 'fMRI_FM':
        pe_dir = scan.get('dir') or ''
        if pe_dir in ('AP', 'PA'):
            scan['series_type'] = 'fMRI_FM_{0}'.format(pe_dir)
            info = 'GE fMRI_FM_{0}; '.format(pe_dir)
        ok, msg = _check_matrix(fields, [90])
        info += msg
        # both dir-AP and dir-PA must exist for same acq (GE integrated sequence)
        acq = scan.get('acq') or ''
        paired = [s for s in all_scans
                  if s['modality'] == 'fmap'
                  and s['suffix'] == 'epi'
                  and (s.get('acq') or '') == acq
                  and s['nii_path'] != scan['nii_path']]
        dirs = set((s.get('dir') or '') for s in paired + [scan])
        if 'AP' in dirs and 'PA' in dirs:
            info += 'Both AP and PA field maps present; '
            compliant = 'Yes' if ok else 'No'
        else:
            info += 'Missing AP or PA field map; '
            compliant = 'No'

    elif base_type in ('dMRI', 'dMRI_FM'):
        if base_type == 'dMRI_FM':
            pe_dir = scan.get('dir') or ''
            if pe_dir in ('AP', 'PA'):
                scan['series_type'] = 'dMRI_FM_{0}'.format(pe_dir)
                info = 'GE dMRI_FM_{0}; '.format(pe_dir)
        ok, msg = _check_matrix(fields, [140])
        info += msg
        compliant = 'Yes' if ok else 'No'
        sd_patterns = _GE_SD_PATTERNS.get(base_type, [])
        if sd_patterns:
            sd_match = 'Yes' if _sd_match(sd, sd_patterns) else 'No'

    pos_bad, pos_msg = _check_patient_position(fields)
    if pos_bad:
        info += pos_msg
        compliant = 'Bad'

    return {
        'series_type':               scan['series_type'],
        'compliant':            compliant,
        'pulse_sequence_name_match': psn_match,
        'series_description_match':  sd_match,
        'info':                      info,
    }

#############################################################################
# Siemens compliance
#############################################################################

_SIEMENS_PSN_PATTERNS = {
    'T1':      ['tfl3d1_16ns', 'tfl_me3d1_16ns'],
    'T2':      ['spc_200ns'],
    'fMRI':    ['epfSM2d1', 'epfid2d1_90'],
    'fMRI_FM': ['epse2d1'],
    'dMRI':    ['ep_b', 'epse2d1_140'],
    'dMRI_FM': ['ep_b', 'epse2d1_140'],
}

_SIEMENS_SD_PATTERNS = {
    'rsfMRI':          ['rsfMRI'],
    'fMRI_MID_task':   ['Monetary'],
    'fMRI_SST_task':   ['Stop'],
    'fMRI_nBack_task': ['Emotional'],
}

_SIEMENS_NTPOINTS = {
    'rsfMRI':          [383],
    'fMRI_MID_task':   [411],
    'fMRI_SST_task':   [445],
    'fMRI_nBack_task': [370],
}


def check_siemens(scan, fields):
    """Siemens-specific compliance checks."""
    series_type = scan['series_type']
    psn = fields.get('PulseSequenceName') or ''
    sd  = fields.get('SeriesDescription') or ''
    image_type = fields.get('ImageType') or ''

    base_type = series_type.replace('_NORM', '') if '_NORM' in series_type else series_type

    psn_patterns = _SIEMENS_PSN_PATTERNS.get(base_type, [])
    psn_match = 'Yes' if (psn_patterns and _psn_match(psn, psn_patterns)) else (
                'No'  if  psn_patterns else 'NA')

    sd_match  = 'NA'
    compliant = 'No'
    info      = 'SIEMENS {0}; '.format(series_type)

    if base_type in ('T1', 'T2'):
        ok_rows, msg_rows = _check_matrix(fields, [256])
        info += msg_rows
        pe_expected = [255, 256] if base_type == 'T1' else [227, 256]
        ok_pe, msg_pe = _check_phase_steps(fields, pe_expected)
        info += msg_pe
        compliant = 'Yes' if (ok_rows and ok_pe) else 'No'
        if series_type in ('T1_NORM', 'T2_NORM'):
            info += 'Normalized image; '

    elif base_type == 'fMRI':
        task_type = _fmri_task_type(scan)
        # also check NumberOfTemporalPositions from JSON as secondary verification
        ntp = _to_int(fields.get('NumberOfTemporalPositions'))
        if ntp is not None and task_type == 'fMRI':
            for t, ntp_list in _SIEMENS_NTPOINTS.items():
                if ntp in ntp_list:
                    task_type = t
                    break
        scan['series_type'] = task_type
        info = 'SIEMENS {0}; '.format(task_type)
        ok_rows, msg_rows = _check_matrix(fields, [90])
        info += msg_rows
        ok_pe, msg_pe = _check_phase_steps(fields, [90])
        info += msg_pe
        compliant = 'Yes' if (ok_rows and ok_pe) else 'No'
        sd_patterns = _SIEMENS_SD_PATTERNS.get(task_type, [])
        if sd_patterns:
            sd_match = 'Yes' if _sd_match(sd, sd_patterns) else 'No'

    elif base_type == 'fMRI_FM':
        # phase encoding direction from BIDS dir entity
        pe_dir = scan.get('dir') or ''
        if pe_dir in ('AP', 'PA'):
            scan['series_type'] = 'fMRI_FM_{0}'.format(pe_dir)
            info = 'SIEMENS fMRI_FM_{0}; '.format(pe_dir)
        ok_rows, msg_rows = _check_matrix(fields, [90])
        info += msg_rows
        ok_pe, msg_pe = _check_phase_steps(fields, [90])
        info += msg_pe
        compliant = 'Yes' if (ok_rows and ok_pe) else 'No'
        if 'NORM' in image_type:
            info += 'Normalized fMRI_FM image; '
            compliant = 'NA'

    elif base_type in ('dMRI', 'dMRI_FM'):
        pe_dir = scan.get('dir') or ''
        if base_type == 'dMRI_FM' and pe_dir in ('AP', 'PA'):
            scan['series_type'] = 'dMRI_FM_{0}'.format(pe_dir)
            info = 'SIEMENS dMRI_FM_{0}; '.format(pe_dir)
        ok_rows, msg_rows = _check_matrix(fields, [140])
        info += msg_rows
        ok_pe, msg_pe = _check_phase_steps(fields, [105, 140])
        info += msg_pe
        compliant = 'Yes' if (ok_rows and ok_pe) else 'No'

    pos_bad, pos_msg = _check_patient_position(fields)
    if pos_bad:
        info += pos_msg
        compliant = 'Bad'

    return {
        'series_type':               scan['series_type'],
        'compliant':            compliant,
        'pulse_sequence_name_match': psn_match,
        'series_description_match':  sd_match,
        'info':                      info,
    }

#############################################################################
# Philips compliance
#############################################################################

_PHILIPS_PSN_PATTERNS = {
    'T1':      ['T1TFE', 'TFE'],
    'T2':      ['TSE'],
    'fMRI':    ['FEEPI'],
    'fMRI_FM': ['SEEPI'],
    'dMRI':    ['DwiSE'],
    'dMRI_FM': ['DwiSE', 'SEEPI'],
}

_PHILIPS_SD_PATTERNS = {
    'rsfMRI':          ['rsfMRI'],
    'fMRI_MID_task':   ['MID', 'Monetary'],
    'fMRI_SST_task':   ['SST', 'Stop'],
    'fMRI_nBack_task': ['Emotion'],
}

_PHILIPS_NTPOINTS = {
    'rsfMRI':          [375, 383],
    'fMRI_MID_task':   [403, 411],
    'fMRI_SST_task':   [437, 445],
    'fMRI_nBack_task': [350, 370],
}


def check_philips(scan, fields):
    """Philips-specific compliance checks."""
    series_type = scan['series_type']
    psn = fields.get('PulseSequenceName') or ''
    sd  = fields.get('SeriesDescription') or ''

    base_type = series_type.replace('_NORM', '') if '_NORM' in series_type else series_type

    psn_patterns = _PHILIPS_PSN_PATTERNS.get(base_type, [])
    psn_match = 'Yes' if (psn_patterns and _psn_match(psn, psn_patterns)) else (
                'No'  if  psn_patterns else 'NA')

    sd_match  = 'NA'
    compliant = 'No'
    info      = 'PHILIPS {0}; '.format(series_type)

    if base_type == 'T1':
        ok, msg = _check_matrix(fields, [240, 256])
        info += msg
        scan_seq = fields.get('ScanningSequence') or ''
        seq_var  = fields.get('SequenceVariant') or ''
        ok_seq = (scan_seq == 'GR')
        ok_var = (seq_var == 'MP')
        info += 'Correct scanning sequence; ' if ok_seq else 'Incorrect scanning sequence; '
        info += 'Correct sequence variant; ' if ok_var else 'Incorrect sequence variant; '
        compliant = 'Yes' if (ok and ok_seq and ok_var) else 'No'
        if series_type == 'T1_NORM':
            info += 'Normalized T1 image; '

    elif base_type == 'T2':
        ok, msg = _check_matrix(fields, [256])
        info += msg
        scan_seq = fields.get('ScanningSequence') or ''
        seq_var  = fields.get('SequenceVariant') or ''
        ok_seq = (scan_seq == 'SE')
        ok_var = (seq_var == 'SK')
        info += 'Correct scanning sequence; ' if ok_seq else 'Incorrect scanning sequence; '
        info += 'Correct sequence variant; ' if ok_var else 'Incorrect sequence variant; '
        compliant = 'Yes' if (ok and ok_seq and ok_var) else 'No'
        if series_type == 'T2_NORM':
            info += 'Normalized T2 image; '

    elif base_type == 'fMRI':
        task_type = _fmri_task_type(scan)
        ntp = _to_int(fields.get('NumberOfTemporalPositions'))
        if ntp is not None and task_type == 'fMRI':
            for t, ntp_list in _PHILIPS_NTPOINTS.items():
                if ntp in ntp_list:
                    task_type = t
                    break
        scan['series_type'] = task_type
        info = 'PHILIPS {0}; '.format(task_type)
        ok, msg = _check_matrix(fields, [89, 90, 92])
        info += msg
        if ntp is not None:
            ntp_expected = _PHILIPS_NTPOINTS.get(task_type, [])
            if ntp_expected and ntp in ntp_expected:
                info += 'Valid number of temporal positions; '
            elif ntp_expected:
                info += 'Invalid number of temporal positions; '
                compliant = 'No'
        compliant = 'Yes' if ok else 'No'
        sd_patterns = _PHILIPS_SD_PATTERNS.get(task_type, [])
        if sd_patterns:
            sd_match = 'Yes' if _sd_match(sd, sd_patterns) else 'No'

    elif base_type == 'fMRI_FM':
        ok, msg = _check_matrix(fields, [89, 90, 92])
        info += msg
        pe_dir = scan.get('dir') or ''
        if pe_dir in ('AP', 'PA'):
            scan['series_type'] = 'fMRI_FM_{0}'.format(pe_dir)
            info = 'PHILIPS fMRI_FM_{0}; '.format(pe_dir) + msg
            compliant = 'Yes' if ok else 'No'
        else:
            info += 'Phase encoding direction could not be established; '

    elif base_type in ('dMRI', 'dMRI_FM'):
        ok, msg = _check_matrix(fields, [140, 141])
        info += msg
        pe_dir = scan.get('dir') or ''
        if base_type == 'dMRI_FM':
            if pe_dir == 'PA':
                scan['series_type'] = 'dMRI_FM_PA'
                info = 'PHILIPS dMRI_FM_PA; ' + msg
            elif pe_dir == 'AP':
                scan['series_type'] = 'dMRI_FM_AP'
                info = 'PHILIPS dMRI_FM_AP; ' + msg
            else:
                info += 'Phase encoding direction could not be established; '
        compliant = 'Yes' if ok else 'No'

    pos_bad, pos_msg = _check_patient_position(fields)
    if pos_bad:
        info += pos_msg
        compliant = 'Bad'

    return {
        'series_type':               scan['series_type'],
        'compliant':            compliant,
        'pulse_sequence_name_match': psn_match,
        'series_description_match':  sd_match,
        'info':                      info,
    }

#############################################################################
# Dispatch
#############################################################################

def check_compliance(scan, fields, all_scans):
    """Dispatch to manufacturer-specific compliance check."""
    series_type = scan['series_type']

    if series_type == 'Undefined':
        return {
            'series_type':               'Undefined',
            'compliant':            'No',
            'pulse_sequence_name_match': 'NA',
            'series_description_match':  'NA',
            'info':                      'Intended ABCD series was not discovered; ',
        }

    manufacturer = (fields.get('Manufacturer') or '').upper()

    if 'GE' in manufacturer:
        return check_ge(scan, fields, all_scans)
    elif 'SIEMENS' in manufacturer:
        return check_siemens(scan, fields)
    elif 'PHILIPS' in manufacturer:
        return check_philips(scan, fields)

    return {
        'series_type':               series_type,
        'compliant':            'NA',
        'pulse_sequence_name_match': 'NA',
        'series_description_match':  'NA',
        'info':                      'Unknown manufacturer: {0}; '.format(manufacturer),
    }

#############################################################################
# Output
#############################################################################

def build_row(pid, sid, scan, fields, result):
    """Build a TSV row dict for one scan."""
    # the session-relative path, as the rawdata and QC scans TSVs spell it
    filename = os.path.join(scan['modality'], os.path.basename(scan['nii_path']))
    series_type  = result['series_type']
    compliant    = result['compliant']
    manufacturer = fields.get('Manufacturer', '')
    sw_versions  = fields.get('SoftwareVersions', '')

    nvols = get_nifti_nvols(scan['nii_path'])
    nvols_str = '' if nvols is None else str(nvols)

    if series_type in THREE_D_SERIES:
        completed = 'Yes'
    else:
        expected = get_expected_nvols(manufacturer, series_type, sw_versions)
        if expected is None or nvols is None:
            completed = 'NA'
        elif nvols < expected:
            completed = 'No'
        else:
            completed = 'Yes'
            if nvols > expected and compliant not in ('No', 'Bad'):
                compliant = 'No'

    return {
        'participant_id':          pid,
        'session_id':              sid,
        'filename':                filename,
        'mr_y_qc__raw__series_num':            fields.get('SeriesNumber', ''),
        'mr_y_qc__raw__series_desc':       fields.get('SeriesDescription', ''),
        'mr_y_qc__raw__series_type':              series_type,
        'mr_y_qc__raw__series__protcom_indicator':               compliant,
        'mr_y_qc__raw__series__psn__match_indicator':  result['pulse_sequence_name_match'],
        'mr_y_qc__raw__series__desc__match_indicator':  result['series_description_match'],
        'mr_y_qc__raw__series__frame_count':                nvols_str,
        'mr_y_qc__raw__series__complete_indicator':               completed,
        'mr_y_qc__raw__series_notes':          result['info'],
        'mr_y_qc__raw__series__dev_manufact':            fields.get('Manufacturer', ''),
        'mr_y_qc__raw__series__dev_model':  fields.get('ManufacturersModelName', ''),
        'mr_y_qc__raw__series__dev__sftw_ver':        fields.get('SoftwareVersions', ''),
        'mr_y_qc__raw__series__dev_serial':      fields.get('DeviceSerialNumber', ''),
        'mr_y_qc__raw__series_psn':       fields.get('PulseSequenceName', ''),
        'mr_y_qc__raw__series__img_type':               fields.get('ImageType', ''),
        'mr_y_qc__raw__series_position':         fields.get('PatientPosition', ''),
        'mr_y_qc__raw__series__acq__matrix_pe':     fields.get('AcquisitionMatrixPE', ''),
        'mr_y_qc__raw__series__pe_steps':      fields.get('PhaseEncodingSteps', ''),
        'mr_y_qc__raw__series__pe__fov_pct':         fields.get('PercentPhaseFOV', ''),
        'mr_y_qc__raw__series_tr':          fields.get('RepetitionTime', ''),
        'mr_y_qc__raw__series_te':                fields.get('EchoTime', ''),
        'mr_y_qc__raw__series__pixel_bw':          fields.get('PixelBandwidth', ''),
        'mr_y_qc__raw__series__study_dt':               fields.get('StudyDate', ''),
        'mr_y_qc__raw__series__study_uid':        fields.get('StudyInstanceUID', ''),
        'mr_y_qc__raw__series__series_uid':       fields.get('SeriesInstanceUID', ''),
        'mr_y_qc__raw__series__series_t':                  fields.get('SeriesTime', ''),
        'mr_y_qc__raw__series__acq_t':             fields.get('AcquisitionTime', ''),
        'mr_y_qc__raw__series__coil_str':                  fields.get('CoilString', ''),
        'mr_y_qc__raw__series__coil_elements':   fields.get('ReceiveCoilActiveElements', ''),
        'mr_y_qc__raw__series__coil_name':             fields.get('ReceiveCoilName', ''),
    }


# the numeric codes the series TSV writes for the categorical columns, the
#   rows holding the words; a series type with no code of its own is not an
#   ABCD series type
YES_NO_NA_CODES = collections.OrderedDict([('Yes', '1'), ('No', '0'),
                                           ('NA', '444')])
SERIES_TYPE_CODES = collections.OrderedDict([
    ('T1', '1'), ('T1_NORM', '2'), ('T2', '3'), ('T2_NORM', '4'),
    ('rsfMRI', '5'), ('fMRI_MID_task', '6'), ('fMRI_SST_task', '7'),
    ('fMRI_nBack_task', '8'), ('fMRI_FM_AP', '9'), ('fMRI_FM_PA', '10'),
    ('dMRI', '11'), ('dMRI_FM_AP', '12'), ('dMRI_FM_PA', '13'),
    ('Undefined', '100')])
SERIES_CODES = {
    'mr_y_qc__raw__series_type': SERIES_TYPE_CODES,
    'mr_y_qc__raw__series__protcom_indicator': collections.OrderedDict([
        ('Yes', '1'), ('No', '0'), ('Bad', '2'), ('NA', '444')]),
    'mr_y_qc__raw__series__psn__match_indicator': YES_NO_NA_CODES,
    'mr_y_qc__raw__series__desc__match_indicator': YES_NO_NA_CODES,
    'mr_y_qc__raw__series__complete_indicator': YES_NO_NA_CODES,
}


def encode_series_row(row):
    """A series row with its categorical words as their numeric codes."""
    out = dict(row)
    for col, codes in SERIES_CODES.items():
        val = out.get(col)
        if val in codes:
            out[col] = codes[val]
        elif col == 'mr_y_qc__raw__series_type' and val:
            out[col] = SERIES_TYPE_CODES['Undefined']
    return out


def decode_series_value(col, val):
    """The word a series TSV code stands for, the value itself if none."""
    val = str(val).strip()
    for word, code in SERIES_CODES.get(col, {}).items():
        if val == code:
            return word
    return val


def na_row(row):
    """A row with its empty values as n/a, as BIDS requires of a TSV."""
    return dict((k, 'n/a' if v is None or str(v).strip() == '' else v)
                for k, v in row.items())


def write_pc_tsv(rows, tsv_path):
    """Write compliance results to a TSV file (one header row + one row per series)."""
    with open(tsv_path, 'w', newline='') as fid:
        writer = csv.DictWriter(fid, fieldnames=TSV_COLUMNS, restval='n/a',
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(na_row(encode_series_row(r)) for r in rows)


def check_fmap_pairs(scans):
    """Session-level fmap pair validation. Returns list of pair result dicts."""
    fmap_scans = [s for s in scans
                  if s['modality'] == 'fmap' and s['suffix'] == 'epi']
    if not fmap_scans:
        return []

    fmap_meta = {}
    for scan in fmap_scans:
        fmap_meta[scan['nii_path']] = load_metadata(scan['json_path'])

    groups = {}
    for scan in fmap_scans:
        meta = fmap_meta[scan['nii_path']]
        key = meta.get('B0FieldIdentifier') or scan.get('acq') or 'unknown'
        groups.setdefault(key, []).append(scan)

    results = []
    for key, group in groups.items():
        fmap_type = 'fmri'
        for scan in group:
            if 'dMRI' in (scan.get('series_type') or ''):
                fmap_type = 'dmri'
                break

        ap_scans = [s for s in group if (s.get('dir') or '') == 'AP']
        pa_scans = [s for s in group if (s.get('dir') or '') == 'PA']

        if len(ap_scans) != 1 or len(pa_scans) != 1:
            results.append({
                'fmap_type': fmap_type,
                'valid':     False,
                'issues':    ['expected 1 AP + 1 PA, found {0} AP + {1} PA'.format(
                               len(ap_scans), len(pa_scans))],
            })
            continue

        ap_meta = fmap_meta[ap_scans[0]['nii_path']]
        pa_meta = fmap_meta[pa_scans[0]['nii_path']]
        issues  = []

        for field in ['SliceThickness', 'ReconMatrixPE',
                      'PhaseEncodingSteps', 'EffectiveEchoSpacing']:
            v_ap = ap_meta.get(field)
            v_pa = pa_meta.get(field)
            if v_ap is not None and v_pa is not None and v_ap != v_pa:
                issues.append('{0} mismatch ({1} vs {2})'.format(field, v_ap, v_pa))

        v_ap = ap_meta.get('FrequencyEncodingSteps')
        v_pa = pa_meta.get('FrequencyEncodingSteps')
        if v_ap is not None and v_pa is not None and v_ap != v_pa:
            issues.append('FrequencyEncodingSteps mismatch ({0} vs {1})'.format(v_ap, v_pa))

        for field in ['ImageOrientationPatientDICOM', 'ShimSetting']:
            v_ap = ap_meta.get(field)
            v_pa = pa_meta.get(field)
            if v_ap is not None and v_pa is not None:
                if isinstance(v_ap, list) and isinstance(v_pa, list):
                    if v_ap != v_pa:
                        issues.append('{0} mismatch'.format(field))
                elif v_ap != v_pa:
                    issues.append('{0} mismatch'.format(field))

        results.append({
            'fmap_type': fmap_type,
            'valid':     len(issues) == 0,
            'issues':    issues,
        })

    return results


# the series column each session summary reads, the rows being the series
#   TSV rows build_row returns
SERIES_FIELD_COLUMN = {
    'SeriesType':        'mr_y_qc__raw__series_type',
    'Compliant':         'mr_y_qc__raw__series__protcom_indicator',
    'Completed':         'mr_y_qc__raw__series__complete_indicator',
    'StudyDate':         'mr_y_qc__raw__series__study_dt',
    'StudyInstanceUID':  'mr_y_qc__raw__series__study_uid',
    'SeriesTime':        'mr_y_qc__raw__series__series_t',
    'SeriesNumber':      'mr_y_qc__raw__series_num',
    'PulseSequenceName': 'mr_y_qc__raw__series_psn',
    'ImageType':         'mr_y_qc__raw__series__img_type',
}


def series_field(row, field):
    """Return a series row's value for field, as '' where it is absent."""
    val = row.get(SERIES_FIELD_COLUMN[field], '')
    return '' if val is None else val


def add_session_event_vars(sess_row, rows):
    """Fill the variables cmig_check_event_pcqc derived from compliance alone.

    The duplicate count compares the fields that identify an acquisition, so
    two rows agreeing on all of them are the same series delivered twice. The
    duration runs from the start of the first series to 5 minutes after the
    start of the last, approximating the last series' own length, in hours,
    from StudyDate and SeriesTime.
    """
    dates = sorted(series_field(r, 'StudyDate') for r in rows
                   if series_field(r, 'StudyDate'))
    stamps = [series_datetime(series_field(r, 'StudyDate'),
                              series_field(r, 'SeriesTime')) for r in rows
              if series_field(r, 'StudyDate')]
    stamps = [t[1] for t in stamps if t[0] == 0]
    sess_row[RAW + 'event_dur'] = (
        round((max(stamps) - min(stamps)).total_seconds() / 3600.0
              + EVENT_DUR_LAST_MIN / 60.0, 4)
        if stamps else 'n/a')

    studies = set(series_field(r, 'StudyInstanceUID') for r in rows
                  if series_field(r, 'StudyInstanceUID'))
    sess_row[RAW + 'event__study__multi_indicator'] = 1 if len(studies) > 1 else 0
    sess_row[RAW + 'event__study__dt__multi_indicator'] = \
        1 if len(set(dates)) > 1 else 0

    # the normalized series count with their plain counterpart, as check_event
    # folds T1_NORM into T1 before counting
    max_dups = 0
    for stem in sorted(set(SERIES_TYPE_TO_STEM.values())):
        key = RAW + '{0}__dup_count'.format(stem)
        if key not in sess_row:
            continue
        stypes = [t for t, v in SERIES_TYPE_TO_STEM.items() if v == stem]
        ids = ['\t'.join([str(series_field(r, f))
                          for f in ('PulseSequenceName', 'StudyDate', 'SeriesTime',
                                    'SeriesNumber', 'ImageType')])
               for r in rows if series_field(r, 'SeriesType') in stypes]
        ndups = len(ids) - len(set(ids))
        sess_row[key] = ndups
        max_dups = max(max_dups, ndups)
    sess_row[RAW + 'event__dup_max'] = max_dups


def _check_num(val):
    """Return val as a number, or None where it is absent or unparsable."""
    try:
        return float(val) if str(val).strip() not in ('', 'n/a', 'NaN', 'nan') else None
    except (TypeError, ValueError):
        return None


def check_series(rows):
    """The series check_event reads, as cmig_check_event_pcqc prepares them.

    Compliance Yes, No or Bad and a defined type, NORM folded into its plain
    type, a study date present, in acquisition order.
    """
    out = []
    for idx, row in enumerate(rows):
        comp = str(series_field(row, 'Compliant')).strip()
        stype = str(series_field(row, 'SeriesType'))
        if comp not in ('Yes', 'No', 'Bad') or 'Undefined' in stype:
            continue
        if not str(series_field(row, 'StudyDate')).strip():
            continue
        done = str(series_field(row, 'Completed')).strip().lower()
        out.append({
            'type': {'T1_NORM': 'T1', 'T2_NORM': 'T2'}.get(stype, stype),
            'comp': comp.lower(),
            'done': 1 if done == 'yes' else 0 if done == 'no' else None,
            'date': str(series_field(row, 'StudyDate')),
            'time': str(series_field(row, 'SeriesTime')),
            'frame': _check_num(row.get('mr_y_qc__raw__series__frame_count', '')),
            'tr': _check_num(row.get('mr_y_qc__raw__series_tr', '')),
            'notes': row.get('mr_y_qc__raw__series_notes', '') or '',
            'filename': row.get('filename', ''),
            'idx': idx,
        })
    return sorted(out, key=lambda r: series_datetime(r['date'], r['time']))


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


def check_slots(series):
    """Number each series by its acquisition, series sharing a time being one."""
    slots, last, n = [], None, -1
    for row in series:
        key = (row['date'], row['time'])
        if key != last:
            n += 1
            last = key
        slots.append(n)
    return slots


def check_fmap(series, slots, pos, stype):
    """The field map pair of the series at pos, which cmig_check_event_pcqc reads.

    Returns fm_missing, fm_wronglocus, fm_invalid, the pair's two filenames,
    its compliance (1 when both are Yes), bad flag (1 when either is Bad),
    notes and completeness, all but the first three None where no pair is
    found.
    """
    ptype = 'fmri' if stype in CHECK_FMRI_STYPES else stype
    types = [r['type'].lower() for r in series]
    fwd = [i for i, t in enumerate(types) if t == ptype + '_fm_pa']
    rev = [i for i, t in enumerate(types) if t == ptype + '_fm_ap']
    fm = {'fm_missing': 1, 'fm_wronglocus': 1, 'fm_invalid': 1,
          'fm_pa_filename': None, 'fm_ap_filename': None, 'fm_pc_score': None,
          'fm_pc_bad': None, 'fm_pc_notes': None, 'fm_complete': None}
    if not (fwd or rev):
        return fm
    fm['fm_missing'] = 0
    here = slots[pos]

    def nearest(positions):
        if not positions:
            return None
        before = [p for p in positions if slots[p] <= here]
        if before:
            return max(before, key=lambda p: slots[p])
        return min(positions, key=lambda p: slots[p])

    f_pos, r_pos = nearest(fwd), nearest(rev)
    if f_pos is None or r_pos is None:
        return fm
    fm['fm_invalid'] = 0
    lo, hi = sorted((slots[f_pos], slots[r_pos]))
    placed = hi <= here and (here - hi) <= 2 and (hi - lo) <= 1
    if placed and lo < hi and hi < here:
        placed = slots[f_pos] < slots[r_pos]
    if placed:
        fm['fm_wronglocus'] = 0
    f_row, r_row = series[f_pos], series[r_pos]
    fm['fm_pa_filename'] = f_row['filename']
    fm['fm_ap_filename'] = r_row['filename']
    fm['fm_pc_score'] = 1 if (f_row['comp'] == 'yes' and r_row['comp'] == 'yes') else 0
    fm['fm_pc_bad'] = 1 if (f_row['comp'] == 'bad' or r_row['comp'] == 'bad') else 0
    fm['fm_pc_notes'] = str(f_row['notes']) + str(r_row['notes'])
    fm['fm_complete'] = 1 if (f_row['done'] and r_row['done']) else 0
    return fm


def add_series_fmap(rows):
    """Add each series' field map pair and its checks to the series rows."""
    series = check_series(rows)
    slots = check_slots(series)
    for i, stype in enumerate(CHECK_STYPES):
        if stype not in CHECK_FMAP_STYPES:
            continue
        for pos, row in enumerate(series):
            if row['type'] != CHECK_SERIES_TYPES[i]:
                continue
            for key, val in check_fmap(series, slots, pos, stype).items():
                rows[row['idx']][SERIES_FM_COLUMN[key]] = val


def check_event(rows):
    """The released variables that need protocol compliance alone.

    The same checks cmig_check_event_pcqc applies, without the QC, so a count
    here equals the one pcqc works out from the same series.
    """
    series = check_series(rows)
    slots = check_slots(series)
    out = {}
    total = {}
    for i, stype in enumerate(CHECK_STYPES):
        stem = CHECK_STEM[stype]
        is_fmap = stype in CHECK_FMAP_STYPES
        c = dict((k, 0) for k, _ in CHECK_COUNT_NAMES + CHECK_FMAP_NAMES)
        comp_sum = {'frame': 0, 'tr': 0}
        comp_n = {'frame': 0, 'tr': 0}
        positions = [p for p, r in enumerate(series)
                     if r['type'] == CHECK_SERIES_TYPES[i]]
        for pos in positions:
            row = series[pos]
            pc_score = 1 if row['comp'] == 'yes' else 0
            pc_bad = 1 if row['comp'] == 'bad' else 0
            completed = row['done']
            c['total_ser'] += 1
            if not pc_score:
                c['ser_failpc'] += 1
            if not completed:
                c['ser_incomp'] += 1
            if is_fmap:
                fm = check_fmap(series, slots, pos, stype)
                fm_pc = fm['fm_pc_score']
                if not fm['fm_missing']:
                    c['fm_wronglocus'] += fm['fm_wronglocus']
                    c['fm_invalid'] += fm['fm_invalid']
                    if not fm_pc:
                        c['fm_failpc'] += 1
                    if not fm['fm_complete']:
                        c['fm_incomp'] += 1
                else:
                    c['fm_missing'] += 1
                    fm_pc = 0
            else:
                fm_pc = 1
            if pc_score and fm_pc:
                c['total_passpc'] += 1
            if completed and not pc_bad and pc_score:
                c['comp_ser'] += 1
                for m in ('frame', 'tr'):
                    if row[m] is not None:
                        comp_sum[m] += row[m]
                        comp_n[m] += 1
        for key, name in CHECK_COUNT_NAMES:
            out['mr_y_qc__raw__{0}__{1}'.format(stem, name)] = c[key]
        if is_fmap:
            for key, name in CHECK_FMAP_NAMES:
                out['mr_y_qc__raw__{0}__{1}'.format(stem, name)] = c[key]
        # the mean over complete and compliant scans divides by their number,
        #   as pcqc does, and is empty where no such scan carries the measure
        for m in CHECK_COMP_MEASURES.get(stype, []):
            out['mr_y_qc__raw__{0}__{1}__comp_mean'.format(stem, m)] = \
                (comp_sum[m] / c['comp_ser']
                 if c['total_ser'] and c['comp_ser'] and comp_n[m] else None)
        total[stype] = c['total_ser']

    # the field map pairs, judged on the field maps' own series, paired in
    #   acquisition order; a pair counts when both are complete and compliant
    types = [r['type'].lower() for r in series]
    for mod in CHECK_PAIR_MODALITIES:
        for d in ('ap', 'pa'):
            out['mr_y_qc__raw__{0}__fm__{1}__series_count'.format(mod, d)] = \
                sum(1 for t in types if t == '{0}_fm_{1}'.format(mod, d))
        fwd = sorted([p for p, t in enumerate(types) if t == mod + '_fm_pa'],
                     key=lambda p: slots[p])
        rev = sorted([p for p, t in enumerate(types) if t == mod + '_fm_ap'],
                     key=lambda p: slots[p])
        n_comp = 0
        for f_pos, r_pos in zip(fwd, rev):
            pair = [series[f_pos], series[r_pos]]
            if all(r['done'] for r in pair) and all(r['comp'] == 'yes' for r in pair):
                n_comp += 1
        out['mr_y_qc__raw__{0}__fm__pass__protcom__comp_count'.format(mod)] = n_comp

    out['mr_y_qc__raw__event__full__protocol_indicator'] = \
        1 if all(total[s] >= CHECK_MIN_NSCANS[i]
                 for i, s in enumerate(CHECK_STYPES)) else 0
    return out


def _first_value(body, col):
    for row in body:
        v = row.get(col, '')
        if v and v.lower() != 'n/a':
            return v
    return None


def _parse_dt(val):
    for fmt in ('%Y-%m-%dT%H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f'):
        try:
            return datetime.datetime.strptime(val, fmt)
        except (TypeError, ValueError):
            continue
    return None


def events_bold(sesdir, fname_events):
    """The rawdata func/ bold of a func/ events file."""
    run_stem = re.sub(r'(_desc-[^_]+)?_events\.tsv$', '',
                      os.path.basename(fname_events))
    return os.path.join(sesdir, 'func', run_stem + '_bold.nii.gz')


def eprime_session_vars(pid, sid, raw_rootdir, exclude_noncompliant_flag=False):
    """The E-Prime session variables of each task, from the rawdata events.

    func/ events are runs with a scan, beh/ events runs without. The scan
    time is the bold json's AcquisitionDateTime (the first stored volume),
    not SeriesTime as the release had it. Old-style
    (release 7.0) events have no timing_ok or start_time, so the valid count
    and the time-based variables are n/a for a task with any old-style func/
    run. match_indicator tests every func/ run, tdiff_min/_max are a true min
    and max over them, and start_t and eprime_dt come from the first run, the
    first func/ run if any, as eprime_missing.csv takes the scanned series
    first.

    A desc-invalid run (wrong trial count, or no responses) counts in
    scan__run_count but not in scan__run__valid_count or match_indicator,
    and start_t and eprime_dt come from valid runs only. A run with both a
    valid file and a stale desc-invalid copy counts once, as the valid one.
    A run whose bold is not usable in the rawdata _scans.tsv (as the run
    selection takes it) does not count in scan__run__valid_count either, nor,
    with exclude_noncompliant_flag, one with compliant 0 there.
    """
    sesdir = os.path.join(raw_rootdir, pid, sid)
    usable = read_scans_usable(sesdir, pid, sid)
    compliant = (read_scans_column(sesdir, pid, sid, 'compliant')
                 if exclude_noncompliant_flag else None)
    out = {}
    for task in EPRIME_TASKS:
        stem = '{0}_{1}_task-{2}'.format(pid, sid, task)
        col = dict((n, EPRIME_COLUMN[n].format(RAW + 'tfmri__{0}__'.format(task)))
                   for n in EPRIME_NAMES)
        runs = {}
        invalid = set()
        for mod in ('func', 'beh'):
            # one file per run, the valid one when a stale desc-invalid copy
            #   is beside it, as the run selection's find_runs() takes it
            by_run = {}
            for f in sorted(glob.glob(os.path.join(sesdir, mod, stem + '_run-*_events.tsv'))):
                m = re.search(r'_run-(\d+)(_desc-invalid)?_events\.tsv$', f)
                if not m:
                    continue
                if m.group(1) in by_run and m.group(2):
                    continue
                by_run[m.group(1)] = f
            runs[mod] = []
            for rnum in sorted(by_run):
                f = by_run[rnum]
                header, body = read_events(f)
                runs[mod].append((f, header, body))
                if 'desc-invalid' in os.path.basename(f):
                    invalid.add(f)
        nfunc, nbeh = len(runs['func']), len(runs['beh'])
        valid_func = [r for r in runs['func'] if r[0] not in invalid]
        old = any(is_old_style(h) for _, h, _ in runs['func'])
        exists = int(bool(nfunc or nbeh))
        vals = {
            'scan__run_count': nfunc,
            'noscan__run_count': nbeh,
            'behav__run_count': sum(
                1 for mod in ('func', 'beh') for _, _, b in runs[mod]
                if 'behavioral' in (_first_value(b, 'task_name') or '').lower()),
            'file_indicator': exists, 'read_indicator': exists,
            'content_indicator': exists, 'corrtask_indicator': exists,
        }
        if old:
            for n in ('scan__run__valid_count', 'match_indicator', 'start_t', 'dt',
                      'dt_indicator', 'tdiff_min', 'tdiff_max'):
                vals[n] = 'n/a'
        else:
            vals['scan__run__valid_count'] = sum(
                1 for f, _, b in valid_func
                if _first_value(b, 'timing_ok') == '1' and
                run_usable(usable, sesdir, events_bold(sesdir, f)) and
                run_compliant(compliant, sesdir, events_bold(sesdir, f)))
            tdiffs = []
            match = int(len(valid_func) > 0)
            for f, _, b in runs['func']:
                start = _parse_dt(_first_value(b, 'start_time'))
                # the bold of this run
                run_stem = re.sub(r'(_desc-[^_]+)?_events\.tsv$', '',
                                  os.path.basename(f))
                acq = bold_acq_datetime(os.path.join(sesdir, 'func',
                                                     run_stem + '_bold.nii.gz'))
                if start is None or acq is None:
                    if f not in invalid:
                        match = 0
                    continue
                d = (start - acq).total_seconds() / 60.0
                tdiffs.append(d)
                if abs(d) > EPRIME_MAX_DELAY and f not in invalid:
                    match = 0
            vals['match_indicator'] = match
            vals['tdiff_min'] = round(min(tdiffs), 2) if tdiffs else 'n/a'
            vals['tdiff_max'] = round(max(tdiffs), 2) if tdiffs else 'n/a'
            first = None
            for mod in ('func', 'beh'):
                for f, _, b in runs[mod]:
                    if f in invalid:
                        continue
                    first = first or _parse_dt(_first_value(b, 'start_time'))
                if first:
                    break
            vals['dt_indicator'] = int(first is not None)
            vals['start_t'] = first.strftime('%H:%M:%S') if first else 'n/a'
            vals['dt'] = first.strftime('%Y-%m-%d') if first else 'n/a'
        for n in EPRIME_NAMES:
            out[col[n]] = vals[n]
    return out


def write_session_pc_tsv(outdir, pid, sid, rows, pair_results, raw_rootdir,
                         exclude_noncompliant_flag=False):
    """Write session-level pc summary to {pid}_{sid}_desc-session_pc.tsv."""
    sess_row = {'participant_id': pid, 'session_id': sid}
    for col in SESSION_PC_COLUMNS:
        if col not in ('participant_id', 'session_id'):
            sess_row[col] = 0

    for row in rows:
        stype = series_field(row, 'SeriesType')
        stem = SERIES_TYPE_TO_STEM.get(stype)
        if stem is None:
            continue
        compliant = series_field(row, 'Compliant')
        completed = series_field(row, 'Completed')
        if stype in ('T1_NORM', 'T2_NORM'):
            sess_row[RAW + stem + '__norm__series_count'] += 1
        elif stem in NOFM_STEMS:
            # the series' own compliance, without the released pass count's
            #   field map condition
            if compliant == 'Yes':
                sess_row[RAW + stem + '__pass__protcom__nofm_count'] += 1
            if completed == 'Yes' and compliant in ('No', 'Bad'):
                sess_row[RAW + stem + '__fail__protcom__comp_count'] += 1
        elif '__fm__' in stem:
            if compliant == 'Yes':
                sess_row[RAW + stem + '__pass__protcom_count'] += 1
            elif compliant in ('No', 'Bad'):
                sess_row[RAW + stem + '__fail__protcom_count'] += 1

    add_session_event_vars(sess_row, rows)

    # E-Prime variables from the rawdata events files
    sess_row.update(eprime_session_vars(pid, sid, raw_rootdir,
                                        exclude_noncompliant_flag))

    # the released variables that need compliance alone, empty ones as n/a
    for key, val in check_event(rows).items():
        sess_row[key] = 'n/a' if val is None else val

    for pair in pair_results:
        ftype = pair['fmap_type']
        if pair['valid']:
            sess_row[RAW + '{0}__fm__pair_count'.format(ftype)] += 1
        else:
            sess_row[RAW + '{0}__fm__pair__invalid_count'.format(ftype)] += 1

    fpath = os.path.join(outdir, '{0}_{1}_desc-session_pc.tsv'.format(pid, sid))
    with open(fpath, 'w', newline='') as fid:
        writer = csv.DictWriter(fid, fieldnames=SESSION_PC_COLUMNS, restval='n/a',
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerow(na_row(sess_row))


def write_info_tsv(outdir, pid, sid, scans):
    """Write session-level scanner metadata to {pid}_{sid}_desc-session_info.tsv."""
    meta = {}
    for mod in sorted(set(s['modality'] for s in scans)):
        for scan in scans:
            if scan['modality'] == mod and scan['json_path']:
                meta = load_metadata(scan['json_path'])
                if meta:
                    break
        if meta:
            break

    info = extract_info_fields(meta)
    row = {
        'participant_id':        pid,
        'session_id':            sid,
        'mr_y_adm__info_sessid': '{0}_{1}'.format(pid, sid),
    }
    row.update(info)

    fpath = os.path.join(outdir, '{0}_{1}_desc-session_info.tsv'.format(pid, sid))
    with open(fpath, 'w', newline='') as fid:
        writer = csv.DictWriter(fid, fieldnames=INFO_COLUMNS, restval='n/a',
                                delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerow(na_row(row))


def write_status(outdir, pid, sid, errcode, errmsg, startdate=None):
    """Write the pc status row and print its block.

    The row and the printed block come from cmig_write_status.py, which the
    MATLAB modules' cmig_write_status.m mirrors. This module printed no block
    before, so cmig_check_batchdir.py could not see a pc job at all.
    """
    fpath = cmig_status.cmig_write_status(
        pid, sid, outdir, 'desc-pc', errcode, errmsg, startdate=startdate)
    if errcode == 0:
        print('{0}: finished'.format(APPNAME))
        print('{0}: completed successfully'.format(APPNAME))
    # pure python with no compiled app, so there is no app_version to record;
    # software_version identifies the code, since it ships in MMPS
    cmig_sing_utils.add_status_versions(fpath, app_name=APPNAME)


def check_prior_status(outdir, pid, sid, forceflag):
    """Return True if a prior successful run exists and forceflag is 0."""
    fpath = cmig_status.status_fpath(outdir, pid, sid, 'desc-pc')
    if not os.path.isfile(fpath):
        return False
    try:
        with open(fpath, 'r') as fid:
            lines = fid.read().splitlines()
        if len(lines) < 2:
            return False
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        status = dict(zip(header, values))
        if int(status.get('errcode', 1)) == 0 and not forceflag:
            return True
    except (IOError, ValueError):
        pass
    return False

#############################################################################
# CLI
#############################################################################

def build_parser():
    parser = argparse.ArgumentParser(
        description='Run ABCD protocol compliance checking on BIDS rawdata.'
    )
    parser.add_argument('--participant_id', required=True,
                        help='Participant ID (e.g., sub-XXXXXXXX)')
    parser.add_argument('--session_id', required=True,
                        help='Session ID (e.g., ses-00A)')
    parser.add_argument('--raw_rootdir', required=True,
                        help='Root directory of raw BIDS data (read-only input)')
    parser.add_argument('--output_rootdir', required=True,
                        help='Root directory for PC output (created if absent)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='Display status messages (default: 0)')
    parser.add_argument('--forceflag', type=int, default=1, choices=[0, 1],
                        help='Overwrite existing output (default: 1)')
    parser.add_argument('--exclude_noncompliant_flag', type=int, default=0,
                        choices=[0, 1],
                        help='scan__run__valid_count also leaves out runs with '
                             'compliant 0 in the rawdata _scans.tsv (default: 0)')
    return parser


def parse_args():
    return build_parser().parse_args()


def run(args):
    pid     = args.participant_id
    sid     = args.session_id
    outdir  = os.path.join(args.output_rootdir, pid, sid)
    startdate = cmig_status.module_start()

    # skip if prior successful run exists
    if check_prior_status(outdir, pid, sid, args.forceflag):
        if args.verbose:
            print('{0}: prior successful run found, skipping'.format(APPNAME))
        return 0, ''

    # validate input
    if not os.path.isdir(args.raw_rootdir):
        errmsg = 'raw_rootdir not found: {0}'.format(args.raw_rootdir)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg

    # find all NIfTI files and JSON sidecars
    scans = get_bids_scans(pid, sid, args.raw_rootdir, verbose=args.verbose)
    if not scans:
        errmsg = 'no scans found for {0} {1}'.format(pid, sid)
        write_status(outdir, pid, sid, 1, errmsg, startdate)
        return 1, errmsg

    # determine series type for each scan from BIDS filename entities
    for scan in scans:
        scan['series_type'] = get_series_type(scan)

    # session-level fmap pair validation (before compliance checks modify series_type)
    pair_results = check_fmap_pairs(scans)

    # run compliance checks and build output rows
    rows = []
    for scan in scans:
        meta   = load_metadata(scan['json_path'])
        fields = extract_fields(meta)
        result = check_compliance(scan, fields, scans)
        rows.append(build_row(pid, sid, scan, fields, result))

    os.makedirs(outdir, exist_ok=True)

    # write session-level info TSV
    write_info_tsv(outdir, pid, sid, scans)

    # write session-level pc TSV
    write_session_pc_tsv(outdir, pid, sid, rows, pair_results, args.raw_rootdir,
                         bool(args.exclude_noncompliant_flag))

    # write series-level pc TSV, with each series' field map pair
    add_series_fmap(rows)
    tsv_fname = '{0}_{1}_desc-series_pc.tsv'.format(pid, sid)
    tsv_path  = os.path.join(outdir, tsv_fname)
    write_pc_tsv(rows, tsv_path)

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
