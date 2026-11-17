#!/usr/bin/env python3
"""
cmig_classify_abcd_series.py

Classify an ABCD imaging series from the DICOM header of one file.

classify_abcd_series(fname, verbose=False)
  Read one DICOM file and return the series type, an errcode and an
  errmsg. Intended to be called during sourcedata generation, before
  dcm2niix and before any BIDS sidecar exists, so that the type it
  assigns is the source of the series type recorded downstream rather
  than a copy of one taken from elsewhere.

classify_metadata(ds)
  The same decision over a pydicom Dataset already in hand.

classify_undefined(ds)
  The localizers, calibration scans, vNav setters and physio
  recordings, which are not acquisitions at all.

classify_by_description(ds)
  The last resort, for a series the vendor tree leaves undefined.

get_metadata_value(ds, keyword)
  A header value, reaching into the enhanced object's functional
  groups where the top-level element is empty.

Ported from the dispatch blocks of MMPS abcd_pc_GE.m,
abcd_pc_Philips.m and abcd_pc_Siemens.m, which pick one of the
per-type check functions. Only the choice of type is ported here. The
protocol compliance tests those functions then run stay where they
are, in cmig_pc_bids_abcd.py, and so do the file counts and the
InstanceNumber gap test, which need every file of the series where
this reads one.

Written for python 3.6 and pydicom 1.4.

Created:  09/18/2026 by Don Hagler
Last mod: 11/04/2026 by Don Hagler
"""

import os
import struct

import pydicom


# Every type classify_metadata can return. fMRI remains (Siemens,
# Philips) and Undefined_fMRI (GE) only where neither the timing nor the
# description names rsfMRI or a task.
SERIES_TYPES = [
    'T1',
    'T1_NORM',
    'T1_NORM_DIS',
    'T1_Derived',
    'T2',
    'T2_NORM',
    'T2_NORM_DIS',
    'T2_Derived',
    'dMRI',
    'dMRI_FM',
    'dMRI_FM_AP',
    'dMRI_FM_PA',
    'fMRI',
    'fMRI_FM',
    'fMRI_FM_AP',
    'fMRI_FM_PA',
    'rsfMRI',
    'fMRI_MID_task',
    'fMRI_SST_task',
    'fMRI_nBack_task',
    'dMRI_SB',
    'dMRI_FM_SB',
    'dMRI_FM_AP_SB',
    'dMRI_FM_PA_SB',
    'fMRI_SB',
    'Extra_dMRI',
    'Localizer',
    'SMARTPLAN alignment',
    'T2_PROJECTION',
    'Calibration',
    'vNav_setter',
    'PMU',
    'Undefined_dMRI',
    'Undefined_fMRI',
    'Undefined_FM',
    'Undefined',
]

# GE separates a diffusion field map from a diffusion series by the
# number of images in the acquisition, not by the sequence label.
GE_DMRI_FM_IMAGES = (81, 162, 648)
GE_DMRI_IMAGES = (8343, 8424)
GE_NEWVER_IMAGES = (16929, 8464)

# Philips uses the acquisition matrix for both distinctions.
PHILIPS_DMRI_MATRIX = (140, 141)

# GE sequence labels, in Private_0019_109c. muxepi2 is tested before
# muxepi because the second matches the first.
GE_LABELS = [
    (['research/mprage_promo',
      'research/ABCD/mprage_promo',
      'research/ABCD/efgre3d'],       'T1'),
    (['research/3dfse_promo',
      'research/ABCD/3dfse_promo',
      'research/ABCD/3dfse'],         'T2'),
    (['research/ABCD/muxepi2',
      'research/PCGC/muxepi2'],       'dMRI'),
    (['research/ABCD/muxepi',
      'research/PCGC/muxepi'],        'fMRI'),
    (['research/ABCD/epi_pepolar',
      'research/PCGC/epi_pepolar'],   'fMRI_FM'),
]

# Philips sequence types, in Private_2001_1020. T1TFE is tested before
# TFE for the same reason, though both give T1.
PHILIPS_TYPES = [
    ('T1TFE', 'T1'),
    ('TFE',   'T1'),
    ('TSE',   'T2'),
    ('DwiSE', 'dMRI'),
    ('SEEPI', 'FM'),
    ('FEEPI', 'fMRI'),
]

# Philips splits its one field map branch on the acquisition matrix.
PHILIPS_FM_MATRIX = {
    (140, 141): 'dMRI_FM',
    (92, 89):   'fMRI_FM',
}


# The fMRI series each vendor records, once the run is told from the
# task it ran. Siemens reads an acquisition duration string, GE a
# numeric one, Philips the number of temporal positions.
SIEMENS_FMRI_DURATION = [
    ('05:11', 'rsfMRI'),
    ('05:22', 'fMRI_MID_task'),
    ('05:33', 'fMRI_MID_task'),
    ('06:00', 'fMRI_SST_task'),
    ('05:00', 'fMRI_nBack_task'),
]
SIEMENS_FMRI_NTP = {
    383: 'rsfMRI',
    411: 'fMRI_MID_task',
    445: 'fMRI_SST_task',
    370: 'fMRI_nBack_task',
}
GE_FMRI_DURATION = {
    312800384: 'rsfMRI',
    302400384: 'fMRI_nBack_task',
    362400448: 'fMRI_SST_task',
    335200416: 'fMRI_MID_task',
}
PHILIPS_FMRI_NTP = {
    375: 'rsfMRI',      383: 'rsfMRI',
    437: 'fMRI_SST_task', 445: 'fMRI_SST_task',
    350: 'fMRI_nBack_task', 370: 'fMRI_nBack_task',
    403: 'fMRI_MID_task', 411: 'fMRI_MID_task',
}

# The last resort on every vendor, the series description. Philips
# spells three of the four differently.
FMRI_DESCRIPTION = [
    ('rsfMRI', 'rsfMRI'),
    # the E11 rest runs of 2016 are named ABCD_fMRI_rest
    ('rest', 'rsfMRI'),
    ('Monetary', 'fMRI_MID_task'),
    ('Stop', 'fMRI_SST_task'),
    ('Emotional', 'fMRI_nBack_task'),
]
PHILIPS_FMRI_DESCRIPTION = [
    ('rsfMRI', 'rsfMRI'),
    # the E11 rest runs of 2016 are named ABCD_fMRI_rest
    ('rest', 'rsfMRI'),
    ('MID', 'fMRI_MID_task'),
    ('Emotion', 'fMRI_nBack_task'),
    ('SST', 'fMRI_SST_task'),
]

def _getstr(ds, keyword, default=''):
    """Return a header value as a string, or default if absent."""
    val = getattr(ds, keyword, None)
    if val is None:
        return default
    if isinstance(val, (list, tuple)) or type(val).__name__ == 'MultiValue':
        return '\\'.join(str(v) for v in val)
    return str(val)


def _get_private(ds, group, element, default=''):
    """Return a private tag as a string, or default if absent."""
    try:
        val = ds[group, element].value
    except (KeyError, TypeError):
        return default
    if val is None:
        return default
    if isinstance(val, bytes):
        try:
            val = val.decode('utf-8', 'replace')
        except Exception:
            return default
    return str(val)


def _acquisition_matrix(ds):
    """Return the non-zero entries of AcquisitionMatrix, as ints."""
    val = getattr(ds, 'AcquisitionMatrix', None)
    if val is None:
        return []
    try:
        return [int(v) for v in val if int(v) != 0]
    except (TypeError, ValueError):
        return []


def get_image_type(ds):
    """ImageType, reaching into the Siemens per-frame private tags.

    An enhanced multiframe object carries the generic
    ORIGINAL\\PRIMARY\\FMRI\\NONE at the top level and the informative
    value inside PerFrameFunctionalGroupsSequence, which is what
    separates a diffusion series from its field map and one fMRI
    branch from another.
    """
    frames = getattr(ds, 'PerFrameFunctionalGroupsSequence', None)
    if frames is not None and len(frames) > 0:
        for outer, inner in ((0x0021, 0x11fe), (0x0021, 0x10fe)):
            try:
                seq = frames[0][outer, inner].value
            except (KeyError, TypeError):
                continue
            if seq is None or len(seq) == 0:
                continue
            for tag in ((0x0021, 0x1175), (0x0021, 0x1075)):
                val = _get_private(seq[0], tag[0], tag[1])
                if val:
                    return val
    return _getstr(ds, 'ImageType')



# Values the enhanced multiframe object moves out of the top level and
# into a functional group, as the accessors of abcd_pc_Siemens.m reach
# for them. get_image_type above is the one the classification itself
# needs; these are for whatever writes the sidecar, the protocol
# compliance checks being their readers. PercentPhaseFieldOfView is
# read from its own element of the MR FOV/Geometry macro, where
# get_PercentPhaseFieldOfView in abcd_pc_Siemens.m reads the phase
# encoding steps.
_SHARED_GROUPS = {
    'NumberOfPhaseEncodingSteps':
        ('MRFOVGeometrySequence', 'MRAcquisitionPhaseEncodingStepsInPlane'),
    'PercentPhaseFieldOfView':
        ('MRFOVGeometrySequence', 'PercentPhaseFieldOfView'),
    'RepetitionTime':
        ('MRTimingAndRelatedParametersSequence', 'RepetitionTime'),
    'PixelBandwidth':
        ('MRImagingModifierSequence', 'PixelBandwidth'),
}


def _from_group(frames, group, keyword):
    """One value out of a functional group, or None."""
    if frames is None or len(frames) == 0:
        return None
    seq = getattr(frames[0], group, None)
    if seq is None or len(seq) == 0:
        return None
    return getattr(seq[0], keyword, None)


def get_metadata_value(ds, keyword):
    """A header value, reaching into the functional groups if need be.

    The enhanced object leaves the top-level element empty and puts
    the value in a shared or per-frame functional group, which is why
    a sidecar written from the top level alone loses it.
    """
    value = getattr(ds, keyword, None)
    if value not in (None, ''):
        return value
    if keyword == 'ImageType':
        return get_image_type(ds)
    if keyword in _SHARED_GROUPS:
        group, inner = _SHARED_GROUPS[keyword]
        return _from_group(
            getattr(ds, 'SharedFunctionalGroupsSequence', None),
            group, inner)
    if keyword == 'EchoTime':
        return _from_group(
            getattr(ds, 'PerFrameFunctionalGroupsSequence', None),
            'MREchoSequence', 'EffectiveEchoTime')
    if keyword == 'AcquisitionDuration':
        return _get_private(ds, 0x0051, 0x100a) or None
    if keyword == 'AcquisitionTime':
        stamp = _getstr(ds, 'AcquisitionDateTime')
        return stamp[8:] if len(stamp) > 8 else None
    return None

def _sequence_name(ds):
    """SequenceName, falling back to PulseSequenceName."""
    name = _getstr(ds, 'SequenceName')
    if not name:
        name = _getstr(ds, 'PulseSequenceName')
    return name




def _pepolar_enhanced(ds):
    """Phase-encode polarity from the per-frame private tags."""
    frames = getattr(ds, 'PerFrameFunctionalGroupsSequence', None)
    if frames is None or len(frames) == 0:
        return None
    for outer, tag in ((0x11fe, 0x111c), (0x10fe, 0x101c)):
        try:
            seq = frames[0][0x0021, outer].value
        except (KeyError, TypeError):
            continue
        if seq is None or len(seq) == 0:
            continue
        val = _get_private(seq[0], 0x0021, tag)
        if val:
            try:
                return int(float(val))
            except ValueError:
                return None
    return None


def _csa_value(block, key):
    """One named value out of a Siemens CSA2 header block.

    The block is a run of tag records: a 64-byte name, then vm, a
    4-byte vr, syngodt, the item count and a check field, each a
    little-endian int, then that many items, each four ints of which
    the second is the byte length, followed by the bytes padded out to
    a multiple of four. Walking it is what keeps a name from matching
    a byte that merely looks like its value.
    """
    pos = block.find(key + b'\x00')
    if pos < 0:
        return None
    head = pos + 64
    if head + 20 > len(block):
        return None
    nitems = struct.unpack_from('<i', block, head + 12)[0]
    offset = head + 20
    for _ in range(max(nitems, 0)):
        if offset + 16 > len(block):
            return None
        length = struct.unpack_from('<i', block, offset + 4)[0]
        offset += 16
        if length <= 0:
            offset += (4 - length % 4) % 4
            continue
        if offset + length > len(block):
            return None
        text = block[offset:offset + length].split(b'\x00')[0].strip()
        if text:
            return text.decode('latin-1')
        offset += length + (4 - length % 4) % 4
    return None


def _pepolar_csa(ds):
    """The same from PhaseEncodingDirectionPositive in the CSA header.

    The value sits in the Siemens CSA image header at (0029,1010),
    which MMPS parses through mmil_spm_dicom_headers.
    """
    try:
        block = ds[0x0029, 0x1010].value
    except (KeyError, TypeError):
        return None
    if isinstance(block, str):
        block = block.encode('latin-1', 'replace')
    value = _csa_value(block, b'PhaseEncodingDirectionPositive')
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _pe_string_siemens(ds):
    """'AP', 'PA' or None, as mmil_pepolar_Siemens returns them."""
    pepolar = _pepolar_enhanced(ds)
    if pepolar is None:
        pepolar = _pepolar_csa(ds)
    if pepolar == 0:
        return 'PA'
    if pepolar == 1:
        return 'AP'
    return None


def _pe_string_philips(ds):
    """The same from the series description, as check_FM_Philips does.

    An uppercase P names the posterior-to-anterior acquisition and an
    uppercase A the reverse, so 'DTI Fieldmap P' is PA. The test is
    case-sensitive: the lowercase p of Fieldmap must not match.
    """
    description = _getstr(ds, 'SeriesDescription')
    if 'P' in description:
        return 'PA'
    if 'A' in description:
        return 'AP'
    return None



def _mark_anat(series_type, ds, manufacturer):
    """Append the normalised and derived variants of a T1 or T2.

    Siemens says so in the image type, DERIVED\\SECONDARY for a
    derived image and NORM with or without DIS for the rest; GE says
    it in ScanOptions as FILTERED_GEMS. Philips distinguishes
    neither.
    """
    if 'GE' in manufacturer:
        if 'FILTERED_GEMS' in _getstr(ds, 'ScanOptions'):
            return series_type + '_NORM'
        return series_type
    if 'SIEMENS' not in manufacturer.upper():
        return series_type
    image_type = str(get_image_type(ds))
    if 'DERIVED' in image_type and 'SECONDARY' in image_type:
        return series_type + '_Derived'
    if 'NORM' in image_type:
        if 'DIS' in image_type:
            return series_type + '_NORM_DIS'
        return series_type + '_NORM'
    return series_type

def _mark_pe(series_type, pe_string):
    """Append the phase encoding direction where one was found."""
    if pe_string:
        return series_type + '_' + pe_string
    return series_type

def _fmri_by_description(ds, patterns):
    """Name an fMRI series from its description, or None."""
    description = _getstr(ds, 'SeriesDescription')
    for pattern, series_type in patterns:
        if pattern in description:
            return series_type
    return None


def _refine_fmri_siemens(ds):
    """rsfMRI or one of the three tasks, for a Siemens fMRI series."""
    duration = _get_private(ds, 0x0051, 0x100a)
    for pattern, series_type in SIEMENS_FMRI_DURATION:
        if pattern in duration:
            return series_type
    if not duration or duration == 'X':
        try:
            ntp = int(getattr(ds, 'NumberOfTemporalPositions', None))
        except (TypeError, ValueError):
            ntp = None
        if ntp in SIEMENS_FMRI_NTP:
            return SIEMENS_FMRI_NTP[ntp]
    return _fmri_by_description(ds, FMRI_DESCRIPTION) or 'fMRI'


def _refine_fmri_ge(ds):
    """The same for GE, off its numeric acquisition duration."""
    try:
        duration = int(_get_private(ds, 0x0019, 0x105a) or 0)
    except (TypeError, ValueError):
        duration = 0
    if duration in GE_FMRI_DURATION:
        return GE_FMRI_DURATION[duration]
    # check_fMRI_GE in abcd_pc_GE.m leaves an unnamed series Undefined_fMRI
    return _fmri_by_description(ds, FMRI_DESCRIPTION) or 'Undefined_fMRI'


def _refine_fmri_philips(ds):
    """The same for Philips, off the number of temporal positions.

    The description fallback follows what check_fMRI_Philips means
    rather than what it does: there the guard tests Undefined_fMRI,
    which that function never sets, so the branch its own comment
    introduces cannot run. Over the merged table it would have named
    one series in some 131000, the temporal positions resolving
    everywhere else.
    """
    try:
        ntp = int(getattr(ds, 'NumberOfTemporalPositions', None))
    except (TypeError, ValueError):
        ntp = None
    if ntp in PHILIPS_FMRI_NTP:
        return PHILIPS_FMRI_NTP[ntp]
    return _fmri_by_description(ds, PHILIPS_FMRI_DESCRIPTION) or 'fMRI'


# Series whose description marks them as not an ABCD acquisition at
# all, from the data reduction at the top of abcd_classify_incoming.m
# and its SBRef and derived-image guards. Matched case-insensitively.
# delrec names a delayed reconstruction, which arrives beside the
# acquisition it was made from and would otherwise be counted twice.
DESCRIPTION_EXCLUDE = [
    'patient', 'setter', '<mpr collection>', 'loc',
    'sbref', 'dmri_pa_fa', 'dmri_pa_colfa', 'delrec',
]

# Description patterns, in the order get_SeriesType in
# abcd_classify_incoming.m applies them. Each entry is a list of
# substrings any of which matches, a list that must all be absent,
# and the type. Order matters: the diffusion tests come before the
# task ones, and the plain fMRI field map tests come last.
DESCRIPTION_RULES = [
    (['Sag_MPRAGE_T1', 'ABCD_T1', '3D T1'], [], 'T1'),
    (['ABCD_T2', 'Sag_CUBE_T2', '3D T2'], [], 'T2'),
    (['rsfMRI', 'rest'], ['Field'], 'rsfMRI'),
    (['rsfMRI'], [], 'fMRI_FM'),
    (['REV_POL: DTI'], [], 'dMRI_FM'),
    (['DTI 1', 'DTI 2', 'ABCD_dMRI'], ['DistortionMap', 'TRACE'], 'dMRI'),
    (['DTI Fieldmap P', 'dMRI_DistortionMap_P'], [], 'dMRI_FM_PA'),
    (['DTI Fieldmap A', 'dMRI_DistortionMap_A'], [], 'dMRI_FM_AP'),
    (['n-back', 'nback'], [], 'fMRI_nBack_task'),
    (['MID', 'Monetary_Incentive'], [], 'fMRI_MID_task'),
    (['SST', 'task_Stop', 'Stop_Task_Run'], [], 'fMRI_SST_task'),
    (['DistortionMap_P', 'Fieldmap P'], [], 'fMRI_FM_PA'),
    (['DistortionMap_A', 'Fieldmap A'], [], 'fMRI_FM_AP'),
]



def classify_undefined(ds):
    """Name a series the trees reject outright, or None.

    check_undefined in each abcd_pc_*.m names four kinds of series
    that are not an acquisition: a localizer, a calibration scan, the
    vNav setter and a physio recording. They are tested before the
    description is read for anything else, which is what stops
    "Cal:rsfMRI_Run_#1" being taken for a resting state run.
    """
    description = _getstr(ds, 'SeriesDescription')
    lowered = description.lower()
    if 'loc' in lowered:
        return 'Localizer'
    if 'calibration' in lowered or 'Cal:' in description:
        return 'Calibration'
    if 'setter' in lowered:
        return 'vNav_setter'
    if 'pmu' in lowered or 'physiolog' in lowered:
        return 'PMU'
    return None

def classify_by_description(ds):
    """Name a series the vendor tree left undefined, or None.

    The tree reads the sequence and the image type, and a series
    whose acquisition matrix is empty or whose sequence tag is
    unrecognised reaches none of it. Its description usually still
    says what it is. The patterns are abcd_classify_incoming.m's,
    which classifies incoming data the same way, extended with the
    Philips spelling of a field map, "fMRI Fieldmap P", that
    function knowing only the Siemens DistortionMap form.

    A projection image is excluded as well. abcd_pc_Philips.m tests
    PROJECTION IMAGE in its T1 and T2 branches, refusing compliance
    and typing one T2_PROJECTION, so a reformat is already held apart
    from the volume it was made from; the description alone cannot
    tell them apart.
    """
    description = _getstr(ds, 'SeriesDescription').strip()
    if not description:
        return None
    lowered = description.lower()
    for token in DESCRIPTION_EXCLUDE:
        if token in lowered:
            return None
    image_type = str(get_image_type(ds))
    if 'DERIVED' in image_type or 'PROJECTION IMAGE' in image_type:
        return None
    for present, absent, series_type in DESCRIPTION_RULES:
        if not any(token in description for token in present):
            continue
        if any(token in description for token in absent):
            continue
        if series_type.startswith('fMRI_FM') and 'fMRI' not in description:
            continue
        return series_type
    return None

def classify_ge(ds):
    """Series type from the GE sequence label.

    An empty acquisition matrix does not stop the dispatch.
    abcd_pc_GE.m calls check_undefined for one but does not return,
    so the sequence label is still read and its answer overwrites the
    json check_undefined just wrote.
    """
    label = _get_private(ds, 0x0019, 0x109c)
    if not label:
        return 'Undefined'
    for patterns, series_type in GE_LABELS:
        for pattern in patterns:
            if pattern in label:
                if series_type == 'dMRI':
                    return _classify_ge_diffusion(ds)
                if series_type == 'fMRI':
                    return _refine_fmri_ge(ds)
                if series_type in ('T1', 'T2'):
                    return _mark_anat(series_type, ds, 'GE')
                return series_type
    return 'Undefined'


def _classify_ge_diffusion(ds):
    """Split a GE diffusion series from its field map.

    The sequence label is muxepi2 for both, so the count of images in
    the acquisition is what separates them. On software other than
    DV25 or DV26 two counts are ambiguous and REV in the series
    description resolves them. As in abcd_pc_GE.m, any version other
    than DV25 or DV26, an empty one included, counts as new; ABCD has
    no GE version older than DV25.
    """
    try:
        images = int(getattr(ds, 'ImagesInAcquisition', None))
    except (TypeError, ValueError):
        return 'Undefined_dMRI'
    software = _getstr(ds, 'SoftwareVersion')
    if not software:
        software = _getstr(ds, 'SoftwareVersions')
    newver = 'DV25' not in software and 'DV26' not in software
    if newver and images in GE_NEWVER_IMAGES:
        if 'REV' in _getstr(ds, 'SeriesDescription'):
            images = 81
        else:
            images = 8424
    description = _getstr(ds, 'SeriesDescription')
    if newver and 'ORIG' not in description and 'REV' not in description:
        # an extra copy with online correction, which abcd_pc_GE.m
        # marks NA to be ignored
        return 'Extra_dMRI'
    if images in GE_DMRI_FM_IMAGES:
        return 'dMRI_FM'
    if images in GE_DMRI_IMAGES:
        return 'dMRI'
    # the counts vary by software line, so fall back on the description
    if 'Map' in description:
        return 'dMRI_FM'
    if 'DTI' in description:
        return 'dMRI'
    return 'Undefined_dMRI'


def classify_philips(ds):
    """Series type from the Philips sequence type."""
    matrix = _acquisition_matrix(ds)
    if not matrix:
        return 'Undefined'
    seq_type = _get_private(ds, 0x2001, 0x1020)
    if not seq_type:
        return 'Undefined'
    for pattern, series_type in PHILIPS_TYPES:
        if pattern in seq_type:
            if series_type == 'FM':
                if len(matrix) < 2:
                    return 'Undefined_FM'
                key = (matrix[0], matrix[1])
                fm = PHILIPS_FM_MATRIX.get(key)
                if fm is None:
                    return 'Undefined_FM'
                return _mark_pe(fm, _pe_string_philips(ds))
            if series_type == 'dMRI':
                return _classify_philips_diffusion(ds, matrix)
            if series_type == 'fMRI':
                return _refine_fmri_philips(ds)
            if series_type == 'T1':
                return _classify_philips_t1(ds, matrix)
            if series_type == 'T2':
                return _classify_philips_t2(ds)
            return series_type
    return 'Undefined'


def _classify_philips_t1(ds, matrix):
    """Split a Philips T1 from the localizers that share its sequence.

    check_T1_Philips calls a series T1 only on four flags together: a
    256 or 256x240 acquisition matrix, ScanningSequence GR,
    SequenceVariant MP and no projection image. A localizer and the
    SMARTPLAN alignment come through the same T1TFE sequence tag and
    fail those, and the description is what names them. The comment
    in that function says as much.
    """
    ok = (all(v == 256 for v in matrix) or
          (len(matrix) >= 2 and matrix[0] == 256 and matrix[1] == 240))
    if (ok and _getstr(ds, 'ScanningSequence') == 'GR'
            and _getstr(ds, 'SequenceVariant') == 'MP'
            and 'PROJECTION IMAGE' not in str(get_image_type(ds))):
        return _mark_anat('T1', ds, 'Philips')
    description = _getstr(ds, 'SeriesDescription')
    if description.strip() == '3D T1':
        return 'T1'
    if 'Loc' in description or 'loc' in description or 'LOC' in description:
        return 'Localizer'
    if 'SMARTPLAN' in description:
        return 'SMARTPLAN alignment'
    return 'Undefined'


def _classify_philips_t2(ds):
    """A Philips T2, or the projection image made from it.

    check_T2_Philips types a projection T2_PROJECTION rather than
    leaving it undefined, where the T1 branch only refuses it
    compliance.
    """
    if 'PROJECTION IMAGE' in str(get_image_type(ds)):
        return 'T2_PROJECTION'
    return 'T2'

def _classify_philips_diffusion(ds, matrix):
    """Split a Philips diffusion series from its field map.

    Both arrive as DwiSE. The acquisition matrix has to be the
    diffusion one, and Fieldmap in the series description is what
    marks the field map.
    """
    if len(matrix) < 2 or (matrix[0], matrix[1]) != PHILIPS_DMRI_MATRIX:
        return 'Undefined_dMRI'
    if 'Fieldmap' in _getstr(ds, 'SeriesDescription'):
        return _mark_pe('dMRI_FM', _pe_string_philips(ds))
    return 'dMRI'


def classify_siemens(ds):
    """Series type from the Siemens sequence name and image type.

    epse2d1_140 is the syngo MR XA30 sequence for both a diffusion
    series and its field map, where the ep_b names are gone. MB in the
    image type separates them, marking the multiband acquisition as
    MOSAIC does on E11; MMPS counts the files instead, which a
    one-file read has no way to do and which an SBRef defeats anyway,
    its NumberOfTemporalPositions being its parent series'.
    """
    description = _getstr(ds, 'SeriesDescription')
    sequence = _sequence_name(ds)
    image_type = get_image_type(ds)

    vnav = 'vNav_setter' in description
    t1_names = ('tfl3d1_16', 'tfl_me3d1_16', 'tfl_vn3d1_176')
    t2_names = ('spc_200', 'spc_96')

    if not vnav and any(n in sequence for n in t1_names):
        return _mark_anat('T1', ds, 'SIEMENS')
    if not vnav and any(n in sequence for n in t2_names):
        return _mark_anat('T2', ds, 'SIEMENS')

    ep_b = 'ep_b' in sequence
    mosaic = 'MOSAIC' in image_type
    trace = 'TRACE' in image_type
    multiband = 'MB' in image_type
    diffusion = 'DIFFUSION' in image_type
    epse140 = 'epse2d1_140' in sequence

    if (ep_b and not mosaic and not trace) or (epse140 and not multiband):
        named = _mark_pe('dMRI_FM', _pe_string_siemens(ds))
        return _mark_sb(named, multiband, diffusion)
    if (ep_b and mosaic) or epse140:
        return _mark_sb('dMRI', multiband, diffusion)
    if 'epse2d1' in sequence:
        return _mark_pe('fMRI_FM', _pe_string_siemens(ds))
    if 'epfSM2d1' in sequence or 'epfid2d1_90' in sequence:
        marked = _mark_sb('fMRI', multiband, mosaic)
        if marked != 'fMRI':
            return marked
        return _refine_fmri_siemens(ds)
    return 'Undefined'


def _mark_sb(series_type, multiband, accompanying):
    """Append _SB where the series is a single-band reference.

    Neither the multiband marker nor the one that goes with it --
    DIFFUSION on the diffusion branches, MOSAIC on fMRI -- means the
    reference image acquired beside the run rather than the run.
    """
    if not multiband and not accompanying:
        return series_type + '_SB'
    return series_type


def classify_metadata(ds):
    """Series type for a pydicom Dataset, dispatched by manufacturer.

    Where the vendor tree reaches no answer the description is read,
    which recovers a series the tree could not see -- one whose
    acquisition matrix is empty, say.
    """
    manufacturer = _getstr(ds, 'Manufacturer').upper()
    # GE as a whole word, so it cannot match inside another name
    if 'GE' in manufacturer.split():
        series_type = classify_ge(ds)
    elif 'SIEMENS' in manufacturer:
        series_type = classify_siemens(ds)
    elif 'PHILIPS' in manufacturer:
        series_type = classify_philips(ds)
    else:
        series_type = 'Undefined'
    if series_type.startswith('Undefined'):
        named = classify_undefined(ds)
        if named:
            return named
        return classify_by_description(ds) or series_type
    return series_type


def classify_abcd_series(fname, verbose=False):
    """Classify a series from one of its DICOM files.

    Parameters
    ----------
    fname : str
        Path to any one DICOM file of the series.
    verbose : bool

    Returns
    -------
    series_type : str
        One of SERIES_TYPES.
    errcode : int
        0 on success, 1 on error.
    errmsg : str
    """
    if not os.path.exists(fname):
        return 'Undefined', 1, 'file not found: {}'.format(fname)
    try:
        ds = pydicom.dcmread(fname, stop_before_pixels=True)
    except Exception as exc:
        return 'Undefined', 1, 'cannot read {}: {}'.format(fname, exc)

    series_type = classify_metadata(ds)
    if verbose:
        print('classify_abcd_series: {} -> {}'.format(fname, series_type))
    return series_type, 0, ''
