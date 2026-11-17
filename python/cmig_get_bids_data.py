#!/usr/bin/env python3
"""
Find BIDS data files for a session, filtered by modality, series type and task.

Python equivalent of cmig_get_bids_data.m, matching its selection rules so that
python callers and the MATLAB modules see the same set of files. Used to decide
which tasks a session has data for, and available to any python module needing
BIDS discovery.

Selection follows the MATLAB:
  - files under {rootdir}/{pid}/{sid}/{modality}/ matching the entity pattern
  - mask files excluded unless that is all there is
  - files intersected with the usable rows of {pid}_{sid}_scans.tsv, when that
    file exists; a session with no usable scans yields nothing

Created:  08/10/2026 by Don Hagler
Last mod: 08/10/2026 by Don Hagler
"""

import os
import re

FNAME_SCANS = '{0}_{1}_scans.tsv'
MASK_PATTERN = r'_mask\.nii'


def _read_scans(sessdir, participant_id, session_id):
    """Return (usable_filenames, found_flag) from the session scans.tsv.

    usable_filenames holds session-relative names of rows with usable == 1.
    found_flag distinguishes an absent scans.tsv, where no filtering applies,
    from one listing no usable scans, where everything is excluded.

    A scans.tsv without a usable column warns and treats all rows as usable,
    matching cmig_get_bids_data.m.
    """
    fname_scans = os.path.join(sessdir, FNAME_SCANS.format(participant_id, session_id))
    if not os.path.isfile(fname_scans):
        return [], False

    with open(fname_scans) as fid:
        lines = fid.read().splitlines()
    if len(lines) < 2:
        return [], False

    header = lines[0].split('\t')
    if 'filename' not in header:
        return [], False
    i_fname = header.index('filename')
    i_usable = header.index('usable') if 'usable' in header else None
    if i_usable is None:
        print('WARNING: no usable column in {0}; treating all scans as usable'.format(
            fname_scans))

    fnames = []
    for line in lines[1:]:
        cols = line.split('\t')
        if len(cols) <= i_fname:
            continue
        if i_usable is not None:
            if len(cols) <= i_usable or cols[i_usable].strip() != '1':
                continue
        fnames.append(cols[i_fname].strip())
    return fnames, True


def get_bids_data(participant_id, session_id, rootdir, modality=None,
                  series_type=None, task=None, acq_type=None):
    """Return (fnames_data, fnames_info) for one session.

    fnames_data are full paths to matching image files, sorted. fnames_info are
    the corresponding json paths, with None where no json exists.

    modality, series_type, task and acq_type each narrow the match when given.
    """
    sessdir = os.path.join(rootdir, participant_id, session_id)
    if not os.path.isdir(sessdir):
        return [], []

    usable, scans_found = _read_scans(sessdir, participant_id, session_id)

    indir = os.path.join(sessdir, modality) if modality else sessdir
    if not os.path.isdir(indir):
        return [], []

    # build the entity pattern in BIDS order, as cmig_get_bids_data.m does
    pat = ''
    if acq_type:
        pat += '.*acq-{0}'.format(acq_type)
    if task:
        pat += '.*task-{0}'.format(task)
    if series_type:
        pat += '.*{0}'.format(series_type)
    regex = re.compile(pat) if pat else None

    fnames_data = []
    for name in sorted(os.listdir(indir)):
        if not name.endswith(('.nii', '.nii.gz')):
            continue
        if regex is not None and not regex.search(name):
            continue
        fnames_data.append(os.path.join(indir, name))

    if not fnames_data:
        return [], []

    # exclude mask files unless that is all there is
    notmask = [f for f in fnames_data if not re.search(MASK_PATTERN, os.path.basename(f))]
    if notmask:
        fnames_data = notmask

    # exclude unusable files; gated on whether scans.tsv was read, so a session
    # with no usable scans excludes everything rather than nothing
    if scans_found:
        keep = set(os.path.normpath(os.path.join(sessdir, u)) for u in usable)
        fnames_data = [f for f in fnames_data if os.path.normpath(f) in keep]
        if not fnames_data:
            return [], []

    fnames_info = []
    for fname in fnames_data:
        fname_json = re.sub(r'\.nii(\.gz)?$', '.json', fname)
        fnames_info.append(fname_json if os.path.isfile(fname_json) else None)

    return fnames_data, fnames_info


def get_bids_tasks(participant_id, session_id, rootdir, modality='func',
                   series_type='bold'):
    """Return the task labels a session has usable data for, sorted.

    Candidates come from the task entities actually present rather than a fixed
    list, so a new or renamed label is not silently dropped.
    """
    fnames_data, _ = get_bids_data(participant_id, session_id, rootdir,
                                   modality=modality, series_type=series_type)
    tasks = set()
    for fname in fnames_data:
        match = re.search(r'_task-([a-zA-Z0-9]+)', os.path.basename(fname))
        if match:
            tasks.add(match.group(1))
    return sorted(tasks)
