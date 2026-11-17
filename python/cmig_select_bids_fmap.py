#!/usr/bin/env python3
"""Choose the field map for a scan among candidate field maps.

The caller filters the candidates (acquisition type, phase-encode direction, usable).
Order, as cmig_select_bids_fmap.m, with the shim as a check where the sidecars record
ShimSetting (Siemens, GE):
1. label: a candidate whose B0FieldIdentifier equals the scan's B0FieldSource, unless its
   ShimSetting differs from the scan's;
2. shim: a candidate with the scan's ShimSetting;
3. label again when no candidate has the scan's shim (shim_match false);
4. time: else any candidate.
Within a step the candidate nearest in AcquisitionDateTime, before or after the scan, is
taken. Unlike cmig_select_bids_fmap.m, a label is overridden when the shims contradict it.

Created:  10/05/2026 by Diliana Pecheva
Last mod: 11/02/2026 by Don Hagler

Renamed from select_fmap.py (abcd_preproc 81210b9)
"""

import datetime


def acq_seconds(info):
    """AcquisitionDateTime in seconds (seconds may lack zero padding); None if unreadable."""
    day, _, clock = str(info.get("AcquisitionDateTime", "")).partition("T")
    try:
        h, m, s = clock.split(":")
        start = datetime.datetime.strptime(day, '%Y-%m-%d')
        return start.timestamp() + int(h) * 3600 + int(m) * 60 + float(s)
    except ValueError:
        return None


def nearest(scan_info, candidates):
    """The candidate acquired nearest to the scan; the first one if the times are unreadable."""
    t = acq_seconds(scan_info)
    times = [acq_seconds(info) for _, info in candidates]
    if t is None or None in times:
        return candidates[0]
    return min(zip(candidates, times), key=lambda ct: abs(ct[1] - t))[0]


def select_fmap(scan_info, candidates):
    """scan_info: the scan's sidecar dict; candidates: list of (path, sidecar dict).
    Returns (path, method, shim_match): method 'label', 'shim' or 'time'; shim_match
    True/False, or None when the shims cannot be compared. None if there are no candidates."""
    if not candidates:
        return None
    shim = scan_info.get("ShimSetting")
    comparable = shim is not None and any(info.get("ShimSetting") is not None for _, info in candidates)

    def same_shim(info):
        return comparable and info.get("ShimSetting") == shim

    label = scan_info.get("B0FieldSource")
    labelled = [c for c in candidates if label and c[1].get("B0FieldIdentifier") == label]
    trusted = [c for c in labelled if not comparable or same_shim(c[1])]
    if trusted:
        path, info = nearest(scan_info, trusted)
        method = "label"
    elif any(same_shim(info) for _, info in candidates):
        path, info = nearest(scan_info, [c for c in candidates if same_shim(c[1])])
        method = "shim"
    elif labelled:
        path, info = nearest(scan_info, labelled)
        method = "label"
    else:
        path, info = nearest(scan_info, candidates)
        method = "time"
    return path, method, same_shim(info) if comparable else None
