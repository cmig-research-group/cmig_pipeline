#!/usr/bin/env python3
"""Choose the DWI run(s) and reverse-PE field map of one ABCD session from its scans.tsv.

Usage: cmig_select_bids_dwi_runs.py IN_DIR/sub-X/ses-Y

Arguments:
  IN_DIR/sub-X/ses-Y
      the session's folder in the DWI input directory (its scans.tsv, dwi/ and fmap/ with their sidecars)

Prints up_nifti, up_nifti2 (Philips second run, else empty), down_nifti, shim_match
(1/0 when the shims could be compared, else empty) and fmap_method (label, shim or time)
as shell assignments; exits 1 with the reason when the session cannot be processed.
Only runs with usable == 1 are used.
Siemens/GE: the DWI is chosen by cmig_select_scan.m's rule (shortest notes when their
lengths differ, then the latest acq_time); its AP field map by cmig_select_bids_fmap.py.
Philips (no ShimSetting): runs are grouped by B0FieldSource, or, without one, by the latest
AP acquired before them; a pair is exactly two usable runs of one group with different
bvecs, chosen by the same rule (longest notes of the two, later run); its field map is a
usable AP labelled for the group, or for an unlabelled group the AP that defines it.

Created:  10/05/2026 by Diliana Pecheva
Prev mod: 10/06/2026 by Diliana Pecheva
Last mod: 11/02/2026 by Don Hagler

Renamed from select_dwi_runs.py (abcd_preproc 81210b9)
"""

import csv
import glob
import json
import os
import shlex
import sys

import numpy as np

from cmig_select_bids_fmap import select_fmap


class NotProcessable(Exception):
    pass


def read_scans(session_dir):
    """(DWI rows, AP field map rows) of the session's scans.tsv, each with its full path."""
    files = glob.glob(os.path.join(session_dir, "*_scans.tsv"))
    if len(files) != 1:
        raise NotProcessable(f"expected one scans.tsv in {session_dir}, found {len(files)}")
    with open(files[0]) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["path"] = os.path.join(session_dir, r["filename"])
    dwi = [r for r in rows if r["filename"].startswith("dwi/") and "_dwi.nii" in r["filename"]]
    ap = [r for r in rows if r["filename"].startswith("fmap/")
          and "_acq-dwi_dir-AP_" in r["filename"] and "_epi.nii" in r["filename"]]
    return dwi, ap


def sidecar(r):
    with open(r["path"].split(".nii")[0] + ".json") as fh:
        return json.load(fh)


def usable(rows):
    return [r for r in rows if r.get("usable") == "1"]


def pick(cands, notes, time):
    """cmig_select_scan.m: shortest notes when their lengths differ, then the latest time."""
    lens = [len(notes(c)) for c in cands]
    if any(lens) and len(set(lens)) > 1:
        cands = [c for c, n in zip(cands, lens) if n == min(lens)]
    return max(cands, key=time)


def select_siemens_ge(dwi, ap):
    d = pick(usable(dwi), lambda r: r["notes"], lambda r: r["acq_time"])
    uap = usable(ap)
    if not uap:
        raise NotProcessable("no usable AP field map")
    path, method, match = select_fmap(sidecar(d), [(a["path"], sidecar(a)) for a in uap])
    return {"up_nifti": d["path"], "up_nifti2": "", "down_nifti": path,
            "shim_match": "" if match is None else str(int(match)), "fmap_method": method}


def philips_fmap(key, first_run, ap):
    """(path, method) of the usable AP field map of a Philips group, or None."""
    if key.startswith("time:"):
        a = next(a for a in ap if a["filename"] == key[len("time:"):])
        return (a["path"], "time") if a.get("usable") == "1" else None
    cands = [(a["path"], sidecar(a)) for a in usable(ap)]
    cands = [c for c in cands if c[1].get("B0FieldIdentifier") == key]
    return select_fmap(sidecar(first_run), cands)[:2] if cands else None


def select_philips(dwi, ap):
    groups = {}
    for r in dwi:
        key = sidecar(r).get("B0FieldSource")
        if not key:
            before = [a for a in ap if a["acq_time"] < r["acq_time"]]
            if not before:
                continue
            key = "time:" + max(before, key=lambda a: a["acq_time"])["filename"]
        groups.setdefault(key, []).append(r)
    pairs = []
    for key, runs in groups.items():
        runs = sorted(usable(runs), key=lambda r: r["acq_time"])
        fmap = philips_fmap(key, runs[0], ap) if len(runs) == 2 else None
        if fmap:
            pairs.append((runs, fmap))
    if not pairs:
        raise NotProcessable("no two usable Philips DWI runs with a usable AP field map of their own")
    (r1, r2), (down, method) = pick(pairs, lambda p: max((r["notes"] for r in p[0]), key=len),
                                    lambda p: p[0][1]["acq_time"])
    b1, b2 = (np.loadtxt(r["path"].split(".nii")[0] + ".bvec") for r in (r1, r2))
    if b1.shape == b2.shape and np.allclose(b1, b2):
        raise NotProcessable(f"Philips runs {r1['filename']} and {r2['filename']} have the same bvecs")
    return {"up_nifti": r1["path"], "up_nifti2": r2["path"], "down_nifti": down,
            "shim_match": "", "fmap_method": method}


def select_runs(session_dir):
    """{up_nifti, up_nifti2, down_nifti, shim_match} for a session folder of the DWI input directory;
    raises NotProcessable with the reason."""
    dwi, ap = read_scans(session_dir)
    if not usable(dwi):
        raise NotProcessable("no usable DWI run")
    if not ap:
        raise NotProcessable("no AP field map")
    if "philips" in str(sidecar(dwi[0]).get("Manufacturer", "")).lower():
        return select_philips(dwi, ap)
    return select_siemens_ge(dwi, ap)


def main():
    try:
        sel = select_runs(sys.argv[1])
    except (NotProcessable, OSError, ValueError, KeyError) as err:
        print(err)
        return 1
    for k, v in sel.items():
        print(f"{k}={shlex.quote(v)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
