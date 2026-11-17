#!/usr/bin/env python3
"""Check the output files of the DWI pipeline and write its status file.

Usage:
  cmig_check_bids_dwi_outputs.py --session --in_dir DIR/sub-X/ses-Y [--status] [--startdate T] [--max_rotation DEG]
  cmig_check_bids_dwi_outputs.py --batch --in_dir DIR [--skipped TSV] [--table FILE] [--status_missing]
                       [--max_rotation DEG]

Options:
  --session
      check one session; --in_dir is that session's output folder
  --batch
      check every sub-*/ses-* of an output dataset; --in_dir is the dataset's root (the
      folder with dataset_description.json and the sub-* folders)
  --in_dir DIR
      the session folder (--session) or the output dataset root (--batch)
  --status
      with --session: write dwi/sub-X_ses-Y_desc-proc_status.tsv
  --startdate YYYY-MM-DDTHH:MM:SS.mmm
      with --status: the job's start time, for the status file's startdate column
  --skipped TSV
      with --batch: the skipped.tsv of the batch (cmig_check_bids_dwi_inputs.py, make_dwi_jobs.py);
      its sessions are added to the table as not processable
  --table FILE
      with --batch: write the table here instead of to stdout
  --status_missing
      with --batch: also write a status file for sessions without one; run only after the
      batch has finished
  --max_rotation DEG
      largest rotation in degrees between a series' gradient directions and its raw data
      copy's before the series counts as misregistered; default: 30

--session prints the problems with the session's files and exits 1 if there are any. The
status file has errcode 0 if nothing is wrong, else 1 with the problems as errmsg.
--batch writes one row per session: state ok or failed with the problems, and, with
--skipped, not_processable with the inputs and selection problem from skipped.tsv.

A series whose gradient directions are rotated more than --max_rotation degrees from its
raw data copy in dwi/proc_dwi/ counts as misregistered to the sMRI image. The checks of
dwi/proc_dwi/ (the raw data copy, processing_report.json) are skipped once the folder has
been deleted (cmig_proc_bids_dwi.py, cleanup_flag 1). The bvecs are read as
BIDS does (relative to the image axes, x negated when the affine's determinant is
positive).

Created:  10/05/2026 by Diliana Pecheva
Prev mod: 10/28/2026 by Diliana Pecheva
Prev mod: 11/02/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Last mod: 11/12/2026 by Don Hagler

Renamed from check_dwi_outputs.py (abcd_preproc 81210b9)
"""

import argparse
import datetime
import glob
import json
import os
import re
import sys

import nibabel as nib
import numpy as np

import cmig_sing_utils
from cmig_check_bids_dwi_inputs import INPUTS, not_processed_msg

APP_NAME = "cmig_proc_bids_dwi"
MODULE = "desc-proc"
PROC_DIR = os.path.join("dwi", "proc_dwi")   # a session's intermediate files and logs, by default
SIDECAR_FIELDS = ("PhaseEncodingDirection", "RegistrationMatrixT1w", "B0UnwarpApplied",
                  "SourceVolumeSize", "VolumeSize", "SpatialReference", "TORTOISEVersion")
SERIES_SUFFIX = "_desc-preproc_dwi"


def split_session(session_dir):
    """(sub-X, ses-Y) from a session folder path."""
    sub, ses = os.path.normpath(session_dir).split(os.sep)[-2:]
    return sub, ses


def status_path(session_dir):
    sub, ses = split_session(session_dir)
    return os.path.join(session_dir, "dwi", f"{sub}_{ses}_{MODULE}_status.tsv")


MAX_ROTATION = 30   # degrees, default of --max_rotation


def world_directions(bvecs, affine):
    """Gradient directions (3 x N) in world coordinates from bvecs relative to the image
    axes, as BIDS reads them: x negated first when det(affine[:3, :3]) > 0."""
    v = np.array(bvecs, float).reshape(3, -1)
    m = np.asarray(affine, float)[:3, :3]
    if np.linalg.det(m) > 0:
        v[0] = -v[0]
    return (m / np.linalg.norm(m, axis=0)) @ v


def gradient_rotation(raw_dirs, out_dirs):
    """Angle in degrees of the rotation best taking the raw directions to the output's
    (3 x N, world coordinates), over the volumes with a direction in both. A direction and
    its opposite are the same, so each pair's sign is matched to the fit as it improves."""
    keep = (np.linalg.norm(raw_dirs, axis=0) > 0.5) & (np.linalg.norm(out_dirs, axis=0) > 0.5)
    a = raw_dirs[:, keep] / np.linalg.norm(raw_dirs[:, keep], axis=0)
    b = out_dirs[:, keep] / np.linalg.norm(out_dirs[:, keep], axis=0)
    r = np.eye(3)
    for _ in range(20):
        sign = np.where(np.sum((r @ a) * b, axis=0) < 0, -1.0, 1.0)
        u, _, vt = np.linalg.svd(b @ (a * sign).T)
        r = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
    return float(np.degrees(np.arccos(np.clip((np.trace(r) - 1) / 2, -1, 1))))


def series_files(series):
    """{kind: path} of the files that go with one published series."""
    base = series[:-len(SERIES_SUFFIX + ".nii.gz")]
    return {"bval": base + SERIES_SUFFIX + ".bval", "bvec": base + SERIES_SUFFIX + ".bvec",
            "json": base + SERIES_SUFFIX + ".json",
            "motion": base + "_desc-motion_timeseries.tsv"}


def session_proc_dir(session_dir, proc_dir=None):
    """The session's intermediate folder: proc_dir if given, else dwi/proc_dwi/."""
    return proc_dir if proc_dir else os.path.join(session_dir, PROC_DIR)


def is_philips(series):
    """True if the series' sidecar names Philips as the manufacturer."""
    try:
        with open(series_files(series)["json"]) as fh:
            return "philips" in str(json.load(fh).get("Manufacturer", "")).lower()
    except (OSError, ValueError):
        return False


def check_series(series, b0dx, fmap_epi, max_rotation=MAX_ROTATION, proc_dir=None):
    """Check one published DWI series against its companion files; returns a list of
    problems, empty if none.

    Checks:
    - the DWI series is 4D;
    - bval, bvec, motion TSV and sidecar VolumeSize have one entry per volume;
    - b=0 volume exists;
    - all the sidecar fields in SIDECAR_FIELDS are present;
    - B0dx and field map are on the DWI grid;
    - the gradient rotation from the raw data copy in proc_dir (default dwi/proc_dwi/), if
      it is still there, is at most max_rotation"""
    files = series_files(series)
    try:
        img = nib.load(series)
        if len(img.shape) != 4:
            return [f"DWI is not 4D: {img.shape}"]
        nvol = img.shape[3]
        bvals = np.loadtxt(files["bval"]).ravel()
        bvecs = np.loadtxt(files["bvec"], ndmin=2)
        with open(files["motion"]) as fh:
            motion_rows = [line for line in fh.read().splitlines()[1:] if line.strip()]
        with open(files["json"]) as fh:
            info = json.load(fh)
    except (OSError, ValueError) as err:
        return [f"unreadable output: {err}"]

    p = []
    volume_size = info.get("VolumeSize")
    counts = {"b-values": len(bvals), "bvec columns": bvecs.shape[1],
              "motion TSV rows": len(motion_rows),
              "sidecar VolumeSize[3]": volume_size[3] if isinstance(volume_size, list)
              and len(volume_size) == 4 else None}
    p += [f"{k} = {v}, DWI has {nvol} volumes" for k, v in counts.items() if v != nvol]
    if bvecs.shape[0] != 3:
        p.append(f"bvec has {bvecs.shape[0]} rows, not 3")
    if not (bvals == 0).any():
        p.append("no b-value exactly 0")
    p += [f"sidecar has no {f}" for f in SIDECAR_FIELDS if f not in info]
    if b0dx.shape != img.shape[:3] or not np.allclose(b0dx.affine, img.affine, atol=1e-4):
        p.append("B0dx is not on the DWI grid")
    if fmap_epi.shape[:3] != img.shape[:3] or not np.allclose(fmap_epi.affine, img.affine, atol=1e-4):
        p.append("corrected AP field map is not on the DWI grid")
    if proc_dir is None:
        proc_dir = os.path.join(os.path.dirname(os.path.dirname(series)), PROC_DIR)
    # the raw data copy keeps the rawdata name, without desc-preproc
    raw = os.path.join(proc_dir, os.path.basename(series).replace(SERIES_SUFFIX, "_dwi"))
    if os.path.isfile(raw) and os.path.isfile(raw[:-len(".nii.gz")] + ".bvec"):
        try:
            raw_bvecs = np.loadtxt(raw[:-len(".nii.gz")] + ".bvec", ndmin=2)
            if raw_bvecs.shape == bvecs.shape:
                angle = gradient_rotation(world_directions(raw_bvecs, nib.load(raw).affine),
                                          world_directions(bvecs, img.affine))
                if angle > max_rotation:
                    p.append(f"gradient directions rotated {angle:.0f} degrees from the raw data "
                             f"(limit {max_rotation:g}): likely misregistered to the sMRI image")
        except (OSError, ValueError) as err:
            p.append(f"unreadable raw data copy in {proc_dir}: {err}")
    return p


def session_files(session_dir):
    """(sub-X_ses-Y stem, output root, the session's DWI series) of a session folder."""
    sub, ses = split_session(session_dir)
    stem = f"{sub}_{ses}"
    root = os.path.dirname(os.path.dirname(os.path.normpath(session_dir)))
    series = sorted(glob.glob(os.path.join(session_dir, "dwi", f"{stem}_run-*{SERIES_SUFFIX}.nii.gz")))
    return stem, root, series


def missing_files(session_dir, proc_dir=None):
    """Check if all necessary output files for a session are present. Returns a list of
    problem messages, one per missing or unexpected file; empty if none are missing.

    Necessary: at least one DWI series in dwi/ and, for each, its .bval, .bvec, .json,
    _desc-motion_timeseries.tsv; no other DWI NIfTI in dwi/; one B0dx (with its .json) and
    one corrected AP field map (with its .json) in fmap/; processing_report.json in
    proc_dir (default dwi/proc_dwi/) while that folder exists."""
    stem, root, series = session_files(session_dir)
    dwi_dir = os.path.join(session_dir, "dwi")
    if not series:
        return [f"no {stem}_run-NN{SERIES_SUFFIX}.nii.gz in dwi/"]
    p = []
    names = [os.path.basename(s) for s in series]
    others = [f for f in sorted(os.listdir(dwi_dir)) if re.search(r"dwi.*\.nii", f) and f not in names]
    if others:
        p.append(f"other files in dwi/ the analysis stream would load: {', '.join(others)}")
    files = [f for s in series for f in series_files(s).values()]
    p += [f"missing {os.path.relpath(f, root)}" for f in files if not os.path.isfile(f)]
    proc_dir = session_proc_dir(session_dir, proc_dir)
    report = os.path.join(proc_dir, "processing_report.json")
    if os.path.isdir(proc_dir) and not os.path.isfile(report):
        p.append(f"missing {report}")
    b0dx = glob.glob(os.path.join(session_dir, "fmap", f"{stem}_acq-dwi_desc-AP*PA*_B0dx.nii.gz"))
    if len(b0dx) != 1:
        p.append(f"expected one B0dx in fmap/, found {len(b0dx)}")
    fmap_epi = glob.glob(os.path.join(session_dir, "fmap",
                                      f"{stem}_acq-dwi_dir-AP_run-*_desc-preproc_epi.nii.gz"))
    if len(fmap_epi) != 1:
        p.append(f"expected one corrected AP field map in fmap/, found {len(fmap_epi)}")
    p += [f"missing {os.path.relpath(f.split('.nii')[0] + '.json', root)}"
          for f in b0dx + fmap_epi if not os.path.isfile(f.split(".nii")[0] + ".json")]
    return p


def check_session(session_dir, max_rotation=MAX_ROTATION, proc_dir=None):
    """Check one session's output folder; proc_dir is its intermediate folder (default
    dwi/proc_dwi/), whose checks are skipped once it is deleted. Returns a list of problem
    messages; empty if there are none.

    Checks, in order:
    1. all necessary output files are present (missing_files); if any is missing, returns
    list of missing files;
    2. each DWI series against its companion files (check_series);
    3. the number of series: two for Philips (one per input run), one otherwise;
    4. no absolute paths in processing_report.json, while proc_dir exists."""
    p = missing_files(session_dir, proc_dir)
    if p:
        return p
    stem, root, series = session_files(session_dir)
    proc_dir = session_proc_dir(session_dir, proc_dir)
    report = os.path.join(proc_dir, "processing_report.json")
    try:
        report_info = {}
        if os.path.isfile(report):
            with open(report) as fh:
                report_info = json.load(fh)
        b = nib.load(glob.glob(os.path.join(session_dir, "fmap", f"{stem}_acq-dwi_desc-AP*PA*_B0dx.nii.gz"))[0])
        e = nib.load(glob.glob(os.path.join(session_dir, "fmap",
                                            f"{stem}_acq-dwi_dir-AP_run-*_desc-preproc_epi.nii.gz"))[0])
    except (OSError, ValueError) as err:
        return [f"unreadable output: {err}"]
    for s in series:
        p += [f"{os.path.basename(s)}: {q}" for q in check_series(s, b, e, max_rotation, proc_dir)]
    expected = 2 if is_philips(series[0]) else 1
    if len(series) != expected:
        p.append(f"expected {expected} series in dwi/, found {len(series)}")
    absolute = [k for k, v in report_info.items() if isinstance(v, str) and v.startswith("/")]
    if absolute:
        p.append(f"absolute paths in processing_report.json: {', '.join(absolute)}")
    return p


def write_status(session_dir, errcode, errmsg, startdate=None):
    """Write dwi/sub-X_ses-Y_desc-proc_status.tsv and print the cmig_write_status block."""
    sub, ses = split_session(session_dir)
    fpath = status_path(session_dir)
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    now = datetime.datetime.now()
    procdate = now.strftime("%Y-%m-%dT%H:%M:%S.") + f"{now.microsecond // 1000:03d}"
    errmsg = " ".join(str(errmsg).split())
    header = ["participant_id", "session_id", "procdate", "errcode", "errmsg"]
    row = [sub, ses, procdate, str(errcode), errmsg]
    if startdate:
        header.append("startdate")
        row.append(startdate)
    # an empty cell is written as n/a (BIDS), as cmig_write_status.py writes it
    cells = [c if c else "n/a" for c in row]
    with open(fpath, "w") as fh:
        fh.write("\t".join(header) + "\n" + "\t".join(cells) + "\n")
    # the version columns, as in every status file
    cmig_sing_utils.add_status_versions(fpath, app_name=APP_NAME)
    print("")
    print(f"cmig_write_status: module = {MODULE}")
    if startdate:
        print(f"cmig_write_status: startdate = {startdate}")
    for key, val in (("outdir", os.path.dirname(fpath)), ("procdate", procdate),
                     ("errcode", errcode), ("errmsg", errmsg)):
        print(f"cmig_write_status: {key} = {val}")
    return fpath


def stage_seconds(session_dir):
    """Sum of wall_s in dwi/proc_dwi/benchmark/stages.tsv, or None if no stage was timed."""
    fpath = os.path.join(session_dir, PROC_DIR, "benchmark", "stages.tsv")
    try:
        with open(fpath) as fh:
            rows = [line.split("\t") for line in fh.read().splitlines()[1:] if line.strip()]
        return round(sum(float(r[4]) for r in rows)) if rows else None
    except (OSError, ValueError, IndexError):
        return None


def main():
    # --help prints the module docstring, which describes every option
    ap = argparse.ArgumentParser(description=__doc__, usage=argparse.SUPPRESS, add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--session", action="store_true", help=argparse.SUPPRESS)
    mode.add_argument("--batch", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--in_dir", required=True, metavar="DIR", help=argparse.SUPPRESS)
    ap.add_argument("--status", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--startdate", metavar="YYYY-MM-DDTHH:MM:SS.mmm", help=argparse.SUPPRESS)
    ap.add_argument("--skipped", metavar="TSV", help=argparse.SUPPRESS)
    ap.add_argument("--table", metavar="FILE", help=argparse.SUPPRESS)
    ap.add_argument("--status_missing", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--max_rotation", type=float, default=MAX_ROTATION, metavar="DEG",
                    help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.session and (args.skipped or args.table or args.status_missing):
        ap.error("--skipped, --table and --status_missing go with --batch")
    if args.batch and (args.status or args.startdate):
        ap.error("--status and --startdate go with --session")
    if args.skipped and not os.path.isfile(args.skipped):
        ap.error(f"--skipped: no file {args.skipped}")

    if args.session:
        problems = check_session(args.in_dir, args.max_rotation)
        for q in problems:
            print(f"OUTPUT CHECK: {q}")
        if args.status:
            write_status(args.in_dir, int(bool(problems)), "; ".join(problems), args.startdate)
        return int(bool(problems))

    names = list(INPUTS)
    lines = ["\t".join(["participant_id", "session_id", "state", "stage_wall_s", "problems"] + names
                       + ["selection"])]
    rows = {}
    for s in sorted(glob.glob(os.path.join(args.in_dir, "sub-*", "ses-*"))):
        problems = check_session(s, args.max_rotation)
        if args.status_missing and not os.path.isfile(status_path(s)):
            write_status(s, int(bool(problems)), "; ".join(problems))
        secs = stage_seconds(s)
        rows[split_session(s)] = ["failed" if problems else "ok", "" if secs is None else str(secs),
                                  "; ".join(problems)] + [""] * (len(names) + 1)
    if args.skipped:
        with open(args.skipped) as fh:
            skipped = [line.split("\t") for line in fh.read().splitlines()]
        header = skipped[0]
        for r in skipped[1:]:
            rec = dict(zip(header, r))
            have = {n: rec.get(n) == "1" for n in names}
            rows.setdefault((rec["participant_id"], rec["session_id"]),
                            ["not_processable", "", not_processed_msg(have, rec.get("selection", ""))]
                            + [str(int(have[n])) for n in names] + [rec.get("selection", "")])
    for (sub, ses), r in sorted(rows.items()):
        lines.append("\t".join([sub, ses] + r))
    text = "\n".join(lines) + "\n"
    if args.table:
        with open(args.table, "w") as fh:
            fh.write(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
