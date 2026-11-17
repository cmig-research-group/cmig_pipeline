#!/usr/bin/env python3
"""Check whether sessions can be processed by the DWI pipeline.

Usage:
  cmig_check_bids_dwi_inputs.py --in_dir DIR --smri_dir DIR --session sub-X/ses-Y
  cmig_check_bids_dwi_inputs.py --in_dir DIR --smri_dir DIR --batch [--skipped FILE]

Options:
  --in_dir DIR
      the DWI input directory: BIDS root with the DWI, field maps and scans.tsv
  --smri_dir DIR
      root of the sMRI processing output (anat/ desc-preproc T1w or T2w)
  --session sub-X/ses-Y
      check this session: prints why it cannot be processed and exits 1, else exits 0
  --batch
      check every sub-*/ses-* of the sMRI tree and write the skipped table
  --skipped FILE
      with --batch: write the skipped table here instead of to stdout

A session can be processed when it has the INPUTS (an sMRI T1w or T2w, an input DWI and
an AP field map) and cmig_select_bids_dwi_runs.py can choose its runs. The skipped table has one row
per session that cannot: participant_id, session_id, one 1/0 column per input (1 =
present) and the selection problem. make_dwi_jobs.py writes it as skipped.tsv in the
batch folder.

Created:  10/05/2026 by Diliana Pecheva
Prev mod: 10/06/2026 by Diliana Pecheva
Prev mod: 11/02/2026 by Don Hagler
Last mod: 11/11/2026 by Don Hagler

Renamed from check_dwi_inputs.py (abcd_preproc 81210b9)
"""

import argparse
import glob
import os
import sys

from cmig_select_bids_dwi_runs import NotProcessable, select_runs

# inputs without which the pipeline exits early: name -> (root, glob under root/sub-X/ses-Y)
INPUTS = {
    "t1w_or_t2w": ("smri_dir", ["anat/*_desc-preproc_T1w.nii*", "anat/*_desc-preproc_T2w.nii*"]),
    "dwi": ("in_dir", ["dwi/*_run-*_dwi.nii*"]),
    "reverse_pe": ("in_dir", ["fmap/*_acq-dwi_dir-AP_*_epi.nii*"]),
}
SKIPPED_COLUMNS = ["participant_id", "session_id"] + list(INPUTS) + ["selection"]


def sessions(smri_dir):
    """(sub-X, ses-Y) pairs of the sMRI tree, sorted."""
    dirs = sorted(glob.glob(os.path.join(smri_dir, "sub-*", "ses-*")))
    return [tuple(d.split(os.sep)[-2:]) for d in dirs if os.path.isdir(d)]


def input_presence(in_dir, smri_dir, sub, ses):
    """{input name: True if present} for the INPUTS of one session."""
    roots = {"in_dir": in_dir, "smri_dir": smri_dir}
    return {name: any(glob.glob(os.path.join(roots[root], sub, ses, pat)) for pat in pats)
            for name, (root, pats) in INPUTS.items()}


def selection_problem(in_dir, sub, ses):
    """Why cmig_select_bids_dwi_runs.py cannot choose runs for the session, or "" if it can."""
    try:
        select_runs(os.path.join(in_dir, sub, ses))
        return ""
    except (NotProcessable, OSError, ValueError, KeyError) as err:
        return str(err)


def input_check(in_dir, smri_dir, sub, ses):
    """(input_presence, selection_problem) of one session; the run selection is only
    tried when the DWI and field map are present."""
    have = input_presence(in_dir, smri_dir, sub, ses)
    sel = selection_problem(in_dir, sub, ses) if have["dwi"] and have["reverse_pe"] else ""
    return have, sel


def processable(have, sel):
    return all(have.values()) and not sel


def not_processed_msg(presence, selection=""):
    """'not processed: ...' for absent inputs and a failed run selection, or '' if neither."""
    parts = [f"missing {', '.join(n for n, ok in presence.items() if not ok)}"] \
        if not all(presence.values()) else []
    parts += [selection] if selection else []
    return f"not processed: {'; '.join(parts)}" if parts else ""


def skipped_row(sub, ses, have, sel):
    """One line of the skipped table."""
    return "\t".join([sub, ses] + [str(int(have[k])) for k in INPUTS] + [sel])


def skipped_table(rows):
    """The skipped table's text from skipped_row lines."""
    return "\n".join(["\t".join(SKIPPED_COLUMNS)] + rows) + "\n"


def main():
    # --help prints the module docstring, which describes every option
    ap = argparse.ArgumentParser(description=__doc__, usage=argparse.SUPPRESS, add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    ap.add_argument("--in_dir", required=True, metavar="DIR", help=argparse.SUPPRESS)
    ap.add_argument("--smri_dir", required=True, metavar="DIR", help=argparse.SUPPRESS)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--session", metavar="sub-X/ses-Y", help=argparse.SUPPRESS)
    mode.add_argument("--batch", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--skipped", metavar="FILE", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.skipped and not args.batch:
        ap.error("--skipped goes with --batch")

    if args.session:
        sub, ses = os.path.normpath(args.session).split(os.sep)[-2:]
        have, sel = input_check(args.in_dir, args.smri_dir, sub, ses)
        msg = not_processed_msg(have, sel)
        print(f"{sub} {ses}: {msg or 'can be processed'}")
        return int(bool(msg))

    rows, n = [], 0
    for sub, ses in sessions(args.smri_dir):
        n += 1
        have, sel = input_check(args.in_dir, args.smri_dir, sub, ses)
        if not processable(have, sel):
            rows.append(skipped_row(sub, ses, have, sel))
    if args.skipped:
        with open(args.skipped, "w") as fh:
            fh.write(skipped_table(rows))
        print(f"{n} sessions, {len(rows)} cannot be processed ({args.skipped})")
    else:
        sys.stdout.write(skipped_table(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
