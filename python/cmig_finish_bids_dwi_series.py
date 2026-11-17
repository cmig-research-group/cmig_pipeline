#!/usr/bin/env python3
"""Finish a published DWI series (or field map): TR in the NIfTI header, tidy b-values,
bvecs in the BIDS convention.

Usage: cmig_finish_bids_dwi_series.py --nifti FILE.nii[.gz] --tr SEC [--bval FILE --bvec FILE]

Options:
  --nifti FILE.nii[.gz]
      the series; if 4D, its header gets pixdim[4] = TR and units mm and s (data
      unchanged; .nii.gz is written without name or time, as gzip -n); a 3D image is left
      as it is
  --tr SEC
      RepetitionTime of the series, in seconds
  --bval FILE
      b-values: those below B0_MAX become exactly 0, the others are rounded to whole numbers
  --bvec FILE
      gradient directions relative to the image axes; with --bval. When the image's
      affine has a positive determinant (e.g. LPS), their x components are negated, as
      BIDS requires; run once per series

Exits 1 on any error.

Created:  10/07/2026 by Diliana Pecheva
Last mod: 11/02/2026 by Don Hagler

Renamed from finish_dwi_series.py (abcd_preproc 81210b9)
"""

import argparse
import gzip
import shutil
import sys

import nibabel as nib
import numpy as np

B0_MAX = 5          # b-values below this are written as 0


def set_tr(path, tr):
    """Write pixdim[4] = tr and units mm/s into the header of a 4D image, copying the data
    bytes; True if written, False for a 3D image (left unchanged)."""
    gz = path.endswith(".gz")
    opener = gzip.open if gz else open
    with opener(path, "rb") as fh:
        hdr = nib.Nifti1Header(binaryblock=fh.read(348))
        rest = fh.read()          # extension flag, extensions and data, unchanged
    if hdr["dim"][0] < 4 or hdr["dim"][4] < 2:
        return False
    hdr["pixdim"][4] = tr
    hdr.set_xyzt_units("mm", "sec")
    tmp = path + ".tmp"
    with open(tmp, "wb") as raw:
        if gz:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as out:
                out.write(hdr.binaryblock)
                out.write(rest)
        else:
            raw.write(hdr.binaryblock)
            raw.write(rest)
    shutil.move(tmp, path)
    return True


def tidy_bvals(bvals):
    b = np.asarray(bvals, float)
    return np.where(b < B0_MAX, 0, np.round(b)).astype(int)


def bids_bvecs(bvecs, affine):
    """bvecs (3 x N) relative to the image axes, x negated when det(affine[:3, :3]) > 0."""
    v = np.array(bvecs, float)
    if np.linalg.det(np.asarray(affine)[:3, :3]) > 0:
        v[0] = -v[0]
    return v


def main():
    # --help prints the module docstring, which describes every option
    ap = argparse.ArgumentParser(description=__doc__, usage=argparse.SUPPRESS, add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    ap.add_argument("--nifti", required=True, metavar="FILE", help=argparse.SUPPRESS)
    ap.add_argument("--tr", required=True, type=float, metavar="SEC", help=argparse.SUPPRESS)
    ap.add_argument("--bval", metavar="FILE", help=argparse.SUPPRESS)
    ap.add_argument("--bvec", metavar="FILE", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if bool(args.bval) != bool(args.bvec):
        ap.error("--bval and --bvec go together")
    try:
        set_tr(args.nifti, args.tr)
        if args.bval:
            bvals = tidy_bvals(np.loadtxt(args.bval).ravel())
            np.savetxt(args.bval, bvals[None, :], fmt="%d")
            bvecs = np.loadtxt(args.bvec).reshape(3, -1)
            if bvecs.shape[1] != bvals.size:
                raise ValueError(f"{args.bvec} has {bvecs.shape[1]} directions, {args.bval} {bvals.size}")
            out = bids_bvecs(bvecs, nib.load(args.nifti).affine)
            if not np.array_equal(out, bvecs):
                np.savetxt(args.bvec, out, fmt="%.6f")
    except (OSError, ValueError) as err:
        print(f"cmig_finish_bids_dwi_series.py: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
