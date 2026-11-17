#!/usr/bin/env python3
"""FLIRT matrices for the atlas route of the b=0 -> T1w registration.

Usage:
  cmig_dwi_atlas_b0_start.py com B0 ATLAS_T2 ATLAS_MASK OUT.mat
      b=0 -> atlas T2 from the headers, shifted so the b=0's brain centre lands on the centre
      of ATLAS_MASK (voxels above 0.5)
  cmig_dwi_atlas_b0_start.py rigid MAT B0 OUT.mat
      MAT (a FLIRT matrix from B0) reduced to its rotation, with the translation that keeps the
      b=0's brain centre where MAT puts it

The b=0's brain centre is the mean position of its voxels above the 75th percentile (first
volume). FLIRT matrices map the input's FSL coordinates (voxels scaled to mm, x flipped when
the image is not radiological) to the reference's.

Created:  10/12/2026 by Diliana Pecheva
Last mod: 11/02/2026 by Don Hagler

Renamed from atlas_b0_start.py (abcd_preproc 81210b9)
"""

import sys

import nibabel as nib
import numpy as np


def fsl_scale(img):
    """voxel -> FSL coordinates."""
    s = np.diag(list(img.header.get_zooms()[:3]) + [1.0])
    if np.linalg.det(img.affine[:3, :3]) > 0:
        f = np.eye(4)
        f[0, 0] = -1
        f[0, 3] = img.shape[0] - 1
        s = s @ f
    return s


def b0_centre(b0):
    """The b=0's brain centre in voxel coordinates (homogeneous)."""
    d = np.asarray(b0.dataobj, np.float32)
    d = d[..., 0] if d.ndim == 4 else d
    return np.r_[np.argwhere(d > np.percentile(d, 75)).mean(0), 1]


def com_start(b0, atlas, mask):
    shift = np.eye(4)
    target = atlas.affine @ np.r_[np.argwhere(np.asarray(mask.dataobj) > 0.5).mean(0), 1]
    shift[:3, 3] = (target - b0.affine @ b0_centre(b0))[:3]
    return fsl_scale(atlas) @ np.linalg.inv(atlas.affine) @ shift @ b0.affine @ np.linalg.inv(fsl_scale(b0))


def rigid(mat, b0):
    c = fsl_scale(b0) @ b0_centre(b0)
    u, _, vt = np.linalg.svd(mat[:3, :3])
    r = np.eye(4)
    r[:3, :3] = u @ vt
    r[:3, 3] = (mat @ c)[:3] - r[:3, :3] @ c[:3]
    return r


def main(argv):
    if len(argv) == 6 and argv[1] == "com":
        m = com_start(nib.load(argv[2]), nib.load(argv[3]), nib.load(argv[4]))
    elif len(argv) == 5 and argv[1] == "rigid":
        m = rigid(np.loadtxt(argv[2]), nib.load(argv[3]))
    else:
        sys.exit(__doc__)
    np.savetxt(argv[-1], m, fmt="%.10f")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
