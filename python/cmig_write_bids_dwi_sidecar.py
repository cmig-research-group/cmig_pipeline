#!/usr/bin/env python3
"""Write the JSON sidecar for the processed DWI.

Usage: cmig_write_bids_dwi_sidecar.py --in_json RUN01.json --source RUN01.nii[.gz]
           [--source RUN02.nii[.gz]] --dwi OUT_dwi.nii.gz --gradwarp 0|1
           [--outlier_map PROC_outlier_map.txt] [--slice_resids PROC_slice_resids.txt]
           [--bval OUT_dwi.bval] [--spatial_reference URI] [--shim_match 0|1]
           [--fmap_selection label|shim|time] [--report processing_report.json]
           [--kind dwi|fmap] [--sources URI ...] [--intended_for URI ...] -o OUT_dwi.json

Options:
  --in_json JSON
      input sidecar of the source scan; its COPY_FIELDS are carried over
  --source NIFTI
      input image the output is made from (SourceVolumeSize, SourceAffine); repeat for
      several, one entry per image
  --dwi NIFTI
      the processed image the sidecar describes (VolumeSize, Affine, SliceThickness)
  --gradwarp 0|1
      1 if gradient nonlinearity correction was applied (GradUnwarpApplied)
  --outlier_map TXT
      TORTOISE outlier map (volumes x slices): adds outlier counts as
      MMIL_Process_DTI_v2.m's censor_info, per shell with --bval
  --slice_resids TXT
      TORTOISE slice residuals (volumes x slices): adds per-shell statistics; needs
      --bval
  --bval BVAL
      b-values of the processed image: adds VolumeCount, Shells, ShellCount and
      ShellVolumeCount
  --spatial_reference URI
      SpatialReference: the image whose space the output is in (BIDS URI)
  --shim_match 0|1
      ShimMatchesFieldMap: 1 if the field map's ShimSetting equals the scan's
  --fmap_selection label|shim|time
      FieldMapSelection: how cmig_select_bids_fmap.py chose the field map
  --report JSON
      TORTOISE's processing_report.json: adds TORTOISEVersion, ReversePEVolumeCount
      (down_nvols) and ReversePEShellCount (down_nshells); TORTOISEVersion is the image's
      versions.json tortoise entry, as in GeneratedBy and the status tortoise_version
      (cmig_sing_utils.read_versions), else the report's TORTOISE_ver without spaces
  --kind dwi|fmap
      what the processed image is, for its Description: a processed DWI series or the
      corrected reverse phase-encoding field map
  --sources URI [URI ...]
      Sources: the input images (BIDS URIs, e.g. bids:rawdata:sub-X/ses-Y/dwi/...)
  --intended_for URI [URI ...]
      IntendedFor, replacing the input sidecar's (whose names are the rawdata scans)
  -o, --output JSON
      sidecar to write

Built like cmig_filter_bids_info.m: COPY_FIELDS are copied from the input
sidecar when present, then the pipeline's fields are merged over them.
Affines are vox2ras to RAS+ world coordinates (the NIfTI convention) with
1-based voxel indices, as the MMPS streams write them. PhaseEncodingDirection is restated
for the output's voxel axes (the source axis taken to the world with SourceAffine, then
the nearest output axis). Per-shell fields are lists aligned with Shells (b rounded to 50,
as autoqc_dMRI_censor_info.m).

Created:  10/05/2026 by Diliana Pecheva
Prev mod: 10/07/2026 by Diliana Pecheva
Prev mod: 11/02/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/13/2026 by Don Hagler

Renamed from write_dwi_sidecar.py (abcd_preproc 81210b9)
"""

import argparse
import json

import nibabel as nib
import numpy as np

import cmig_sing_utils

# taglist of cmig_filter_bids_info.m (Last Mod 08/26/2026)
# Description by --kind
DESCRIPTIONS = {
    "dwi": "DWI processed with TORTOISE: denoising, Gibbs-ringing correction, motion, eddy-current "
           "and slice-to-volume correction with outlier replacement, susceptibility distortion "
           "correction (DRBUDDI) and, where applicable, gradient nonlinearity correction; "
           "resampled into the space of the sMRI image",
    "fmap": "Reverse phase-encoding b=0 field map, corrected by TORTOISE together with the DWI; "
            "on the DWI grid",
}
COPY_FIELDS = (
    "AcquisitionDateTime", "AcquisitionTime", "B0FieldIdentifier",
    "B0FieldSource", "BodyPart", "ConversionSoftware",
    "ConversionSoftwareVersion", "DeviceSerialNumber", "DwellTime",
    "EchoTime", "EffectiveEchoSpacing", "FlipAngle",
    "IntendedFor", "InversionTime", "MRAcquisitionType",
    "MagneticFieldStrength", "Manufacturer", "ManufacturersModelName",
    "Modality", "MultibandAccelerationFactor", "NonlinearGradientCorrection",
    "NonlinearGradientCorrectionType", "ParallelReductionFactorInPlane",
    "ParallelReductionFactorOutOfPlane",
    "PartialFourier", "PartialFourierDirection", "PatientID",
    "PatientName",
    "PatientPosition", "PhaseEncodingDirection", "ProtocolName",
    "PulseSequenceDetails", "ReceiveCoilActiveElements", "ReceiveCoilName",
    "RepetitionTime", "ScanOptions", "ScanningSequence",
    "SequenceName", "SequenceVariant", "SeriesDescription",
    "SeriesInstanceUID", "SeriesNumber", "SeriesTime",
    "SliceThickness", "SoftwareVersions", "SpacingBetweenSlices",
    "SpoilingState", "SpoilingType",
    "StudyDate", "StudyDescription", "StudyInstanceUID",
    "TaskName", "TotalReadoutTime",
)
# BIDS-named field <- vendor counterpart, used when the BIDS one is absent or empty
FILL_FIELDS = {"ReceiveCoilName": "CoilString", "SequenceName": "PulseSequenceName"}

# sidecar field <- processing_report.json key (--report); the report itself is not kept
REPORT_FIELDS = {"TORTOISEVersion": "TORTOISE_ver", "ReversePEVolumeCount": "down_nvols",
                 "ReversePEShellCount": "down_nshells"}



def tortoise_version(report):
    """TORTOISE's version as GeneratedBy and the status tortoise_version name it (the
    versions.json tortoise entry, the release tag, e.g. TORTOISEV4.3.0); without
    versions.json (the host), the report's TORTOISE_ver ("TORTOISE V4.3.0") without spaces."""
    version = cmig_sing_utils.read_versions().get("tortoise")
    if version:
        return version
    version = report.get("TORTOISE_ver")
    return "".join(str(version).split()) if version not in (None, "") else None


ONE_BASED = np.array([[1, 0, 0, -1], [0, 1, 0, -1], [0, 0, 1, -1], [0, 0, 0, 1]], float)


def mmps_affine(img):
    """vox2ras for 1-based voxel indices, as a nested list."""
    return (img.affine @ ONE_BASED).tolist()


def restate_direction(direction, src_affine, out_affine):
    """A BIDS axis direction (i, i-, j, ...) of the source image restated for the output
    image's voxel axes: the source axis in world coordinates, then the nearest output axis."""
    axis, sign = "ijk".index(direction[0]), -1.0 if direction.endswith("-") else 1.0
    world = sign * np.asarray(src_affine, float)[:3, axis]
    out = np.asarray(out_affine, float)[:3, :3]
    dots = (world / np.linalg.norm(world)) @ (out / np.linalg.norm(out, axis=0))
    k = int(np.argmax(np.abs(dots)))
    return "ijk"[k] + ("" if dots[k] > 0 else "-")


def volume_size(img):
    """[x, y, z, t]; t = 1 for a 3D image."""
    shape = list(img.shape[:4])
    return [int(n) for n in shape + [1] * (4 - len(shape))]


OUTLIER_FIELDS = ("OutlierSliceVolumeCount", "OutlierVolumeCount", "OutlierSliceCount",
                  "MaxOutlierVolumesPerSlice", "MaxOutlierSlicesPerVolume")


def shell_labels(bvals):
    """b-values rounded to the nearest 50, as autoqc_dMRI_censor_info.m groups shells."""
    return 50 * np.round(np.asarray(bvals, float) / 50)


def outlier_counts(o):
    """The five censor_info counts of a volumes x slices flag matrix, in OUTLIER_FIELDS order."""
    return (int(o.sum()), int(o.any(axis=1).sum()), int(o.any(axis=0).sum()),
            int(o.sum(axis=0).max()), int(o.sum(axis=1).max()))


def outlier_summary(outlier_map, shells=None):
    """Counts from TORTOISE's outlier map (volumes x slices, 1 = slice replaced),
    for the whole series and, given per-volume shell labels, per shell."""
    o = np.asarray(outlier_map) > 0
    info = dict(zip(OUTLIER_FIELDS, outlier_counts(o)))
    info["OutlierVolumes"] = [int(i) for i in np.flatnonzero(o.any(axis=1))]  # 0-based volume indices
    if shells is not None:
        per_shell = [outlier_counts(o[shells == b]) for b in np.unique(shells)]
        for k, field in enumerate(OUTLIER_FIELDS):
            info[field + "PerShell"] = [c[k] for c in per_shell]
    return info


def residual_summary(slice_resids, shells):
    """Per-shell statistics of TORTOISE's per-volume slice residual: the mean over
    the slices it scores (unscored slices are 0 in every volume)."""
    r = np.asarray(slice_resids, float)
    v = r[:, (r != 0).any(axis=0)].mean(axis=1)
    info = {f"SliceResidual{k}PerShell": [] for k in ("Mean", "Median", "Std", "P25", "P75")}
    for b in np.unique(shells):
        x = v[shells == b]
        p25, p75 = np.percentile(x, [25, 75])
        info["SliceResidualMeanPerShell"].append(float(x.mean()))
        info["SliceResidualMedianPerShell"].append(float(np.median(x)))
        info["SliceResidualStdPerShell"].append(float(x.std(ddof=1)) if len(x) > 1 else None)  # N-1, as MATLAB
        info["SliceResidualP25PerShell"].append(float(p25))
        info["SliceResidualP75PerShell"].append(float(p75))
    return info


def check_rows(name, a, nvol):
    if a.shape[0] != nvol:
        raise ValueError(f"{name} has {a.shape[0]} rows, the DWI {nvol} volumes")


def build_sidecar(in_info, sources, dwi, gradwarp, outlier_map=None, slice_resids=None,
                  bvals=None, spatial_reference=None, shim_match=None, fmap_selection=None,
                  kind=None, report=None, source_uris=None, intended_for=None):
    """Sidecar dict from the input sidecar, the source images and the output DWI."""
    info = {f: in_info[f] for f in COPY_FIELDS if f in in_info}
    for field, src in FILL_FIELDS.items():
        if info.get(field) in (None, "", []) and in_info.get(src) not in (None, "", []):
            info[field] = in_info[src]
    src_size = [volume_size(s) for s in sources]
    src_aff = [mmps_affine(s) for s in sources]
    one = len(sources) == 1
    info.update({
        "SliceThickness": float(dwi.header.get_zooms()[2]),  # output voxel size along k
        "RegistrationMatrixT1w": np.eye(4).tolist(),         # output is on the sMRI T1w grid
        "SourceAffine": src_aff[0] if one else src_aff,
        "Affine": mmps_affine(dwi),
        "SourceVolumeSize": src_size[0] if one else src_size,
        "VolumeSize": volume_size(dwi),
        "B0UnwarpApplied": True,
        "GradUnwarpApplied": gradwarp,
        "SkullStripped": False,
        "BiasCorrectionApplied": False,
    })
    if info.get("PhaseEncodingDirection"):
        info["PhaseEncodingDirection"] = restate_direction(
            info["PhaseEncodingDirection"], sources[0].affine, dwi.affine)
    if spatial_reference:
        info["SpatialReference"] = spatial_reference
    if shim_match is not None:
        info["ShimMatchesFieldMap"] = bool(shim_match)
    if fmap_selection:
        info["FieldMapSelection"] = fmap_selection
    if kind:
        info["Description"] = DESCRIPTIONS[kind]
    if source_uris:
        info["Sources"] = list(source_uris)
    if intended_for:
        info["IntendedFor"] = list(intended_for)
    if report is not None:
        info.update({f: report[k] for f, k in REPORT_FIELDS.items() if k in report})
        if "TORTOISEVersion" in info:
            info["TORTOISEVersion"] = tortoise_version(report)
    nvol = info["VolumeSize"][3]
    shells = None
    if bvals is not None:
        check_rows("bval", np.asarray(bvals).reshape(-1, 1), nvol)
        shells = shell_labels(bvals)
        info["VolumeCount"] = int(nvol)   # the number of frames
        info["Shells"] = [int(b) for b in np.unique(shells)]
        info["ShellCount"] = len(info["Shells"])
        info["ShellVolumeCount"] = [int((shells == b).sum()) for b in np.unique(shells)]
    if outlier_map is not None:
        check_rows("outlier map", outlier_map, nvol)
        info.update(outlier_summary(outlier_map, shells))
    if slice_resids is not None:
        if shells is None:
            raise ValueError("slice residuals are summarised per shell and need the b-values")
        check_rows("slice residuals", slice_resids, nvol)
        info.update(residual_summary(slice_resids, shells))
    return info


def main():
    # --help prints the module docstring, which describes every option
    ap = argparse.ArgumentParser(description=__doc__, usage=argparse.SUPPRESS, add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    ap.add_argument("--in_json", required=True, metavar="JSON",
                    help=argparse.SUPPRESS)
    ap.add_argument("--source", action="append", required=True, metavar="NIFTI",
                    help=argparse.SUPPRESS)
    ap.add_argument("--dwi", required=True, metavar="NIFTI",
                    help=argparse.SUPPRESS)
    ap.add_argument("--gradwarp", type=int, choices=(0, 1), required=True,
                    help=argparse.SUPPRESS)
    ap.add_argument("--outlier_map", metavar="TXT",
                    help=argparse.SUPPRESS)
    ap.add_argument("--slice_resids", metavar="TXT",
                    help=argparse.SUPPRESS)
    ap.add_argument("--bval", metavar="BVAL",
                    help=argparse.SUPPRESS)
    ap.add_argument("--spatial_reference", metavar="URI",
                    help=argparse.SUPPRESS)
    ap.add_argument("--shim_match", type=int, choices=(0, 1),
                    help=argparse.SUPPRESS)
    ap.add_argument("--fmap_selection", choices=("label", "shim", "time"),
                    help=argparse.SUPPRESS)
    ap.add_argument("--report", metavar="JSON", help=argparse.SUPPRESS)
    ap.add_argument("--kind", choices=tuple(DESCRIPTIONS), help=argparse.SUPPRESS)
    ap.add_argument("--sources", nargs="+", metavar="URI", help=argparse.SUPPRESS)
    ap.add_argument("--intended_for", nargs="+", metavar="URI", help=argparse.SUPPRESS)
    ap.add_argument("-o", "--output", required=True, metavar="JSON", help=argparse.SUPPRESS)
    args = ap.parse_args()
    with open(args.in_json) as fh:
        in_info = json.load(fh)
    load2d = lambda fname: np.loadtxt(fname, ndmin=2) if fname else None
    bvals = np.loadtxt(args.bval).ravel() if args.bval else None
    report = None
    if args.report:
        with open(args.report) as fh:
            report = json.load(fh)
    info = build_sidecar(in_info, [nib.load(s) for s in args.source],
                         nib.load(args.dwi), bool(args.gradwarp), load2d(args.outlier_map),
                         load2d(args.slice_resids), bvals, args.spatial_reference, args.shim_match,
                         args.fmap_selection, args.kind, report, args.sources,
                         args.intended_for)
    with open(args.output, "w") as fh:
        json.dump(info, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
