#!/usr/bin/env python3
"""Derive gradwarp settings from a BIDS sidecar JSON fields:
NonlinearGradientCorrection, NonlinearGradientCorrectionType, TablePosition,
Manufacturer, ManufacturersModelName, MagneticFieldStrength.

Values are emitted as KEY=VALUE, not as command-line flags, so any pipeline can
consume them under whatever flag names it uses.

Usage:
    cmig_get_bids_gradwarpinfo.py [--gw_dir=DIR] [--json] <sidecar.json>

    cmig_get_bids_gradwarpinfo.py sub-X_ses-Y_run-01_dwi.json          # KEY=VALUE lines
    cmig_get_bids_gradwarpinfo.py --json sub-X_ses-Y_run-01_dwi.json   # values as JSON

    eval "$(cmig_get_bids_gradwarpinfo.py "${up_base}.json")"           # set vars in bash

Options:
  sidecar
      path to the JSON sidecar
  --gw_dir GW_DIR
      folder holding the .mgh grad fields (default: the GW_DIR environment variable)
  --json
      emit a JSON object with the derivation intermediates

From Python, to get the values as a dict:
    import cmig_get_bids_gradwarpinfo as gw
    out = gw.resolve_sidecar("sub-X_ses-Y_run-01_dwi.json")
    out["gw_mode"], out["isocenter"], out["gw_field_L"]    # check out["gradwarp"] first

Emits
    gw_mode=2d/3d/1d  : the gradwarp mode to be applied to the data
                        1d: through-plane (z) gradwarp only
                        2d: in-plane gradwarp only
                        3d: full 3D gradwarp
    isocenter=0 or 1  : whether the image center is the isocenter of the magnet
                        0: the image center is not assumed to be the isocenter
                        1: the image center is assumed to be the isocenter
    gw_field_L=<path> : full path to the gradient nonlinearity field L, MGH format
    gw_field_P=<path> : full path to the gradient nonlinearity field P, MGH format
    gw_field_H=<path> : full path to the gradient nonlinearity field H, MGH format


Exits 0 with the settings on stdout when gradwarp applies; 0 with an empty stdout
and a NOTE on stderr when no gradwarp applies (Philips data, or the scanner
already applied 3D correction); 1 with an ERROR on stderr when the sidecar does
not determine a configuration.

Created:  10/05/2026 by Diliana Pecheva
Last mod: 11/02/2026 by Don Hagler

Renamed from gradwarp_args.py (abcd_preproc 81210b9)
"""

import argparse
import json
import os
import shlex
import sys
import warnings

from cmig_read_bids_sidecar import read_sidecar

# Every sidecar field resolve() reads. Absent fields are left out rather than
# set to None, so the `in` tests below still see them as absent.
SIDECAR_FIELDS = (
    "Manufacturer", "ManufacturersModelName", "MagneticFieldStrength",
    "NonlinearGradientCorrection", "NonlinearGradientCorrectionType",
    "TablePosition",
)

# folder of the gradient nonlinearity .mgh fields; set in the image, and by the job files
# for host-native runs
GW_DIR_DEFAULT = os.environ.get("GW_DIR", "")


class GradwarpError(Exception):
    """Raised when the sidecar does not determine a gradwarp configuration."""


# Reasons gradwarp does not apply. Kept distinct: "already corrected" means the
# data is fine, "unavailable" means it stays uncorrected.
NO_GRADWARP_NOTES = {
    "philips": "Philips data: no gradwarp correction to apply",
    "scanner_3d": "scanner applied 3D gradient correction; no gradwarp needed",
}


# Types 7, 8, 19 and 20 are not reachable from get_gwtype; the table is kept
# complete so it can be checked against the original by direct comparison.
# 19 and 20 have no files on disk, so reaching them raises on the missing file.
GRAD_STEM = {
    0: "siemens_sonata",
    1: "siemens_allegra",
    2: "ge_brm",
    3: "ge_crm",
    4: "siemens_avanto",
    5: "siemens_axxess",
    6: "siemens_quantum",
    7: "ge_whole",
    8: "ge_zoom",
    9: "ge_xrmb",
    10: "ge_xrmw",
    11: "siemens_skyra",
    12: "siemens_connectome",
    13: "siemens_prisma",
    14: "ge_HRMw",
    15: "ge_HRMB2",
    16: "siemens_verio",
    17: "siemens_vida",
    18: "siemens_CimaX",
    19: "siemens_AeraXJ",
    20: "siemens_AeraXQ",
}

# Exact-match models. Kept exact deliberately: a model string such as
# "Connectome_Skyra" matches neither "connectome" nor "skyra" and yields an
# unknown type rather than a guess.
SIEMENS_EXACT = {
    "sonata": 0, "trio": 0, "sonatavision": 0,
    "allegra": 1,
    "avanto": 4, "triotim": 4,
    "espree": 5, "axxess": 5,
    "symphony": 6, "symphonyvision": 6, "symphonytim": 6,
    "skyra": 11,
    "connectome": 12,
    "verio": 16, "biograph_mmr": 16,
}
SIEMENS_SUBSTR = {"prisma": 13, "vida": 17, "cimax": 18}

GE_EXACT = {
    "brm": 2,
    "crm": 3,
    "discovery mr450": 9, "discovery mr750": 9,
    "discovery mr750w": 10,
    "signa premier": 14,
    "signa uhp": 15,
}


def _require(info, key):
    if key not in info:
        raise GradwarpError(f"sidecar has no {key} field")
    return info[key]


# NonlinearGradientCorrectionType records what the SCANNER already did; the
# value here is what is left for us to correct. The inversion is deliberate.
# 1D/2D/3D are the MATLAB spellings; ND/DIS2D/DIS3D are the DICOM ImageType
# distortion tokens that dcm2niix propagates.
#   ND    = no distortion correction applied
#   DIS2D = 2D in-plane correction applied
#   DIS3D = 3D correction applied
# None means the scanner corrected everything and no gradwarp is needed.
UNWARPFLAG_BY_TYPE = {
    "ND":    0,      # nothing done      -> unwarp 3D
    "1D":    2,      # scanner did z     -> in-plane left
    "2D":    1,      # scanner did plane -> z left
    "DIS2D": 1,
    "3D":    None,   # fully corrected
    "DIS3D": None,
}


def get_unwarpflag(info):
    """0/1/2 = unwarp needed; None = no gradwarp to apply.

    Philips sidecars carry no NonlinearGradientCorrection field. There are no
    Philips grad fields to apply either, so this warns and returns None rather
    than raising. (resolve() short-circuits Philips before reaching here, so
    this path only matters when the function is called directly.)

    Any other sidecar missing the field raises: absence tells us nothing about
    what the scanner did, and assuming "nothing" would double-correct data that
    was already corrected.
    """
    if "NonlinearGradientCorrection" not in info:
        if "philips" in str(info.get("Manufacturer", "")).lower():
            warnings.warn(
                "Philips sidecar has no NonlinearGradientCorrection field; "
                "no gradwarp correction to apply", stacklevel=2)
            return None
        raise GradwarpError("sidecar has no NonlinearGradientCorrection field")
    if not info["NonlinearGradientCorrection"]:
        return 0
    raw = _require(info, "NonlinearGradientCorrectionType")
    kind = str(raw).strip().upper()
    # Membership test, not .get(): None is a meaningful value in the table, so
    # a default would conflate "fully corrected" with "unrecognised token".
    if kind not in UNWARPFLAG_BY_TYPE:
        raise GradwarpError(f"unrecognised NonlinearGradientCorrectionType {raw!r}")
    return UNWARPFLAG_BY_TYPE[kind]


def get_isoctrflag(info):
    """0 when scanning at isocenter, 1 otherwise. Passed to --isocenter as-is."""
    pos = _require(info, "TablePosition")
    if not isinstance(pos, (list, tuple)):
        pos = [pos]
    return 0 if all(v == 0 for v in pos) else 1


def get_gwtype(info):
    """Gradient type number, or None if the scanner is not recognised.

    The manufacturer tests are substring matches, so "ge" matches
    "GE MEDICAL SYSTEMS"; siemens is tested first, so the two cannot collide.
    """
    manu = str(_require(info, "Manufacturer")).lower()
    model = str(_require(info, "ManufacturersModelName")).lower().strip()

    if "siemens" in manu:
        if model in SIEMENS_EXACT:
            return SIEMENS_EXACT[model]
        for needle, gwtype in SIEMENS_SUBSTR.items():
            if needle in model:
                return gwtype
    elif "ge" in manu:
        if model in GE_EXACT:
            return GE_EXACT[model]
        if model == "signa creator" and _require(info, "MagneticFieldStrength") == 3:
            return 9
    return None


def gw_mode_for(unwarpflag):
    """0/1/2 -> the --gw_mode value.

    None has no mode: it means the scanner corrected everything and no gradwarp
    flag is emitted at all. resolve() returns before reaching this.
    """
    if unwarpflag not in (0, 1, 2):
        raise GradwarpError(
            f"no gw_mode for unwarpflag {unwarpflag!r}"
            + ("; None means no gradwarp applies" if unwarpflag is None else ""))
    return {0: "3d", 1: "1d", 2: "2d"}[unwarpflag]


def grad_files(gwtype, gw_dir):
    """(L, P, H) grad field paths for a gradient type, verified to exist."""
    if not gw_dir:
        raise GradwarpError("no gradwarp field folder: set GW_DIR or pass --gw_dir")
    if gwtype not in GRAD_STEM:
        raise GradwarpError(f"no grad files defined for gwtype {gwtype!r}")
    stem = GRAD_STEM[gwtype]
    paths = tuple(
        os.path.join(gw_dir, f"{stem}_d{axis}.mgh") for axis in "LPH"
    )
    for path in paths:
        if not os.path.isfile(path):
            raise GradwarpError(f"grad file not found: {path}")
    return paths


def resolve(info, gw_dir):
    """Sidecar dict -> {"gradwarp": False, "reason": ...} or the full settings.

    Manufacturer is checked before anything else because Philips sidecars do
    not carry NonlinearGradientCorrection at all, so requiring that field first
    would fail them. There are no Philips grad fields to apply in any case.

    unwarpflag is then resolved before gwtype, so a scan the scanner already
    corrected never needs its scanner recognised.
    """
    manu = str(_require(info, "Manufacturer")).lower()
    if "philips" in manu:
        return {"gradwarp": False, "reason": "philips"}

    unwarpflag = get_unwarpflag(info)
    if unwarpflag is None:
        return {"gradwarp": False, "reason": "scanner_3d"}

    isocenter = get_isoctrflag(info)

    gwtype = get_gwtype(info)
    if gwtype is None:
        raise GradwarpError(
            "unknown gradient type for "
            f"Manufacturer={info.get('Manufacturer')!r} "
            f"ManufacturersModelName={info.get('ManufacturersModelName')!r}")

    field_l, field_p, field_h = grad_files(gwtype, gw_dir)
    return {
        "gradwarp": True,
        "gwtype": gwtype,
        "unwarpflag": unwarpflag,
        "gw_mode": gw_mode_for(unwarpflag),
        "isocenter": isocenter,
        "gw_field_L": field_l,
        "gw_field_P": field_p,
        "gw_field_H": field_h,
    }


def resolve_sidecar(path, gw_dir=GW_DIR_DEFAULT):
    """Sidecar path -> settings dict, for use from Python rather than the CLI.

    Unlike main(), prints nothing and returns the values. Raises GradwarpError
    if the sidecar does not determine a configuration.
    """
    return resolve(read_sidecar(path, SIDECAR_FIELDS, required=False), gw_dir)


SETTING_KEYS = ("gw_mode", "isocenter", "gw_field_L", "gw_field_P", "gw_field_H")


def format_settings(settings):
    """KEY=VALUE lines, shell-quoted so `eval` is safe.

    Deliberately not flag names: each consumer builds its own. A pipeline that
    spells these differently reads the values and renames them itself.
    """
    return "\n".join(
        f"{key}={shlex.quote(str(settings[key]))}" for key in SETTING_KEYS)


def main(argv=None):
    # --help prints the module docstring, which describes every option
    ap = argparse.ArgumentParser(description=__doc__, usage=argparse.SUPPRESS, add_help=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    ap.add_argument("sidecar", help=argparse.SUPPRESS)
    ap.add_argument("--gw_dir", default=GW_DIR_DEFAULT,
                    help=argparse.SUPPRESS)
    ap.add_argument("--json", action="store_true",
                    help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    try:
        info = read_sidecar(args.sidecar, SIDECAR_FIELDS, required=False)
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"ERROR: cannot read {args.sidecar}: {exc}\n")
        return 1

    try:
        settings = resolve(info, args.gw_dir)
    except GradwarpError as exc:
        sys.stderr.write(f"ERROR: {args.sidecar}: {exc}\n")
        return 1

    if not settings["gradwarp"]:
        sys.stderr.write(
            f"NOTE: {args.sidecar}: {NO_GRADWARP_NOTES[settings['reason']]}\n")
        if args.json:
            print(json.dumps(settings, indent=2))
        return 0

    print(json.dumps(settings, indent=2) if args.json else format_settings(settings))
    return 0


if __name__ == "__main__":
    sys.exit(main())
