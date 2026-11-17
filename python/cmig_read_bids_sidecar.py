#!/usr/bin/env python3
"""Read fields from a BIDS JSON sidecar.

From Python:
    from cmig_read_bids_sidecar import read_sidecar, SidecarError
    vals = read_sidecar("sub-X_ses-Y_run-01_dwi.json", ["Manufacturer", "TablePosition"])

Usage: cmig_read_bids_sidecar.py SIDECAR.json FIELD [FIELD ...]

Arguments:
  SIDECAR.json
      the BIDS JSON sidecar to read
  FIELD
      a top-level field name; one or more

From bash, one Field=value line per field, shell-quoted so eval is safe:
    eval "$(cmig_read_bids_sidecar.py sub-X_ses-Y_run-01_dwi.json Manufacturer ManufacturersModelName)"

On the command line lists are printed space-separated and booleans as
true/false. Exits 1 with an ERROR on stderr if the file cannot be read or a
field is missing.

Created:  10/05/2026 by Diliana Pecheva
Last mod: 11/02/2026 by Don Hagler

Renamed from read_sidecar.py (abcd_preproc 81210b9)
"""

import json
import re
import shlex
import sys


class SidecarError(Exception):
    """Raised when a required field is missing from a sidecar."""


def read_sidecar(fname, fields, required=True):
    """Return {field: value} for the given top-level fields of a JSON sidecar.

    required=True raises SidecarError naming every missing field. With
    required=False missing fields are left out of the dict, so callers can
    still tell absent from present with `in`. File and JSON errors (OSError,
    ValueError) propagate unchanged.
    """
    with open(fname) as fh:
        info = json.load(fh)
    if not isinstance(info, dict):
        raise ValueError(f"{fname} is not a JSON object")
    missing = [f for f in fields if f not in info]
    if required and missing:
        raise SidecarError(f"{fname} has no {', '.join(missing)}")
    return {f: info[f] for f in fields if f in info}


def _shell_value(value):
    if isinstance(value, list):
        return " ".join(_shell_value(v) for v in value)
    if isinstance(value, bool) or value is None:
        return json.dumps(value)
    return str(value)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) < 2:
        sys.stderr.write(__doc__)
        return 1
    fname, fields = argv[0], argv[1:]
    bad = [f for f in fields if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", f)]
    if bad:
        sys.stderr.write(f"ERROR: not usable as shell variable names: {', '.join(bad)}\n")
        return 1
    try:
        values = read_sidecar(fname, fields)
    except (OSError, ValueError) as exc:
        sys.stderr.write(f"ERROR: cannot read {fname}: {exc}\n")
        return 1
    except SidecarError as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1
    for f in fields:
        print(f"{f}={shlex.quote(_shell_value(values[f]))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
