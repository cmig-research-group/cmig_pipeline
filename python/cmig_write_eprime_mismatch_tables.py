#!/usr/bin/env python3
"""
cmig_write_eprime_mismatch_tables.py

Write the REDCap eprime-scan delay exemptions (eprime_mismatch_ok_{task} in
the DAIC tracking table) as release-layout tables, so that
cmig_mkjobs_bids_pipeline.py reads them from the same layout as the release's
own tables.

Writes {outdir}/mr_y_qc__raw__tfmri__{task}.tsv for task mid, sst and nback,
with participant_id, session_id and
mr_y_qc__raw__tfmri__{task}__eprime__tdiff__ign_indicator, the release name
of the flag: one row per BIDS session in the tracking table, the value as the
table holds it (1, 0 or n/a), a session with several rows taking 1 when any
row sets it. Rows whose id_redcap or redcap_event_name does not map to BIDS
are skipped.

Usage
-----
  cmig_write_eprime_mismatch_tables.py --tracking_csv DAL_ABCD_QC_daic_tracking.csv \\
      --outdir /path/to/output [--verbose 0|1]

Created:  10/20/2026 by Don Hagler
Last mod: 10/20/2026 by Don Hagler
"""
import argparse
import csv
import os
import sys

from cmig_convert_abcd_eprime import bids_sub, bids_ses

SCRIPTNAME = 'cmig_write_eprime_mismatch_tables'
TASKS = ('mid', 'sst', 'nback')
TABLE = 'mr_y_qc__raw__tfmri__{0}'
COLUMN = 'mr_y_qc__raw__tfmri__{0}__eprime__tdiff__ign_indicator'


def parse_args():
    parser = argparse.ArgumentParser(
        description='Write the REDCap eprime_mismatch_ok flags as release-layout tables.')
    parser.add_argument('--tracking_csv', required=True,
                        help='DAIC tracking table (DAL_ABCD_QC_daic_tracking.csv)')
    parser.add_argument('--outdir', required=True,
                        help='directory for the mr_y_qc__raw__tfmri__{task}.tsv files')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1])
    return parser.parse_args()


def read_tracking(fname_csv):
    """{(participant_id, session_id): {task: '1', '0' or ''}} and the number of
    rows skipped for an unmapped id or event."""
    flags = {}
    nskip = 0
    with open(fname_csv) as fid:
        for row in csv.DictReader(fid):
            try:
                key = (bids_sub(row['id_redcap']), bids_ses(row['redcap_event_name']))
            except (KeyError, ValueError):
                nskip += 1
                continue
            entry = flags.setdefault(key, dict((t, '') for t in TASKS))
            for task in TASKS:
                val = (row.get('eprime_mismatch_ok_' + task) or '').strip()
                if val == '1' or (val == '0' and entry[task] == ''):
                    entry[task] = val
    return flags, nskip


def write_tables(flags, outdir, verbose=True):
    os.makedirs(outdir, exist_ok=True)
    for task in TASKS:
        fname = os.path.join(outdir, TABLE.format(task) + '.tsv')
        nset = 0
        with open(fname, 'w') as fid:
            writer = csv.writer(fid, delimiter='\t', lineterminator='\n')
            writer.writerow(['participant_id', 'session_id', COLUMN.format(task)])
            for (pid, sid) in sorted(flags):
                val = flags[(pid, sid)][task]
                nset += val == '1'
                writer.writerow([pid, sid, val if val else 'n/a'])
        if verbose:
            print('{0}: wrote {1}: {2} sessions, {3} with the flag set'.format(
                SCRIPTNAME, fname, len(flags), nset))


def run(args):
    if not os.path.isfile(args.tracking_csv):
        return 1, 'tracking csv not found: {0}'.format(args.tracking_csv)
    flags, nskip = read_tracking(args.tracking_csv)
    if args.verbose and nskip:
        print('{0}: {1} rows skipped, id or event not mapped to BIDS'.format(SCRIPTNAME, nskip))
    write_tables(flags, args.outdir, bool(args.verbose))
    return 0, ''


def main():
    errcode, errmsg = run(parse_args())
    if errmsg:
        print('{0}: ERROR: {1}'.format(SCRIPTNAME, errmsg), file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
