#!/usr/bin/env python3
"""
List the sessions in a BIDS rawdata tree, with the modality directories each
holds, into a TSV (sessions.tsv) that later steps read instead of rawdata.

One pass over rawdata, directory listings only: one listing per participant
and one per session. A session's _scans.tsv is written last by the rawdata
generation, so its mtime (scans_mtime, empty where there is none) marks a
session as ready.

An existing output file is merged with the new listing: rows found again are
updated, new rows added, and rows no longer found are kept with present 0.
With --sublist, only the listed participants' rows can be marked absent.
The file is written to <output_file>.new and then renamed, so an interrupted
run leaves the previous file as it was.

Columns:
  participant_id, session_id, anat, func, dwi, fmap, beh (0/1: directory
  present), scans_mtime, listed_date, present (0/1)

Created:  10/29/2026 by Don Hagler
Last mod: 10/29/2026 by Don Hagler
"""

import argparse
import csv
import datetime
import os
import sys

MODALITIES = ['anat', 'func', 'dwi', 'fmap', 'beh']
COLUMNS = (['participant_id', 'session_id'] + MODALITIES +
           ['scans_mtime', 'listed_date', 'present'])
TIME_FORMAT = '%Y-%m-%dT%H:%M:%S'


def _subdirs(path, prefix):
    """Sorted names of the subdirectories of path that start with prefix."""
    try:
        entries = list(os.scandir(path))
    except OSError:
        return []
    return sorted(e.name for e in entries
                  if e.name.startswith(prefix) and e.is_dir())


def list_sessions(raw_rootdir, sublist=None, seslist=None, verbose=False):
    """Return {(participant_id, session_id): row} for the sessions in rawdata."""
    listed_date = datetime.datetime.now().strftime(TIME_FORMAT)
    pids = sorted(sublist) if sublist else _subdirs(raw_rootdir, 'sub-')
    rows = {}
    for n, pid in enumerate(pids):
        if verbose and n % 1000 == 0:
            print('{0} of {1} participants'.format(n, len(pids)))
        for sid in _subdirs(os.path.join(raw_rootdir, pid), 'ses-'):
            if seslist and sid not in seslist:
                continue
            sesdir = os.path.join(raw_rootdir, pid, sid)
            try:
                names = set(os.listdir(sesdir))
            except OSError:
                continue
            row = {'participant_id': pid, 'session_id': sid}
            for mod in MODALITIES:
                row[mod] = int(mod in names)
            fname_scans = '{0}_{1}_scans.tsv'.format(pid, sid)
            row['scans_mtime'] = ''
            if fname_scans in names:
                mtime = os.path.getmtime(os.path.join(sesdir, fname_scans))
                row['scans_mtime'] = datetime.datetime.fromtimestamp(
                    mtime).strftime(TIME_FORMAT)
            row['listed_date'] = listed_date
            row['present'] = 1
            rows[(pid, sid)] = row
    return rows


def read_sessions(fname):
    """{(participant_id, session_id): row} from an existing sessions.tsv."""
    if not os.path.isfile(fname):
        return {}
    with open(fname) as fid:
        return {(r['participant_id'], r['session_id']): r
                for r in csv.DictReader(fid, delimiter='\t')}


def merge_sessions(old, new, sublist=None, seslist=None):
    """The new rows, plus the old rows not found again, marked present 0;
    an old row outside the participants or sessions listed is kept as is."""
    merged = dict(new)
    for key, row in old.items():
        if key in new:
            continue
        pid, sid = key
        row = dict(row)
        if (not sublist or pid in sublist) and (not seslist or sid in seslist):
            row['present'] = 0
        merged[key] = row
    return merged


def write_sessions(fname, rows):
    """Write the rows, sorted, to fname by way of fname.new."""
    fname_new = fname + '.new'
    with open(fname_new, 'w') as fid:
        writer = csv.DictWriter(fid, fieldnames=COLUMNS, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow(rows[key])
    os.replace(fname_new, fname)


def parse_args():
    parser = argparse.ArgumentParser(
        description='List the sessions in a BIDS rawdata tree, with their '
                    'modality directories, into sessions.tsv.')
    parser.add_argument('--raw_rootdir', required=True,
                        help='BIDS rawdata root directory')
    parser.add_argument('--output_file', required=True,
                        help='Output TSV (sessions.tsv), merged with an existing one')
    parser.add_argument('--sublist', nargs='+', default=None,
                        help='Participant IDs to list (default: all)')
    parser.add_argument('--seslist', nargs='+', default=None,
                        help='Session IDs to list (default: all)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1],
                        help='Display progress (default: 1)')
    return parser.parse_args()


def run(args):
    if not os.path.isdir(args.raw_rootdir):
        errmsg = 'raw_rootdir not found: {0}'.format(args.raw_rootdir)
        print('ERROR: {0}'.format(errmsg), file=sys.stderr)
        return 1, errmsg
    new = list_sessions(args.raw_rootdir, args.sublist, args.seslist,
                        bool(args.verbose))
    old = read_sessions(args.output_file)
    rows = merge_sessions(old, new, args.sublist, args.seslist)
    outdir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(outdir, exist_ok=True)
    write_sessions(args.output_file, rows)
    nready = sum(1 for r in new.values() if r['scans_mtime'])
    nabsent = sum(1 for r in rows.values() if str(r['present']) == '0')
    print('{0} sessions listed, {1} with _scans.tsv; {2} rows no longer found; '
          'wrote {3}'.format(len(new), nready, nabsent, args.output_file))
    return 0, None


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
