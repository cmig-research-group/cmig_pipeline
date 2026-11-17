#!/usr/bin/env python3
"""
Device bookkeeping for running the cmig pipeline over many sessions spread
across several storage devices: devices.tsv (one row per device) and
assignments.tsv (one row per participant, naming its device).

A participant is assigned, with all its sessions, to the enabled device with
the lowest projected fraction of its size after the assignment, among those
that stay under max_frac (projected) and hard_frac (df use, plus the space
held by running jobs). Projected use is base_used_gb (df use when the device
was added) plus nsess_assigned x avg_session_gb.

Sizes are in GiB (1024**3 bytes). Tables are written to <file>.new and then
renamed, so an interrupted write leaves the previous file as it was.

devices.tsv columns:
  device, rootdir, enabled (0/1), size_gb, base_used_gb, nsess_assigned,
  npart_assigned, projected_gb, projected_frac, df_used_gb, df_frac, df_date

assignments.tsv columns:
  participant_id, device, assigned_date, nsess_at_assign

Without a mode flag the CLI prints the device table; --add_device adds a row,
--update_df measures df of each device, and --dry_assign simulates assigning
the eligible participants of --sessions, writing nothing.

A session is eligible when present, with each required modality, in seslist,
with _scans.tsv (ready_only), and in the ready file where one is given (a
participant_id,session_id csv or tsv, e.g. a rawdata readiness check's list).

Created:  10/29/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Last mod: 11/16/2026 by Don Hagler
"""

import argparse
import csv
import datetime
import os
import sys

GIB = 1024 ** 3
TIME_FORMAT = '%Y-%m-%dT%H:%M:%S'
DEVICE_COLUMNS = ['device', 'rootdir', 'enabled', 'size_gb', 'base_used_gb',
                  'nsess_assigned', 'npart_assigned', 'projected_gb',
                  'projected_frac', 'df_used_gb', 'df_frac', 'df_date']
ASSIGN_COLUMNS = ['participant_id', 'device', 'assigned_date',
                  'nsess_at_assign']
INT_FIELDS = ['enabled', 'nsess_assigned', 'npart_assigned']
GB_FIELDS = ['size_gb', 'base_used_gb', 'projected_gb', 'df_used_gb']
FRAC_FIELDS = ['projected_frac', 'df_frac']


def _now():
    return datetime.datetime.now().strftime(TIME_FORMAT)


def _fmt(val, ndec=None):
    """A value as written to a TSV: '' for None, ndec decimals for a float."""
    if val is None:
        return ''
    if ndec is not None and isinstance(val, float):
        return '{0:.{1}f}'.format(val, ndec)
    return str(val)


def read_tsv(fname):
    """List of dicts from a TSV with a header row; [] if the file is missing."""
    if not os.path.isfile(fname):
        return []
    with open(fname) as fid:
        return [dict(r) for r in csv.DictReader(fid, delimiter='\t')]


def write_tsv(fname, rows, columns):
    """Write rows (dicts) to fname by way of fname.new; None is written as ''."""
    outdir = os.path.dirname(os.path.abspath(fname))
    os.makedirs(outdir, exist_ok=True)
    fname_new = fname + '.new'
    with open(fname_new, 'w') as fid:
        writer = csv.DictWriter(fid, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _fmt(row.get(k)) for k in columns})
    os.replace(fname_new, fname)


def _num(val, fun):
    if val is None or str(val).strip() == '':
        return None
    return fun(float(val)) if fun is int else fun(val)


def read_devices(fname):
    """devices.tsv rows with the numeric fields converted ('' -> None)."""
    devices = read_tsv(fname)
    for dev in devices:
        for key in INT_FIELDS:
            dev[key] = _num(dev.get(key), int)
        for key in GB_FIELDS + FRAC_FIELDS:
            dev[key] = _num(dev.get(key), float)
    return devices


def write_devices(fname, devices):
    """Write devices.tsv: GiB with 1 decimal, fractions with 3."""
    rows = []
    for dev in devices:
        row = dict(dev)
        for key in GB_FIELDS:
            row[key] = _fmt(dev.get(key), 1)
        for key in FRAC_FIELDS:
            row[key] = _fmt(dev.get(key), 3)
        rows.append(row)
    write_tsv(fname, rows, DEVICE_COLUMNS)


def read_assignments(fname):
    """{participant_id: row} from assignments.tsv."""
    return {r['participant_id']: r for r in read_tsv(fname)}


def write_assignments(fname, assignments):
    """Write assignments.tsv, sorted by participant_id."""
    write_tsv(fname, [assignments[k] for k in sorted(assignments)],
              ASSIGN_COLUMNS)


def _existing_path(path):
    """path, or its nearest existing parent."""
    path = os.path.abspath(path)
    while not os.path.exists(path):
        parent = os.path.dirname(path)
        if parent == path:
            break
        path = parent
    return path


def measure_df(path):
    """(size_gb, used_gb) of the file system holding path, or None."""
    try:
        st = os.statvfs(path)
    except OSError:
        return None
    size_gb = st.f_blocks * st.f_frsize / GIB
    used_gb = (st.f_blocks - st.f_bfree) * st.f_frsize / GIB
    return size_gb, used_gb


def update_df(devices, date=None):
    """Set df_used_gb, df_frac, df_date of each device; returns the names of
    the devices whose df failed (their values are left unchanged)."""
    date = date or _now()
    failed = []
    for dev in devices:
        df = measure_df(_existing_path(dev['rootdir']))
        if df is None or not df[0]:
            failed.append(dev['device'])
            continue
        size_gb, used_gb = df
        dev['df_used_gb'] = used_gb
        # fraction of size_gb in the table, as projected_frac uses
        dev['df_frac'] = used_gb / (dev.get('size_gb') or size_gb)
        dev['df_date'] = date
    return failed


def _update_one(dev, avg_session_gb):
    if dev.get('size_gb') is None or dev.get('base_used_gb') is None:
        return
    dev['projected_gb'] = (dev['base_used_gb'] +
                           (dev.get('nsess_assigned') or 0) * avg_session_gb)
    dev['projected_frac'] = dev['projected_gb'] / dev['size_gb']


def update_projected(devices, avg_session_gb):
    """projected_gb = base_used_gb + nsess_assigned x avg_session_gb;
    projected_frac = projected_gb / size_gb."""
    for dev in devices:
        _update_one(dev, avg_session_gb)


def choose_device(devices, nsess, avg_session_gb, max_frac, hard_frac,
                  running=None, transient_gb=0.0):
    """Name of the device for a participant with nsess sessions, or None.

    running is {device: number of running jobs}, each holding transient_gb.
    """
    running = running or {}
    best = None
    for dev in devices:
        if not dev.get('enabled') or not dev.get('size_gb'):
            continue
        size_gb = dev['size_gb']
        projected_gb = dev.get('projected_gb')
        if projected_gb is None:
            projected_gb = ((dev.get('base_used_gb') or 0.0) +
                            (dev.get('nsess_assigned') or 0) * avg_session_gb)
        after_frac = (projected_gb + nsess * avg_session_gb) / size_gb
        if after_frac >= max_frac:
            continue
        if dev.get('df_frac') is not None:
            df_gb = (dev['df_used_gb'] +
                     running.get(dev['device'], 0) * transient_gb)
            if df_gb / size_gb >= hard_frac:
                continue
        key = (after_frac, projected_gb, dev['device'])
        if best is None or key < best:
            best = key
    return best[2] if best else None


def assign(devices, assignments, participant_id, nsess, avg_session_gb,
           max_frac, hard_frac, running=None, transient_gb=0.0, date=None):
    """Device of participant_id, assigning one (and updating devices and
    assignments in memory) if it has none yet; None if no device fits."""
    if participant_id in assignments:
        return assignments[participant_id]['device']
    name = choose_device(devices, nsess, avg_session_gb, max_frac, hard_frac,
                         running, transient_gb)
    if name is None:
        return None
    assignments[participant_id] = {
        'participant_id': participant_id, 'device': name,
        'assigned_date': date or _now(), 'nsess_at_assign': nsess}
    for dev in devices:
        if dev['device'] == name:
            dev['nsess_assigned'] = (dev.get('nsess_assigned') or 0) + nsess
            dev['npart_assigned'] = (dev.get('npart_assigned') or 0) + 1
            _update_one(dev, avg_session_gb)
            break
    return name


def add_device(devices, name, rootdir, enabled=1):
    """Append a device row with size_gb and base_used_gb from df;
    returns (errcode, errmsg)."""
    if any(d['device'] == name for d in devices):
        return 1, 'device {0} already in the table'.format(name)
    df = measure_df(_existing_path(rootdir))
    if df is None:
        return 1, 'df failed for {0}'.format(rootdir)
    size_gb, used_gb = df
    dev = {'device': name, 'rootdir': rootdir, 'enabled': int(enabled),
           'size_gb': size_gb, 'base_used_gb': used_gb,
           'nsess_assigned': 0, 'npart_assigned': 0,
           'projected_gb': used_gb, 'projected_frac': used_gb / size_gb,
           'df_used_gb': used_gb, 'df_frac': used_gb / size_gb,
           'df_date': _now()}
    devices.append(dev)
    return 0, None


def read_ready_file(fname):
    """Set of (participant_id, session_id) from a csv or tsv with those
    columns (tab-delimited if the header has a tab)."""
    with open(fname) as fid:
        header = fid.readline()
        fid.seek(0)
        delim = '\t' if '\t' in header else ','
        rows = list(csv.DictReader(fid, delimiter=delim))
    if rows and not {'participant_id', 'session_id'} <= set(rows[0]):
        raise ValueError('{0}: no participant_id and session_id '
                         'columns'.format(fname))
    return {(r['participant_id'].strip(), r['session_id'].strip())
            for r in rows}


def eligible_sessions(sessions_rows, seslist=None, required=('anat',),
                      ready_only=True, ready_pairs=None):
    """{participant_id: [session_id, ...]} of the sessions present, with each
    required modality, in seslist (if given), with _scans.tsv when
    ready_only, and in ready_pairs (if given, a set of (participant_id,
    session_id))."""
    out = {}
    for row in sessions_rows:
        if str(row.get('present', '')).strip() != '1':
            continue
        if any(str(row.get(m, '')).strip() != '1' for m in required):
            continue
        if seslist and row['session_id'] not in seslist:
            continue
        if ready_only and not row.get('scans_mtime'):
            continue
        if ready_pairs is not None and \
                (row['participant_id'], row['session_id']) not in ready_pairs:
            continue
        out.setdefault(row['participant_id'], []).append(row['session_id'])
    for pid in out:
        out[pid].sort()
    return out


def print_devices(devices):
    hdr = ['device', 'enabled', 'size_gb', 'base_used_gb', 'npart', 'nsess',
           'projected_gb', 'proj_frac', 'df_used_gb', 'df_frac', 'df_date']
    print('\t'.join(hdr))
    for d in devices:
        print('\t'.join([d['device'], _fmt(d.get('enabled')),
                         _fmt(d.get('size_gb'), 1),
                         _fmt(d.get('base_used_gb'), 1),
                         _fmt(d.get('npart_assigned')),
                         _fmt(d.get('nsess_assigned')),
                         _fmt(d.get('projected_gb'), 1),
                         _fmt(d.get('projected_frac'), 3),
                         _fmt(d.get('df_used_gb'), 1),
                         _fmt(d.get('df_frac'), 3),
                         _fmt(d.get('df_date'))]))


def dry_assign(devices, assignments, eligible, avg_session_gb, max_frac,
               hard_frac):
    """Simulate assigning the unassigned participants of eligible in sorted
    order (devices and assignments are changed in memory); prints the result."""
    before = {d['device']: dict(d) for d in devices}
    nassigned, nnone, nsess_none, nprev = 0, 0, 0, 0
    for pid in sorted(eligible):
        if pid in assignments:
            nprev += 1
            continue
        nsess = len(eligible[pid])
        if assign(devices, assignments, pid, nsess, avg_session_gb, max_frac,
                  hard_frac) is None:
            nnone += 1
            nsess_none += nsess
        else:
            nassigned += 1
    print('{0} eligible participants, {1} already assigned, {2} assigned, '
          '{3} found no device ({4} sessions)'.format(
              len(eligible), nprev, nassigned, nnone, nsess_none))
    print('\t'.join(['device', 'npart', 'nsess', 'projected_gb', 'proj_frac',
                     'npart_after', 'nsess_after', 'projected_gb_after',
                     'proj_frac_after']))
    for d in devices:
        b = before[d['device']]
        print('\t'.join([d['device'], _fmt(b.get('npart_assigned')),
                         _fmt(b.get('nsess_assigned')),
                         _fmt(b.get('projected_gb'), 1),
                         _fmt(b.get('projected_frac'), 3),
                         _fmt(d.get('npart_assigned')),
                         _fmt(d.get('nsess_assigned')),
                         _fmt(d.get('projected_gb'), 1),
                         _fmt(d.get('projected_frac'), 3)]))


def parse_args():
    parser = argparse.ArgumentParser(
        description='Device bookkeeping: devices.tsv, df of each device, '
                    'and dry assignment of participants to devices.')
    parser.add_argument('--devices', required=True,
                        help='Device table (devices.tsv)')
    parser.add_argument('--add_device', nargs=2, metavar=('NAME', 'ROOTDIR'),
                        default=None,
                        help='Add a device, with size and use from df')
    parser.add_argument('--update_df', action='store_true',
                        help='Measure df of every device and write the table')
    parser.add_argument('--dry_assign', action='store_true',
                        help='Simulate assigning the eligible participants of '
                             '--sessions not yet in --assignments; writes nothing')
    parser.add_argument('--sessions', default=None,
                        help='Session table (sessions.tsv) for --dry_assign')
    parser.add_argument('--assignments', default=None,
                        help='Assignment table (assignments.tsv)')
    parser.add_argument('--avg_session_gb', type=float, default=12.0,
                        help='Projected GiB per session (default: 12)')
    parser.add_argument('--max_frac', type=float, default=0.70,
                        help='Projected fraction a device may not reach '
                             '(default: 0.70)')
    parser.add_argument('--hard_frac', type=float, default=0.90,
                        help='df fraction above which a device takes no new '
                             'participants (default: 0.90)')
    parser.add_argument('--seslist', nargs='+', default=None,
                        help='Session IDs eligible (default: all)')
    parser.add_argument('--ready_only', type=int, default=1, choices=[0, 1],
                        help='Only sessions with _scans.tsv (default: 1)')
    parser.add_argument('--ready_file', default=None,
                        help='participant_id,session_id csv or tsv: only '
                             'these sessions are eligible (default: all)')
    return parser.parse_args()


def _error(errmsg):
    print('ERROR: {0}'.format(errmsg), file=sys.stderr)
    return 1, errmsg


def run(args):
    devices = read_devices(args.devices)
    if args.add_device:
        name, rootdir = args.add_device
        errcode, errmsg = add_device(devices, name, rootdir)
        if errcode:
            return _error(errmsg)
        write_devices(args.devices, devices)
        print('added {0} to {1}'.format(name, args.devices))
        print_devices(devices)
    elif args.update_df:
        failed = update_df(devices)
        write_devices(args.devices, devices)
        print_devices(devices)
        if failed:
            return _error('df failed for {0}'.format(' '.join(failed)))
    elif args.dry_assign:
        if not args.sessions or not os.path.isfile(args.sessions):
            return _error('sessions file not found: {0}'.format(args.sessions))
        if not devices:
            return _error('no devices in {0}'.format(args.devices))
        assignments = (read_assignments(args.assignments)
                       if args.assignments else {})
        ready_pairs = None
        if args.ready_file:
            if not os.path.isfile(args.ready_file):
                return _error('ready file not found: {0}'.format(
                    args.ready_file))
            ready_pairs = read_ready_file(args.ready_file)
        eligible = eligible_sessions(read_tsv(args.sessions), args.seslist,
                                     ready_only=bool(args.ready_only),
                                     ready_pairs=ready_pairs)
        update_projected(devices, args.avg_session_gb)
        dry_assign(devices, assignments, eligible, args.avg_session_gb,
                   args.max_frac, args.hard_frac)
    else:
        if not devices:
            return _error('no devices in {0}'.format(args.devices))
        print_devices(devices)
    return 0, None


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
