#!/usr/bin/env python3
"""
Write an edited copy of the DAIC merged manual raw QC csv
(DAL_ABCD_QC_merged_pcqcinfo.csv) whose consensus QC and revdisp carry the
values cmig_copy_manual_rawqc_bids_abcd.py computes from the per-reviewer
QC_XX columns (build_detail_rows, compute_consensus), which apply the
severe-score override through consensus_qc: a reviewer's QC 1 with a QU of
SEVERE_SCORE or more counts as 0.

Only rows of the SERIES_TYPE_MAP series types with a pGUID, and only where
the series was reviewed, are considered; only the QC and revdisp fields are
rewritten, keeping each field's surrounding whitespace. Every other byte is
copied unchanged.

Allowed edits:
  override     a reviewer QC 1 turned 0 by a severe QU changes the consensus
  ge_fieldmap  a GE field map (fMRI_FM, dMRI_FM: one series holding both
               directions) whose csv QC/revdisp differ from the consensus of
               its QC_XX columns: the merge of the two directions set each
               reviewer's QC_XX but kept the first direction's QC, nrev and
               revdisp
Any other difference stops the run without writing the csv.

Writes to --outdir: edits.tsv (one row per edited record), unexpected.tsv
(if any), override_sessions.tsv (with --id_events_file: the sessions of that
file with a QC edit; also written beside --outfile) and, with
--verify_flag 1, verify.txt.

Created:  11/11/2026 by Don Hagler
Last mod: 11/11/2026 by Don Hagler
"""

import argparse
import csv
import io
import os
import re
import sys
from collections import Counter

import cmig_copy_manual_rawqc_bids_abcd as cp

APPNAME = 'cmig_edit_manual_rawqc_abcd'

GE_FM = set(['fMRI_FM', 'dMRI_FM'])
ENC = dict(encoding='utf-8', errors='surrogateescape', newline='')

EDIT_COLUMNS = [
    'record', 'participant_id', 'session_id', 'VisitID', 'SeriesType',
    'SeriesInstanceUID', 'QC_csv', 'QC_new', 'revdisp_csv', 'revdisp_new',
    'n_reviewer_flipped', 'QC_qcxx', 'revdisp_qcxx',
]
EDIT_FLAGS = [
    'QC_edit', 'revdisp_edit', 'QC_override', 'revdisp_override',
    'QC_ge_fieldmap', 'revdisp_ge_fieldmap',
]
SESSION_COLUMNS = ['participant_id', 'session_id', 'n_series',
                   'series_types', 'reason']

csv.field_size_limit(sys.maxsize)

#############################################################################
# raw csv records
#############################################################################

def records(fp):
    # the raw text of each csv record (physical lines joined while a quoted
    # field is open), line endings included
    buf = []
    nq = 0
    for line in fp:
        buf.append(line)
        nq += line.count('"')
        if nq % 2 == 0:
            yield ''.join(buf)
            buf = []
            nq = 0
    if buf:
        yield ''.join(buf)


def split_raw(rec):
    # a raw record as its raw field strings (quotes kept) and the line
    # terminator; ','.join of the fields plus the terminator rebuilds it
    if rec.endswith('\r\n'):
        body, term = rec[:-2], '\r\n'
    elif rec.endswith('\n') or rec.endswith('\r'):
        body, term = rec[:-1], rec[-1]
    else:
        body, term = rec, ''
    fields = []
    cur = []
    inq = False
    for ch in body:
        if ch == '"':
            inq = not inq
            cur.append(ch)
        elif ch == ',' and not inq:
            fields.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    fields.append(''.join(cur))
    return fields, term


def parse(rec):
    rows = list(csv.reader(io.StringIO(rec, newline='')))
    if len(rows) != 1:
        raise RuntimeError('record parses to {0} rows'.format(len(rows)))
    return rows[0]


def set_value(raw, new):
    # the value inside a raw unquoted field replaced, its whitespace kept
    core = raw.strip()
    if '"' in raw or not core:
        raise RuntimeError('unexpected raw field {0!r}'.format(raw))
    i = raw.index(core)
    return raw[:i] + str(new) + raw[i + len(core):]


def fmt(x):
    return repr(x) if isinstance(x, str) else str(x)

#############################################################################
# consensus
#############################################################################

def row_consensus(keep, stype, revs):
    # the copy script's consensus with the override (summ), and without it
    # (raw: the same detail rows with their QU scores cleared, so
    # consensus_qc keeps every verdict)
    det = cp.build_detail_rows('p', 's', {'filename': 'f'}, keep, revs)
    summ = cp.compute_consensus('p', 's', ('f', stype), det)
    plain = [dict(r, **dict((c, None) for c in cp.QU_COLUMNS)) for r in det]
    raw = cp.compute_consensus('p', 's', ('f', stype), plain)
    flipped = sum(1 for r in det
                  if r['QC'] == 1 and cp.consensus_qc(r) == 0)
    return summ, raw, flipped

#############################################################################
# edit
#############################################################################

def edit_csv(args):
    # returns (edits, unexpected, nrec); writes args.outfile unless
    # unexpected differences are found or dry_flag
    types = set(cp.SERIES_TYPE_MAP.values())
    edits = []
    unexpected = []
    nrec = 0
    tmp = args.outfile + '.partial'
    out = None
    with open(args.infile, 'r', **ENC) as fp:
        it = records(fp)
        hrec = next(it)
        header = parse(hrec)
        revs = sorted(set(f[3:] for f in header if re.match(r'^QC_\w+$', f)))
        iq = header.index('QC')
        ir = header.index('revdisp')
        ist = header.index('SeriesType')
        ipg = header.index('pGUID')
        if not args.dry_flag:
            out = open(tmp, 'w', **ENC)
            out.write(hrec)
        for rec in it:
            nrec += 1
            row = parse(rec)
            new_rec = rec
            st = (row[ist] or '').strip()
            if st in types and row[ipg].strip():
                keep = dict((k, (v or '').strip())
                            for k, v in zip(header, row) if k)
                summ, raw, flipped = row_consensus(keep, st, revs)
                if summ['QC'] is not None:
                    qc_csv = cp.norm_int(row[iq])
                    rd_csv = cp.norm_int(row[ir])
                    dq = qc_csv != summ['QC']
                    dr = rd_csv != summ['revdisp']
                    if dq or dr:
                        rqc, rrd = raw['QC'], raw['revdisp']
                        # the cause of each field's edit: override = the
                        # override changes the QC_XX consensus; ge_fieldmap =
                        # the csv value differs from the QC_XX consensus
                        ovq = int(dq and flipped > 0 and rqc != summ['QC'])
                        ovr = int(dr and flipped > 0 and
                                  rrd != summ['revdisp'])
                        geq = int(dq and st in GE_FM and qc_csv != rqc)
                        ger = int(dr and st in GE_FM and rd_csv != rrd)
                        info = (nrec, cp.pguid_to_pid(keep['pGUID']),
                                cp.eventname_to_sid(keep['EventName']),
                                keep.get('VisitID'), st,
                                keep.get('SeriesInstanceUID'),
                                row[iq], summ['QC'], row[ir],
                                summ['revdisp'], flipped, rqc, rrd)
                        # every edited field needs a cause; a row that is not
                        # a GE field map must match the QC_XX consensus
                        ok = qc_csv is not None and rd_csv is not None and \
                            (not dq or ovq or geq) and \
                            (not dr or ovr or ger) and \
                            (st in GE_FM or (qc_csv == rqc and rd_csv == rrd))
                        if not ok:
                            unexpected.append(info)
                        else:
                            edits.append(info + (int(dq), int(dr), ovq, ovr,
                                                 geq, ger))
                            new_rec = rebuild(rec, row, iq, ir, dq, dr,
                                              summ, nrec)
            if out:
                out.write(new_rec)
    if out:
        out.close()
        if unexpected:
            os.remove(tmp)
        else:
            os.rename(tmp, args.outfile)
    return edits, unexpected, nrec


def rebuild(rec, row, iq, ir, dq, dr, summ, nrec):
    # the record with only the QC and revdisp fields replaced, checked by
    # parsing it back
    fields, term = split_raw(rec)
    if ','.join(fields) + term != rec or len(fields) != len(row):
        raise RuntimeError('raw split failed, record {0}'.format(nrec))
    if dq:
        fields[iq] = set_value(fields[iq], summ['QC'])
    if dr:
        fields[ir] = set_value(fields[ir], summ['revdisp'])
    new_rec = ','.join(fields) + term
    new_row = parse(new_rec)
    diff = [j for j in range(len(row)) if row[j] != new_row[j]]
    want = ([iq] if dq else []) + ([ir] if dr else [])
    if len(new_row) != len(row) or sorted(diff) != sorted(want):
        raise RuntimeError('rebuild check failed, record {0}'.format(nrec))
    return new_rec


def write_rows(fname, columns, rows):
    with open(fname, 'w') as f:
        f.write('\t'.join(columns) + '\n')
        for e in rows:
            f.write('\t'.join(fmt(x) for x in e) + '\n')


def report(edits, nrec):
    print('records', nrec)
    print('rows edited', len(edits))
    print('QC edits', sum(e[13] for e in edits),
          dict(Counter((e[6].strip(), e[7]) for e in edits if e[13])))
    print('revdisp edits', sum(e[14] for e in edits),
          '(in rows with QC unchanged:',
          sum(1 for e in edits if e[14] and not e[13]), ')')
    for lab, a, b in (('override', 15, 16), ('ge_fieldmap', 17, 18)):
        print('{0}: rows {1}, QC edits {2}, revdisp edits {3}'.format(
            lab, sum(1 for e in edits if e[a] or e[b]),
            sum(e[a] for e in edits), sum(e[b] for e in edits)))
    print('rows in both categories',
          sum(1 for e in edits if (e[15] or e[16]) and (e[17] or e[18])))
    print('fields with both causes',
          sum(1 for e in edits if e[15] and e[17]) +
          sum(1 for e in edits if e[16] and e[18]))


def override_sessions(edits, id_events_file):
    # the id_events sessions with a QC edit, and counts by reason
    ids = set((r['participant_id'], r['session_id'])
              for r in csv.DictReader(open(id_events_file)))
    ses = {}
    for e in edits:
        if not e[13]:
            continue
        s = ses.setdefault((e[1], e[2]), {'n': 0, 'types': [], 'why': set()})
        s['n'] += 1
        if e[4] not in s['types']:
            s['types'].append(e[4])
        if e[15]:
            s['why'].add('override')
        if e[17]:
            s['why'].add('ge_fieldmap')
    lines = []
    cnt = Counter()
    for key in sorted(ses):
        s = ses[key]
        why = 'both' if len(s['why']) == 2 else list(s['why'])[0]
        if key in ids:
            lines.append((key[0], key[1], s['n'], ','.join(s['types']), why))
            cnt[why] += 1
        cnt['all sessions (any) ' + why] += 1
    ser = Counter()
    for e in edits:
        if e[13]:
            why = 'both' if (e[15] and e[17]) else \
                ('override' if e[15] else 'ge_fieldmap')
            ser[(why, 'in id_events' if (e[1], e[2]) in ids
                 else 'not in id_events')] += 1
    print('QC-edited series by reason', sorted(ser.items()))
    print('id_events sessions by reason', sorted(
        (k, v) for k, v in cnt.items() if not k.startswith('all')))
    print('all sessions by reason', sorted(
        (k, v) for k, v in cnt.items() if k.startswith('all')))
    print('id_events sessions total', len(lines))
    return lines

#############################################################################
# verify
#############################################################################

def verify(args):
    # the edited csv against the original: same header and record count;
    # parsed fields differ only in QC and revdisp, only in the records of
    # edits.tsv, as planned; other records byte-identical; the copy script's
    # consensus on the edited csv matches its QC and revdisp. Returns the
    # number of failed checks
    plan = {}
    with open(os.path.join(args.outdir, 'edits.tsv')) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            plan[int(r['record'])] = r
    c = Counter()
    types = set(cp.SERIES_TYPE_MAP.values())
    with open(args.infile, 'r', **ENC) as fa, \
            open(args.outfile, 'r', **ENC) as fb:
        ia = records(fa)
        ib = records(fb)
        ha = next(ia)
        c['1 header byte-identical'] = int(ha == next(ib))
        header = parse(ha)
        iq = header.index('QC')
        ir = header.index('revdisp')
        revs = sorted(set(f[3:] for f in header if re.match(r'^QC_\w+$', f)))
        n = 0
        while True:
            a = next(ia, None)
            b = next(ib, None)
            if a is None and b is None:
                break
            n += 1
            if a is None or b is None:
                c['1 RECORD COUNT MISMATCH'] += 1
                break
            if a == b:
                c['3 records byte-identical'] += 1
                if n in plan:
                    c['3 PLANNED EDIT BUT IDENTICAL'] += 1
            else:
                c['3 records differing'] += 1
                if n not in plan:
                    c['3 DIFFERING BUT NOT PLANNED'] += 1
            ra = parse(a)
            rb = parse(b)
            if len(ra) != len(rb):
                c['2 FIELD COUNT DIFFERS'] += 1
            diff = [j for j in range(max(len(ra), len(rb)))
                    if j >= len(ra) or j >= len(rb) or ra[j] != rb[j]]
            for j in diff:
                c['2 field differs: ' + header[j]] += 1
                if j not in (iq, ir):
                    c['2 OTHER FIELD DIFFERS'] += 1
            if n in plan:
                p = plan[n]
                want = []
                if p['QC_edit'] == '1':
                    want.append(iq)
                    if rb[iq].strip() != p['QC_new']:
                        c['2 QC VALUE NOT AS PLANNED'] += 1
                    c['2 QC edit {0} -> {1}'.format(
                        ra[iq].strip(), rb[iq].strip())] += 1
                if p['revdisp_edit'] == '1':
                    want.append(ir)
                    if rb[ir].strip() != p['revdisp_new']:
                        c['2 REVDISP VALUE NOT AS PLANNED'] += 1
                    c['2 revdisp edit {0} -> {1}'.format(
                        ra[ir].strip(), rb[ir].strip())] += 1
                if sorted(diff) != sorted(want):
                    c['2 DIFF FIELDS NOT AS PLANNED'] += 1
            elif diff:
                c['2 DIFFERS IN UNPLANNED RECORD'] += 1
            # the copy script's consensus on the edited row
            row = dict(zip(header, rb))
            st = (row.get('SeriesType') or '').strip()
            if st not in types or not (row.get('pGUID') or '').strip():
                continue
            keep = dict((k, (v or '').strip()) for k, v in row.items() if k)
            c['4 rows'] += 1
            summ = row_consensus(keep, st, revs)[0]
            qc_csv = cp.norm_int(keep.get('QC'))
            if summ['QC'] is None:
                c['4 not reviewed'] += 1
            elif qc_csv != summ['QC']:
                c['4 QC DIFFERS'] += 1
            else:
                c['4 QC same'] += 1
            if summ['QC'] is not None and \
                    cp.norm_int(keep.get('revdisp')) != summ['revdisp']:
                c['4 REVDISP DIFFERS'] += 1
    c['1 records'] = n
    c['2 planned edits (edits.tsv rows)'] = len(plan)
    # the upper-case keys count failures
    nfail = sum(v for k, v in c.items() if k[2:].split(':')[0].isupper())
    nfail += 1 - c['1 header byte-identical']
    lines = ['{0:60s} {1}'.format(k, c[k]) for k in sorted(c)]
    lines.append('failed checks: {0}'.format(nfail))
    with open(os.path.join(args.outdir, 'verify.txt'), 'w') as f:
        f.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    return nfail

#############################################################################
# main
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Write an edited copy of the merged manual raw QC csv, '
                    'its consensus QC and revdisp recomputed with the '
                    'severe-score override.')
    parser.add_argument('--infile', required=True,
                        help='unedited DAL_ABCD_QC_merged_pcqcinfo.csv')
    parser.add_argument('--outfile', required=True,
                        help='edited csv to write')
    parser.add_argument('--outdir', required=True,
                        help='directory for edits.tsv, unexpected.tsv, '
                             'override_sessions.tsv and verify.txt')
    parser.add_argument('--id_events_file', default=None,
                        help='csv with participant_id, session_id: write '
                             'override_sessions.tsv for its sessions '
                             '(default: none)')
    parser.add_argument('--verify_flag', type=int, default=0, choices=[0, 1],
                        help='check the edited csv against the original '
                             '(default: 0)')
    parser.add_argument('--dry_flag', type=int, default=0, choices=[0, 1],
                        help='write only the edit tables, not the csv '
                             '(default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='overwrite an existing outfile (default: 0)')
    return parser.parse_args()


def run(args):
    if not os.path.isfile(args.infile):
        print('{0}: ERROR: infile {1} not found'.format(APPNAME, args.infile))
        return 1
    if os.path.abspath(args.infile) == os.path.abspath(args.outfile):
        print('{0}: ERROR: outfile is the infile'.format(APPNAME))
        return 1
    if os.path.exists(args.outfile) and not args.forceflag and \
            not args.dry_flag:
        print('{0}: ERROR: outfile {1} exists (use --forceflag 1)'.format(
            APPNAME, args.outfile))
        return 1
    os.makedirs(args.outdir, exist_ok=True)
    if not args.dry_flag:
        os.makedirs(os.path.dirname(os.path.abspath(args.outfile)),
                    exist_ok=True)
    edits, unexpected, nrec = edit_csv(args)
    write_rows(os.path.join(args.outdir, 'edits.tsv'),
               EDIT_COLUMNS + EDIT_FLAGS, edits)
    report(edits, nrec)
    print('unexpected', len(unexpected))
    if unexpected:
        write_rows(os.path.join(args.outdir, 'unexpected.tsv'),
                   EDIT_COLUMNS, unexpected)
        print('{0}: ERROR: unexpected differences, {1} not written'.format(
            APPNAME, args.outfile))
        return 1
    if args.id_events_file:
        lines = override_sessions(edits, args.id_events_file)
        dirs = [args.outdir]
        if not args.dry_flag:
            dirs.append(os.path.dirname(os.path.abspath(args.outfile)))
        for d in sorted(set(os.path.abspath(x) for x in dirs)):
            with open(os.path.join(d, 'override_sessions.tsv'), 'w') as f:
                f.write('\t'.join(SESSION_COLUMNS) + '\n')
                for x in lines:
                    f.write('\t'.join(str(v) for v in x) + '\n')
    if args.dry_flag:
        return 0
    print('{0}: wrote {1}'.format(APPNAME, args.outfile))
    if args.verify_flag and verify(args):
        print('{0}: ERROR: verify failed'.format(APPNAME))
        return 1
    return 0


def main():
    sys.exit(run(parse_args()))


if __name__ == '__main__':
    main()
