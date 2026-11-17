#!/usr/bin/env python3
"""
Write edited copies of the ABCD manual post QC tables
(mr_y_qc__post__man__<modality>.tsv of release 7.0, and the
_desc-legacy.tsv tables of MMPS abcd_write_legacy_postqc.py), whose
consensus _score and __rev__disp_indicator are recomputed from the
per-reviewer QC_XX and QU_XX_<type> columns of the DAL_ABCD
<type>_qcinfo_indiv.csv files with the severe-score override: a reviewer's
QC 1 with any QU of cmig_copy_manual_rawqc_bids_abcd.SEVERE_SCORE or more
counts as 0. The consensus is that of compute_consensus there (and of
MMIL_Summarize_postQC): the mean QC less 0.01, rounded, so a split fails;
revdisp 1 where the QCs differ.

A table row is edited only where the consensus recomputed from the csv
without the override reproduces the row's _score and __rev_count (the
csvs hold reviews made after release 7.0), and the override changes it.
Rows that are not reproduced are listed, not edited. Only the _score and
__rev__disp_indicator fields change; every other byte is copied.

Writes to --outdir the edited tables (same file names) and to --logdir, per
modality and table: edits_<table>.tsv (the edited rows) and
not_reproduced_<table>.tsv (reviewed rows the csv does not reproduce), and,
with --verify_flag 1, verify.txt.

Created:  11/11/2026 by Don Hagler
Last mod: 11/11/2026 by Don Hagler
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, OrderedDict

import cmig_copy_manual_rawqc_bids_abcd as cp

APPNAME = 'cmig_edit_manual_postqc_abcd'

PREFIX = 'mr_y_qc__post__man__'
# modality: (csv stem, {csv QU type: release score type})
MODALITIES = OrderedDict([
    ('fsurf', ('fsurf', {'Artifact': 'artifact', 'inhomogeneity': 'inhomog',
                         'motion': 'mot', 'PialOver': 'pialover',
                         'WMUnder': 'wmunder'})),
    ('t2', ('sMRI_T2w', {'Artifact': 'artifact', 'inhomogeneity': 'inhomog',
                         'motion': 'mot', 'regT1': 'regist'})),
    ('dmri', ('dMRI', {'cutoff': 'cutoff', 'B0warp': 'distort',
                       'imgQual': 'quality', 'regT1': 'regist',
                       'fseg': 'segmnt'})),
    ('fmri', ('fMRI', {'cutoff': 'cutoff', 'B0warp': 'distort',
                       'imgQual': 'quality', 'regT1': 'regist'})),
])
CSV_PATTERN = 'DAL_ABCD_{0}_qcinfo_indiv.csv'
TABLES = [('release7.0', 'prev_dir', '{0}.tsv'),
          ('legacy', 'legacy_dir', '{0}_desc-legacy.tsv')]
NA = 'n/a'

EDIT_COLUMNS = ['participant_id', 'session_id', 'score_old', 'score_new',
                'disp_old', 'disp_new', 'rev_count', 'n_reviewer_flipped',
                'reviewers_flipped', 'qu_types']
NOTREP_COLUMNS = ['participant_id', 'session_id', 'score', 'rev_count',
                  'disp', 'in_csv', 'csv_score', 'csv_rev_count', 'csv_disp',
                  'override_score', 'override_disp', 'override_changes']

#############################################################################
# csv
#############################################################################

def visit_to_ids(visit):
    # VisitID <site>_INV<id>_<event> -> (sub-<id>, ses-<NN>A), or None
    m = re.match(r'^[^_]+_INV([^_]+)_(.+)$', visit.strip())
    if not m:
        return None
    event = m.group(2)
    if event == 'baseline':
        return 'sub-' + m.group(1), 'ses-00A'
    m2 = re.match(r'^(\d+)year$', event)
    if not m2:
        return None
    return 'sub-' + m.group(1), 'ses-{0:02d}A'.format(int(m2.group(1)))


def reviewer_qc(qc, qus):
    # the reviewer's verdict with the override: QC 1 with a severe QU is 0
    if qc == 1 and any(q is not None and q >= cp.SEVERE_SCORE for q in qus):
        return 0
    return qc


def consensus(qcs):
    # (QC, nrev, revdisp) by compute_consensus of the copy script
    rows = [dict([('QC', q), ('notes', None)] +
                 [(c, None) for c in cp.QU_COLUMNS]) for q in qcs]
    summ = cp.compute_consensus('', '', ('', ''), rows)
    return summ['QC'], summ['nrev'], summ['revdisp']


def read_indiv(fname, types):
    # {(pid, sid): visit} with the csv consensus and the consensus without
    # and with the override, from the per-reviewer columns
    with open(fname, newline='') as f:
        rd = csv.DictReader(f)
        header = rd.fieldnames
        revs = [c[3:] for c in header if re.match(r'^QC_[A-Z]+$', c)]
        visits = {}
        for r in rd:
            ids = visit_to_ids(r['VisitID'])
            if ids is None or ids in visits:
                raise RuntimeError('{0}: bad or repeated VisitID {1}'.format(
                    fname, r['VisitID']))
            plain = []
            over = []
            flipped = []
            qu_types = set()
            for rev in revs:
                qc = cp.norm_int(r['QC_' + rev])
                if qc is None:
                    continue
                qus = dict((t, cp.norm_float(r['QU_{0}_{1}'.format(rev, t)]))
                           for t in types)
                eff = reviewer_qc(qc, list(qus.values()))
                plain.append(qc)
                over.append(eff)
                if eff != qc:
                    flipped.append(rev)
                    qu_types.update(types[t] for t, q in qus.items()
                                    if q is not None and q >= cp.SEVERE_SCORE)
            visits[ids] = {
                'csv': (cp.norm_int(r['QC']), cp.norm_int(r['nrev']),
                        cp.norm_int(r['revdisp'])),
                'plain': consensus(plain),
                'over': consensus(over),
                'flipped': flipped,
                'qu_types': sorted(qu_types),
            }
    return visits

#############################################################################
# tables
#############################################################################

def read_lines(fname):
    # the header fields and each line as (fields, terminator)
    with open(fname, newline='') as f:
        lines = f.readlines()
    out = []
    for line in lines:
        body = line.rstrip('\r\n')
        out.append((body.split('\t'), line[len(body):]))
    return out


def fmt_val(v):
    return NA if v is None else str(v)


def edit_table(fname, fname_out, mod, visits):
    # returns (edits, not_reproduced, counts); writes fname_out
    pre = PREFIX + mod
    lines = read_lines(fname)
    header = lines[0][0]
    nf = len(header)
    iq = header.index(pre + '_score')
    ir = header.index(pre + '__rev__disp_indicator')
    ic = header.index(pre + '__rev_count')
    edits = []
    notrep = []
    c = Counter()
    out = [lines[0]]
    for fields, term in lines[1:]:
        if len(fields) != nf:
            raise RuntimeError('{0}: {1} fields in a row, header {2}'.format(
                fname, len(fields), nf))
        c['rows'] += 1
        key = (fields[0], fields[1])
        v = visits.get(key)
        score, disp, cnt = fields[iq], fields[ir], fields[ic]
        reviewed = score != NA
        if reviewed:
            c['reviewed'] += 1
        if v is None:
            if reviewed:
                c['reviewed, not in csv'] += 1
                notrep.append((key[0], key[1], score, cnt, disp, 0) +
                              ('',) * 6)
            out.append((fields, term))
            continue
        changes = v['over'][0] != v['plain'][0] or \
            v['over'][2] != v['plain'][2]
        if not reviewed:
            if v['plain'][0] is not None:
                c['not reviewed, reviewed in csv'] += 1
            out.append((fields, term))
            continue
        same = (score, cnt) == (fmt_val(v['plain'][0]),
                                fmt_val(v['plain'][1]))
        if same and disp != fmt_val(v['plain'][2]):
            c['reviewed, reproduced, disp differs'] += 1
        if not same:
            c['reviewed, not reproduced'] += 1
            if changes:
                c['reviewed, not reproduced, override changes csv'] += 1
            notrep.append((key[0], key[1], score, cnt, disp, 1) +
                          tuple(fmt_val(x) for x in v['plain']) +
                          (fmt_val(v['over'][0]), fmt_val(v['over'][2]),
                           int(changes)))
            out.append((fields, term))
            continue
        c['reviewed, reproduced'] += 1
        if not changes:
            out.append((fields, term))
            continue
        new = list(fields)
        new[iq] = fmt_val(v['over'][0])
        new[ir] = fmt_val(v['over'][2])
        c['edited'] += 1
        if new[iq] != score:
            c['score {0} -> {1}'.format(score, new[iq])] += 1
        if new[ir] != disp:
            c['disp {0} -> {1}'.format(disp, new[ir])] += 1
        edits.append((key[0], key[1], score, new[iq], disp, new[ir], cnt,
                      len(v['flipped']), ','.join(v['flipped']),
                      ','.join(v['qu_types'])))
        out.append((new, term))
    tmp = fname_out + '.partial'
    with open(tmp, 'w', newline='') as f:
        for fields, term in out:
            f.write('\t'.join(fields) + term)
    os.rename(tmp, fname_out)
    return edits, notrep, c


def write_rows(fname, columns, rows):
    with open(fname, 'w') as f:
        f.write('\t'.join(columns) + '\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')

#############################################################################
# verify
#############################################################################

def verify_table(fname, fname_out, mod, edits, visits):
    # the edited table against the original: same header and line count;
    # unedited lines byte-identical; edited lines differ only in _score and
    # __rev__disp_indicator, as planned; every edited value the override
    # consensus of the csv. Returns a Counter (upper-case keys are failures)
    pre = PREFIX + mod
    c = Counter()
    with open(fname, 'rb') as f:
        a = f.read().split(b'\n')
    with open(fname_out, 'rb') as f:
        b = f.read().split(b'\n')
    if len(a) != len(b):
        c['LINE COUNT DIFFERS'] += 1
        return c
    c['header identical'] = int(a[0] == b[0])
    if a[0] != b[0]:
        c['HEADER DIFFERS'] += 1
    header = a[0].decode().split('\t')
    want = (header.index(pre + '_score'),
            header.index(pre + '__rev__disp_indicator'))
    plan = dict(((e[0], e[1]), e) for e in edits)
    for la, lb in zip(a[1:], b[1:]):
        fa = la.decode().split('\t')
        fb = lb.decode().split('\t')
        key = tuple(fa[:2]) if len(fa) > 1 else None
        if la == lb:
            c['lines identical'] += 1
            if key in plan:
                c['PLANNED BUT IDENTICAL'] += 1
            continue
        c['lines differing'] += 1
        if key not in plan:
            c['DIFFERS BUT NOT PLANNED'] += 1
            continue
        diff = [j for j in range(len(fa)) if fa[j] != fb[j]]
        if len(fa) != len(fb) or any(j not in want for j in diff):
            c['OTHER FIELD DIFFERS'] += 1
        e = plan[key]
        over = visits[key]['over']
        if (fb[want[0]], fb[want[1]]) != (e[3], e[5]) or \
                (fb[want[0]], fb[want[1]]) != (str(over[0]), str(over[2])):
            c['VALUE NOT AS PLANNED'] += 1
    c['planned edits'] = len(plan)
    return c

#############################################################################
# main
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Write edited copies of the release 7.0 and legacy '
                    'manual post QC tables, the consensus recomputed from '
                    'the qcinfo_indiv csvs with the severe-score override.')
    parser.add_argument('--indiv_dir', required=True,
                        help='directory with DAL_ABCD_<type>_qcinfo_indiv.csv')
    parser.add_argument('--prev_dir', required=True,
                        help='release 7.0 {0}<modality>.tsv tables'.format(
                            PREFIX))
    parser.add_argument('--legacy_dir', required=True,
                        help='{0}<modality>_desc-legacy.tsv tables'.format(
                            PREFIX))
    parser.add_argument('--outdir', required=True,
                        help='directory for the edited tables')
    parser.add_argument('--logdir', required=True,
                        help='directory for the edit tables and verify.txt')
    parser.add_argument('--modality', nargs='+', default=list(MODALITIES),
                        choices=list(MODALITIES),
                        help='modalities (default: all)')
    parser.add_argument('--verify_flag', type=int, default=0, choices=[0, 1],
                        help='check the edited tables against the originals '
                             '(default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='overwrite existing edited tables (default: 0)')
    return parser.parse_args()


def run(args):
    for d in (args.indiv_dir, args.prev_dir, args.legacy_dir):
        if not os.path.isdir(d):
            print('{0}: ERROR: {1} not found'.format(APPNAME, d))
            return 1
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.logdir, exist_ok=True)
    jobs = []
    for mod in args.modality:
        for src, opt, pat in TABLES:
            fname = os.path.join(getattr(args, opt), pat.format(PREFIX + mod))
            fname_out = os.path.join(args.outdir, os.path.basename(fname))
            if os.path.abspath(fname) == os.path.abspath(fname_out):
                print('{0}: ERROR: {1} would overwrite its source'.format(
                    APPNAME, fname_out))
                return 1
            if not os.path.isfile(fname):
                print('{0}: ERROR: {1} not found'.format(APPNAME, fname))
                return 1
            if os.path.exists(fname_out) and not args.forceflag:
                print('{0}: ERROR: {1} exists (use --forceflag 1)'.format(
                    APPNAME, fname_out))
                return 1
            jobs.append((mod, src, fname, fname_out))
    vlines = []
    nfail = 0
    visits = {}
    for mod, src, fname, fname_out in jobs:
        stem, types = MODALITIES[mod]
        if mod not in visits:
            fcsv = os.path.join(args.indiv_dir, CSV_PATTERN.format(stem))
            visits[mod] = read_indiv(fcsv, types)
            # the csv consensus columns must be those of its QC_XX columns
            bad = sum(1 for v in visits[mod].values()
                      if v['plain'][0] is not None and v['csv'] != v['plain'])
            print('{0}: {1}: {2} visits, {3} reviewed, {4} where the csv '
                  'consensus differs from its QC_XX columns'.format(
                      APPNAME, os.path.basename(fcsv), len(visits[mod]),
                      sum(1 for v in visits[mod].values()
                          if v['plain'][0] is not None), bad))
            if bad:
                print('{0}: ERROR: csv consensus not reproduced'.format(
                    APPNAME))
                return 1
        edits, notrep, c = edit_table(fname, fname_out, mod, visits[mod])
        tab = '{0}_{1}'.format(mod, src)
        write_rows(os.path.join(args.logdir, 'edits_{0}.tsv'.format(tab)),
                   EDIT_COLUMNS, edits)
        write_rows(os.path.join(args.logdir,
                                'not_reproduced_{0}.tsv'.format(tab)),
                   NOTREP_COLUMNS, notrep)
        print('{0}: {1}: {2}'.format(APPNAME, tab, ', '.join(
            '{0} {1}'.format(k, c[k]) for k in sorted(c))))
        qu = Counter(t for e in edits for t in e[9].split(','))
        print('{0}: {1}: QU types of the edited rows: {2}'.format(
            APPNAME, tab, dict(sorted(qu.items()))))
        if args.verify_flag:
            vc = verify_table(fname, fname_out, mod, edits, visits[mod])
            nf = sum(v for k, v in vc.items() if k.isupper())
            nfail += nf
            vlines += ['{0} {1:40s} {2}'.format(tab, k, vc[k])
                       for k in sorted(vc)]
            vlines.append('{0} failed checks: {1}'.format(tab, nf))
    if args.verify_flag:
        with open(os.path.join(args.logdir, 'verify.txt'), 'w') as f:
            f.write('\n'.join(vlines) + '\n')
        print('\n'.join(vlines))
        if nfail:
            print('{0}: ERROR: verify failed'.format(APPNAME))
            return 1
    return 0


def main():
    sys.exit(run(parse_args()))


if __name__ == '__main__':
    main()
