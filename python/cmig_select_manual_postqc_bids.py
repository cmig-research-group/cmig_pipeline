#!/usr/bin/env python3
"""
Select the sessions for manual post-processing QC (cmig_manual_postqc_bids),
per session and modality (t2, fmri, dmri, fsurf), as a union of layers:

  1 prev_fail   previous failures: overall _score 0 (prev_fail), or _score 1
                with __rev__disp_indicator 1 (prev_disp_pass), in the release
                7.0 tables or the legacy tables (abcd_write_legacy_postqc.py)
  2 rawqc       borderline manual raw QC scores of the scans the proc outputs
                used (their sidecars' SeriesInstanceUID), from the compiled
                manual raw QC series table (--rawqc_series_file): passed
                (score 1) with a QU of 2 or more (rawqc_borderline), or with
                the reviewers split (rawqc_disp_pass); T1w -> fsurf, T2w ->
                t2, bold -> fmri, dwi -> dmri
  3 random      a random sample per modality: --random_fraction of the
                previous-failure rows, at least --random_floor, stratified in
                proportion to event x manufacturer; or the sample in
                --random_sample_file
  4 predict     multivariate prediction of a manual post QC failure per
                modality (cmig_autoqc_select.py): predictors from the
                compiled tables in --tabulated_dir (--predictor_file,
                cmig_autoqc_predictors.py) with the Mahalanobis distances
                (cmig_autoqc_mahal.py, references fit on the sessions that
                are not previous failures), trained on the release 7.0 overall
                QC (--prev_dir) of the sessions with the modality's post QC
                done; those scored above --predict_pthr (uncalibrated) are
                selected (predicted, by p)

Each layer is a function returning rows (participant_id, session_id,
modality, reason, detail); a session and modality picked by several layers
is listed once, with every reason. Rows for visits not in --id_events_file
go to the _dropped file with why.

--filter (as autoqc_select_all's filter) then keeps the rows by their reviews
in --review_rootdir: nrev_done, the reviewers with a QC for the session and
modality, and the consensus _score of the session's
desc-session_manualpostqc file (cmig_manual_postqc_bids, where a reviewer's
pass with an artifact scored 3 counts as a fail). A row is fully reviewed
with nrev_done >= --max_nrev (never with --max_nrev 0).
  new   not fully reviewed (default; a row with one review stays, for the
        next reviewer)
  old   fully reviewed
  any   every row
  pass  fully reviewed, consensus _score 1
  fail  fully reviewed, consensus _score 0
A row without a consensus _score is kept by neither pass nor fail. The rows
the filter leaves out are counted in the summary, not dropped. Of the rows it
keeps, those whose review inputs (those cmig_manual_postqc_bids loads, with a
status errcode 0) are missing go to the _dropped file with why.

Manufacturer (GE, Siemens, Philips) and model come from the compiled
mr_y_adm__info table in --tabulated_dir; the random layer's strata use the
manufacturer.

--max_rows caps the rows selected per modality, over every layer, counting
only the rows --filter keeps. The rows of --prev_selected_file the filter
keeps come first, then the others by priority_score, highest first: +inf
for every layer but predict, a predicted row its p score; ties are broken by
a draw seeded with --random_seed. Rows over the cap go to the _dropped file.

Writes to --outdir:
  manual_postqc_selected_<date>.tsv          one row per session and modality
  manual_postqc_selected_<date>_dropped.tsv  rows left out, and why
  manual_postqc_selected_<date>_summary.txt  counts, with the rows --filter
                                             left out
  manual_postqc_selected_<date>_rawqc_detail.tsv  the scans the rawqc layer
                                                 selected, with their scores
  manual_postqc_selected_<date>_predict/     the predict layer's
      <mod>_model.json (screened metrics, means/SDs, SVD coefs, betas, score
      parms, Mahalanobis references, training ids), <mod>_labels.tsv,
      <mod>_metrics.tsv, <mod>_scores.tsv (every scored session), mahal/
      (the distances and references)

Created:  11/06/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/10/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/18/2026 by Don Hagler
Last mod: 11/22/2026 by Don Hagler
"""

import argparse
import contextlib
import csv
import datetime
import glob
import io
import json
import os
import random
import re
import sys
from collections import Counter, OrderedDict

import cmig_manual_postqc_bids as postqc
import cmig_sing_utils

APPNAME = 'cmig_select_manual_postqc_bids'
INF = float('inf')

# the order of the random draws across modalities
MODALITIES = ['fsurf', 't2', 'dmri', 'fmri']
PREFIX = 'mr_y_qc__post__man__'
KEY = ('participant_id', 'session_id', 'modality')

# reasons in layer order, with the layer number (priority) of each
REASONS = OrderedDict([('prev_fail', 1), ('prev_disp_pass', 1),
                       ('rawqc_borderline', 2), ('rawqc_disp_pass', 2),
                       ('random', 3), ('predicted', 4)])

# priority_score ranks the rows for --max_rows: every layer but predict at
#   +inf, a predicted row at its p score (-inf to inf), so predicted rows
#   rank below the others, by p; ties broken by a seeded random draw
COLUMNS = ['participant_id', 'session_id', 'modality', 'reason', 'priority',
           'priority_score', 'event', 'manufacturer', 'model', 'prev_score',
           'prev_notes', 'detail', 'random_seed', 'nrev_done',
           'pipeline_rootdir']
DROPPED_COLUMNS = ['participant_id', 'session_id', 'modality', 'reason',
                   'why_dropped']

# the compiled manual raw QC series table (cmig_compile_bids_tabdat.py): its
#   column prefix, the QU columns (__<q>_score) and the least QU selected
RAWQC = 'mr_y_qc__raw__man__series'
RAWQC_QU = ['mot', 'susc', 'cutoff', 'line', 'ghost', 'moire', 'famap']
RAWQC_QU_MIN = 2
# the proc sidecars of the scans each review modality's outputs used
USED_SCANS = [('anat', 'T1w', 'fsurf'), ('anat', 'T2w', 't2'),
              ('func', 'bold', 'fmri'), ('dwi', 'dwi', 'dmri')]
# _rawqc_detail.tsv: one row per selected scan, with the table's values
#   (column, table column without the prefix)
RAWQC_DETAIL = ([('series_type', '_type'), ('score', '_score'),
                 ('rev_count', '__rev_count'),
                 ('rev__disp_indicator', '__rev__disp_indicator'),
                 ('score_source', '__score_source')] +
                [(q + '_score', '__{0}_score'.format(q)) for q in RAWQC_QU] +
                [('notes', '_notes')])
RAWQC_DETAIL_COLUMNS = (['participant_id', 'session_id', 'modality', 'scan',
                         'series_uid', 'reason'] +
                        [c for c, _ in RAWQC_DETAIL])

# the compiled table of each session's scanner (--tabulated_dir): manufacturer
#   and model
ADM_INFO = 'mr_y_adm__info'
# the predict layer's predictor list, in $CMIG_PIPELINE_DIR/parms
PREDICTOR_FILE = 'autoqc_predictors.tsv'

#############################################################################
# input tables
#############################################################################

def read_tsv(fname, delimiter='\t'):
    with open(fname, newline='') as f:
        return [dict(r) for r in csv.DictReader(f, delimiter=delimiter)]


def write_tsv(fname, columns, rows):
    with open(fname, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow(dict((c, '' if row.get(c) is None else row.get(c))
                                 for c in columns))


def read_id_events(fname):
    """{(id_redcap, redcap_event_name): (participant_id, session_id)}."""
    ids = {}
    for r in read_tsv(fname, delimiter=','):
        ids[(r['id_redcap'], r['redcap_event_name'])] = (
            r['participant_id'], r['session_id'])
    return ids


def manufacturer_name(m):
    """The short name (GE, Siemens, Philips) of a manufacturer value (e.g.
    'GE MEDICAL SYSTEMS', 'SIEMENS', 'Philips Medical Systems'), another
    value as it is, or 'unknown' if empty."""
    m = clean(m)
    for pre, name in (('siemens', 'Siemens'), ('ge', 'GE'),
                      ('philips', 'Philips')):
        if m.lower().startswith(pre):
            return name
    return m or 'unknown'


def read_adm_info(tabulated_dir):
    """{(pid, sid): {manufacturer, model}} from the compiled ADM_INFO table
    in tabulated_dir (cmig_compile_bids_tabdat.py); the first non-empty
    value of each, where a session has several rows."""
    fname = os.path.join(tabulated_dir, ADM_INFO + '.tsv')
    if not os.path.isfile(fname):
        raise IOError('missing {0}'.format(fname))
    cols = (('manufacturer', ADM_INFO + '__dev_manufact'),
            ('model', ADM_INFO + '__dev_model'))
    info = {}
    with open(fname, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        missing = [c for c in ['participant_id', 'session_id'] +
                   [col for _, col in cols] if c not in (reader.fieldnames or [])]
        if missing:
            raise IOError('{0} lacks column {1}'.format(fname, ', '.join(missing)))
        for r in reader:
            d = info.setdefault((r['participant_id'], r['session_id']), {})
            for k, col in cols:
                if not d.get(k):
                    d[k] = clean(r.get(col))
    for d in info.values():
        d['manufacturer'] = manufacturer_name(d['manufacturer'])
    return info

#############################################################################
# review inputs
#############################################################################

# the input roots each modality needs (as cmig_manual_postqc_bids)
ROOTS = {'t2': ['proc_rootdir'], 'fmri': ['proc_rootdir', 'postqc_rootdir'],
         'dmri': ['proc_rootdir', 'analysis_rootdir'],
         'fsurf': ['fsurf_rootdir']}


def session_dir(root_kind, root, pid, sid):
    if root_kind == 'fsurf_rootdir':
        return os.path.join(root, '{0}_{1}'.format(pid, sid))
    return os.path.join(root, pid, sid)


def find_root(args, root_kind, pid, sid):
    """The first of the roots given for root_kind that holds the session."""
    for root in getattr(args, root_kind) or []:
        if os.path.isdir(session_dir(root_kind, root, pid, sid)):
            return root
    return None


def pipeline_root(args, pid, sid, mod):
    """The cmig_pipeline root holding all of the session's roots for the
    modality (each named as postqc.PIPELINE_SUBDIRS has it), else ''."""
    parents = set()
    for kind in ROOTS[mod]:
        root = find_root(args, kind, pid, sid)
        if root is None:
            return ''
        root = os.path.normpath(root)
        if os.path.basename(root) != postqc.PIPELINE_SUBDIRS[kind]:
            return ''
        parents.add(os.path.dirname(root))
    return parents.pop() if len(parents) == 1 else ''


def status_errmsg(fname):
    """'' if the status file's last row has errcode 0, else why not."""
    if not os.path.isfile(fname):
        return 'no status file {0}'.format(os.path.basename(fname))
    try:
        rows = read_tsv(fname)
    except (IOError, OSError, csv.Error) as e:
        return 'cannot read {0}: {1}'.format(os.path.basename(fname), e)
    if not rows:
        return 'empty status file {0}'.format(os.path.basename(fname))
    if (rows[-1].get('errcode') or '').strip() != '0':
        return 'status errcode {0} in {1}: {2}'.format(
            rows[-1].get('errcode'), os.path.basename(fname),
            rows[-1].get('errmsg') or '')
    return ''


def check_inputs(args, pid, sid, mod):
    """'' if the review inputs of the session and modality exist, with a
    status errcode 0, else why not."""
    stem = '{0}_{1}'.format(pid, sid)
    roots = {}
    for kind in ROOTS[mod]:
        if not getattr(args, kind):
            return 'no --{0} given'.format(kind)
        roots[kind] = find_root(args, kind, pid, sid)
        if roots[kind] is None:
            return 'not in --{0}'.format(kind)
    # the inputs cmig_manual_postqc_bids loads, found as it finds them
    targs = argparse.Namespace(modality=mod, proc_rootdir=roots.get('proc_rootdir'),
                               postqc_rootdir=roots.get('postqc_rootdir'),
                               analysis_rootdir=roots.get('analysis_rootdir'),
                               fsurf_rootdir=roots.get('fsurf_rootdir'))
    with contextlib.redirect_stdout(io.StringIO()):
        inputs, errmsg = postqc.get_inputs(targs, pid, sid)
    if inputs is None:
        return errmsg
    status = []
    if mod == 'fsurf':
        status.append(os.path.join(roots['fsurf_rootdir'], stem,
                                   stem + '_desc-fsurf_status.tsv'))
    else:
        status.append(os.path.join(roots['proc_rootdir'], pid, sid, 'anat',
                                   stem + '_desc-proc_status.tsv'))
    if mod == 'fmri':
        status.append(os.path.join(roots['postqc_rootdir'], pid, sid, 'func',
                                   stem + '_desc-postqc_status.tsv'))
    elif mod == 'dmri':
        status.append(os.path.join(roots['analysis_rootdir'], pid, sid, 'dwi',
                                   stem + '_desc-analyze_status.tsv'))
    for fname in status:
        errmsg = status_errmsg(fname)
        if errmsg:
            return errmsg
    return ''


def processed_sessions(args, mod):
    """(pid, sid) of the sessions with the modality's directory in its
    roots (as cmig_manual_postqc_bids finds them), sorted."""
    sessions = set()
    if mod == 'fsurf':
        for root in args.fsurf_rootdir or []:
            for d in glob.glob(os.path.join(root, 'sub-*_ses-*')):
                name = os.path.basename(d)
                idx = name.find('_ses-')
                if os.path.isdir(d) and idx > 0:
                    sessions.add((name[:idx], name[idx + 1:]))
        return sorted(sessions)
    kind, subdir = {'t2': ('proc_rootdir', 'anat'),
                    'fmri': ('postqc_rootdir', 'func'),
                    'dmri': ('analysis_rootdir', 'dwi')}[mod]
    for root in getattr(args, kind) or []:
        for d in glob.glob(os.path.join(root, 'sub-*', 'ses-*', subdir)):
            if os.path.isdir(d):
                sessions.add((os.path.basename(os.path.dirname(os.path.dirname(d))),
                              os.path.basename(os.path.dirname(d))))
    return sorted(sessions)


def count_reviews(review_rootdir, pid, sid, mod):
    """Reviewers with a QC for the session and modality in the review root,
    from the per-modality files and an old single file alike."""
    if not review_rootdir:
        return 0
    rows = postqc.read_review_rows(postqc.output_stem(review_rootdir, pid, sid))
    return len(set(r.get('reviewer_id') for r in rows
                   if r.get('modality') == mod and
                   (r.get('QC') or '').strip() != ''))


def consensus_score(review_rootdir, pid, sid, mod):
    """The consensus _score (1 pass, 0 fail) of the session and modality in
    its desc-session_manualpostqc file in the review root, else None."""
    if not review_rootdir:
        return None
    rows = postqc.read_tsv_rows(postqc.output_stem(review_rootdir, pid, sid) +
                                '_desc-session_manualpostqc.tsv')
    return num(rows[0].get(PREFIX + mod + '_score')) if rows else None


# --filter values, as autoqc_select_all's filter
FILTERS = ['new', 'old', 'any', 'pass', 'fail']


def filter_keep(filt, nrev_done, max_nrev, score=None):
    """Whether --filter keeps a row: new not fully reviewed (nrev_done <
    max_nrev, every row with max_nrev 0), old fully reviewed, any every row,
    pass / fail fully reviewed with consensus score 1 / 0 (None: neither)."""
    full = max_nrev > 0 and nrev_done >= max_nrev
    if filt == 'any':
        return True
    if filt == 'new':
        return not full
    if filt == 'old':
        return full
    return full and score == {'pass': 1, 'fail': 0}[filt]

#############################################################################
# layers
#############################################################################

def clean(val):
    val = (val or '').strip()
    return '' if val == 'n/a' else val


def layer_prev_fail(args, ctx):
    """Rows of the sessions that failed manual post QC before (prev_fail),
    or passed with reviewers split (prev_disp_pass), in either the release
    7.0 tables (--prev_dir) or the legacy tables (--legacy_dir)."""
    sources = [('release7.0', args.prev_dir, '{0}.tsv')]
    if args.legacy_dir:
        sources.append(('legacy', args.legacy_dir, '{0}_desc-legacy.tsv'))
    rows = []
    for mod in MODALITIES:
        pre = PREFIX + mod
        fail = OrderedDict()
        disp = OrderedDict()
        for src, dirname, pat in sources:
            fname = os.path.join(dirname, pat.format(pre))
            if not os.path.isfile(fname):
                raise IOError('missing {0}'.format(fname))
            for r in read_tsv(fname):
                key = (r['participant_id'], r['session_id'])
                score = clean(r.get(pre + '_score'))
                notes = clean(r.get(pre + '_notes'))
                if score == '0':
                    fail.setdefault(key, []).append((src, notes))
                elif score == '1' and clean(
                        r.get(pre + '__rev__disp_indicator')) == '1':
                    disp.setdefault(key, []).append((src, notes))
        for reason, score, picks in (('prev_fail', 0, fail),
                                     ('prev_disp_pass', 1, disp)):
            for (pid, sid), srcs in sorted(picks.items()):
                if reason == 'prev_disp_pass' and (pid, sid) in fail:
                    continue
                rows.append({'participant_id': pid, 'session_id': sid,
                             'modality': mod, 'reason': reason,
                             'detail': '{0} in {1}'.format(
                                 reason, ','.join(s for s, _ in srcs)),
                             'prev_score': score,
                             'prev_notes': next((n for _, n in srcs if n), '')})
    return rows


def num(val):
    """A table value as a number, or None (empty, n/a, not a number)."""
    try:
        return float(clean(val))
    except ValueError:
        return None


def fmt_num(val):
    return str(int(val)) if val == int(val) else str(val)


def read_rawqc_series(fname):
    """{(pid, sid): [rows]} of the compiled manual raw QC series table: the
    passed series (score 1) with a QU of RAWQC_QU_MIN or more or the
    reviewers split (rev__disp_indicator 1)."""
    picks = {}
    for r in read_tsv(fname):
        if num(r.get(RAWQC + '_score')) != 1:
            continue
        if num(r.get(RAWQC + '__rev__disp_indicator')) == 1 or any(
                (num(r.get(RAWQC + '__{0}_score'.format(q))) or 0) >= RAWQC_QU_MIN
                for q in RAWQC_QU):
            picks.setdefault((r['participant_id'], r['session_id']), []).append(r)
    return picks


def used_series(args, pid, sid):
    """{SeriesInstanceUID: (scan, review modality)} of the scans the
    session's proc outputs used, from their desc-preproc sidecars, else
    those of the older names without it (scan: the sidecar's name without
    the participant, session and desc-preproc, e.g. run-01_T1w)."""
    root = find_root(args, 'proc_rootdir', pid, sid)
    used = {}
    if root is None:
        return used
    stem = '{0}_{1}_'.format(pid, sid)
    for subdir, suffix, mod in USED_SCANS:
        fnames = sorted(glob.glob(os.path.join(
            root, pid, sid, subdir, stem + '*_desc-preproc_' + suffix + '.json')))
        if not fnames:
            fnames = [f for f in sorted(glob.glob(os.path.join(
                root, pid, sid, subdir, stem + '*_' + suffix + '.json')))
                if '_desc-' not in os.path.basename(f)]
        for fname in fnames:
            name = os.path.basename(fname).replace('_desc-preproc_', '_')
            if '_space-' in name:
                continue
            try:
                with open(fname) as f:
                    uid = json.load(f).get('SeriesInstanceUID')
            except (IOError, OSError, ValueError):
                continue
            if uid:
                used[uid] = (name[len(stem):-len('.json')], mod)
    return used


def layer_rawqc(args, ctx):
    """Rows of the sessions whose used scans (the proc sidecars' series)
    passed manual raw QC with a QU of RAWQC_QU_MIN or more
    (rawqc_borderline) or with the reviewers split (rawqc_disp_pass), from
    the compiled manual raw QC series table (--rawqc_series_file). The T1w
    feeds the fsurf review, the T2w t2, bold fmri, dwi dmri. Each selected
    scan goes to ctx['rawqc_detail']."""
    picks = read_rawqc_series(args.rawqc_series_file)
    rows = []
    detail = []
    n_used = 0
    for (pid, sid), series in sorted(picks.items()):
        used = used_series(args, pid, sid)
        for r in series:
            uid = (r.get(RAWQC + '_uid') or '').strip()
            if uid not in used:
                continue
            n_used += 1
            scan, mod = used[uid]
            high = [(q, num(r.get(RAWQC + '__{0}_score'.format(q))))
                    for q in RAWQC_QU]
            high = [(q, v) for q, v in high if v is not None and v >= RAWQC_QU_MIN]
            reasons = []
            if high:
                reasons.append(('rawqc_borderline', '{0} {1}'.format(
                    scan, ' '.join('{0}={1}'.format(q, fmt_num(v))
                                   for q, v in high))))
            if num(r.get(RAWQC + '__rev__disp_indicator')) == 1:
                reasons.append(('rawqc_disp_pass', '{0} reviewers split'
                                .format(scan)))
            for reason, text in reasons:
                rows.append({'participant_id': pid, 'session_id': sid,
                             'modality': mod, 'reason': reason,
                             'detail': text})
            d = {'participant_id': pid, 'session_id': sid, 'modality': mod,
                 'scan': scan, 'series_uid': uid,
                 'reason': ';'.join(x for x, _ in reasons)}
            for c, col in RAWQC_DETAIL:
                d[c] = clean(r.get(RAWQC + col))
            detail.append(d)
    ctx['rawqc_detail'] = detail
    ctx['rawqc_counts'] = (sum(len(v) for v in picks.values()), n_used)
    return rows


def allocate(counts, n):
    """Proportional allocation of n over the cells (an OrderedDict of cell
    sizes), by largest remainder."""
    total = sum(counts.values())
    quota = OrderedDict((k, n * c / total) for k, c in counts.items())
    alloc = OrderedDict((k, min(int(q), counts[k])) for k, q in quota.items())
    left = n - sum(alloc.values())
    order = sorted(quota, key=lambda k: quota[k] - int(quota[k]), reverse=True)
    i = 0
    while left > 0 and i < 10 * len(order):
        k = order[i % len(order)]
        if alloc[k] < counts[k]:
            alloc[k] += 1
            left -= 1
        i += 1
    return alloc


def draw_random_sample(population, n_prev, strata, fraction, floor, seed):
    """(rows, sizes): per modality, in MODALITIES order with one generator,
    max(round(fraction x n_prev), floor) sessions (at most all), drawn from
    population[mod] (sorted (pid, sid)) in proportion to the strata cells
    (strata(pid, sid) gives (event, manufacturer))."""
    rng = random.Random(seed)
    rows = []
    sizes = OrderedDict()
    for mod in MODALITIES:
        pop = sorted(population.get(mod, []))
        n = min(max(int(round(fraction * n_prev.get(mod, 0))), floor), len(pop))
        sizes[mod] = (n_prev.get(mod, 0), n, len(pop))
        if not n:
            continue
        cells = OrderedDict()
        for key in pop:
            cells.setdefault(strata(*key), []).append(key)
        cells = OrderedDict((k, cells[k]) for k in sorted(cells))
        alloc = allocate(OrderedDict((k, len(g)) for k, g in cells.items()), n)
        for k, g in cells.items():
            for i in sorted(rng.sample(range(len(g)), alloc[k])):
                pid, sid = g[i]
                rows.append({'participant_id': pid, 'session_id': sid,
                             'modality': mod, 'reason': 'random',
                             'detail': 'random {0} {1}'.format(*k),
                             'random_seed': seed})
    return rows, sizes


def read_sampled(fnames):
    """{(pid, sid, modality)} sampled before: rows of selection files with
    reason random, or every row of a file without a reason column (a saved
    random sample)."""
    sampled = set()
    for fname in fnames or []:
        for r in read_tsv(fname):
            if 'reason' in r and 'random' not in (r['reason'] or '').split(';'):
                continue
            sampled.add((r['participant_id'], r['session_id'], r['modality']))
    return sampled


def layer_random(args, ctx):
    """Rows of the random sample: read from --random_sample_file, else drawn
    from the sessions whose review inputs exist, sessions sampled in --exclude_sampled_files
    left out."""
    if args.random_sample_file:
        rows = []
        for r in read_tsv(args.random_sample_file):
            seed = r.get('random_seed') or r.get('seed') or ''
            rows.append({'participant_id': r['participant_id'],
                         'session_id': r['session_id'],
                         'modality': r['modality'], 'reason': 'random',
                         'detail': 'random, from {0}'.format(
                             os.path.basename(args.random_sample_file)),
                         'random_seed': seed})
        ctx['random_sizes'] = None
        return rows
    n_prev = Counter(r['modality'] for r in ctx['rows']
                     if r['reason'] in ('prev_fail', 'prev_disp_pass'))
    sampled = read_sampled(args.exclude_sampled_files)
    # processed: the review inputs exist, so no draw is dropped later
    population = {}
    for mod in MODALITIES:
        population[mod] = [(pid, sid) for pid, sid in ctx['processed'](mod)
                           if (pid, sid) in ctx['visits'] and
                           (pid, sid, mod) not in sampled and
                           not check_inputs(args, pid, sid, mod)]
    info = ctx['adm_info']

    def strata(pid, sid):
        return (sid, info.get((pid, sid), {}).get('manufacturer', 'unknown'))

    rows, ctx['random_sizes'] = draw_random_sample(
        population, n_prev, strata, args.random_fraction, args.random_floor,
        args.random_seed)
    return rows


def predict_modality(args, ctx, mod, mahal_file, refs, excl, outdir):
    """(rows, scores) of the predict layer for one modality; writes
    <mod>_model.json, <mod>_labels.tsv, <mod>_metrics.tsv, <mod>_scores.tsv
    to outdir."""
    import math
    import numpy as np
    import cmig_autoqc_predictors as pred
    import cmig_autoqc_select as asel
    ids = pred.postqc_ids(args.tabulated_dir, mod, args.predictor_file)
    ids, X, names, info = pred.load_predictors(
        args.tabulated_dir, mod, args.predictor_file, ids=ids,
        mahal_file=mahal_file, strict=False)
    labels = pred.load_labels(args.prev_dir, mod)
    y = np.array([labels.get(k, float('nan')) for k in ids], dtype=float)
    idx_train = np.nonzero(~np.isnan(y))[0]
    # predictors missing for most training sessions are left out
    keep = pred.screen_missing(X, idx_train, args.predict_max_frac_nan)
    screened_out = [n for j, n in enumerate(names) if j not in keep]
    X = X[:, keep]
    names = [names[j] for j in keep]
    # 8.0 reviews, not the release 7.0 labels, decide what is selectable
    allowed = [i for i, (pid, sid) in enumerate(ids)
               if args.max_nrev <= 0 or count_reviews(
                   args.review_rootdir, pid, sid, mod) < args.max_nrev]
    manuf = [ctx['adm_info'].get(k, {}).get('manufacturer', 'unknown')
             for k in ids]
    lname = pred.LABEL_SPEC[mod]['column']
    parms = {'pthr': args.predict_pthr, 'filter': 'any',
             'scanner_flag': args.predict_scanner_flag,
             'min_n': args.predict_min_n, 'min_nfail': args.predict_min_nfail}
    log = (lambda m: print('{0}: {1}: {2}'.format(APPNAME, mod, m))) \
        if args.verbose else None
    res = asel.train(X, names, y, [{'name': lname, 'type': 'pass1'}],
                     ids=['{0}_{1}'.format(*k) for k in ids], parms=parms,
                     idx_train=idx_train, idx_allowed=allowed,
                     manufacturer=manuf, log=log)
    model = res['model']
    model['modality'] = mod
    model['predictor_file'] = args.predictor_file
    model['predictor_fill'] = info['fill']
    model['predictors_screened_out'] = screened_out
    model['max_frac_nan'] = args.predict_max_frac_nan
    model['n_missing_task'] = info['n_missing_task']
    model['mahal_reference'] = dict(
        (fam, ref) for fam, ref in refs.items()
        if any(n.startswith(fam + '_') for n in names))
    model['mahal_ref_excluded'] = sorted('{0}_{1}'.format(p, s)
                                         for p, s, m in excl if m == mod)
    asel.model_to_json(model, os.path.join(outdir, mod + '_model.json'))
    asel.write_rows(res['labels_rows'], asel.LABEL_COLS,
                    os.path.join(outdir, mod + '_labels.tsv'))
    asel.write_rows(res['metrics_rows'], asel.SCREEN_COLS,
                    os.path.join(outdir, mod + '_metrics.tsv'))
    ps = res['pscore'].get(lname)
    yp = res['yp'].get(lname)
    scores = []
    allowed_set = set(allowed)
    for i, (pid, sid) in enumerate(ids):
        scores.append({'participant_id': pid, 'session_id': sid,
                       'modality': mod, 'label': y[i],
                       'yp': yp[i] if yp is not None else math.nan,
                       'p': ps[i] if ps is not None else math.nan,
                       'allowed': 1 if i in allowed_set else 0,
                       'train': 1 if not math.isnan(y[i]) else 0})
    asel.write_rows(scores, ['participant_id', 'session_id', 'modality',
                             'label', 'yp', 'p', 'allowed', 'train'],
                    os.path.join(outdir, mod + '_scores.tsv'))
    rows = []
    if ps is None:
        status = model['labels'][0]['status']
        ctx.setdefault('predict_counts', OrderedDict())[mod] = (
            len(ids), idx_train.size, len(names), status, 0)
        return rows
    sel = asel.select(res['pscore'], y[:, None], model, idx_allowed=allowed,
                      yp=res['yp'])
    for r in sel['rows']:
        pid, sid = ids[r['row']]
        rows.append({'participant_id': pid, 'session_id': sid,
                     'modality': mod, 'reason': 'predicted',
                     'detail': 'predicted {0} p={1:.3g}'.format(
                         lname[len(PREFIX):], r['p']),
                     'p_score': r['p']})
    ctx.setdefault('predict_counts', OrderedDict())[mod] = (
        len(ids), idx_train.size, len(names), model['labels'][0]['status'],
        len(rows))
    return rows


def layer_predict(args, ctx):
    """Rows of the sessions predicted to fail manual post QC, per modality:
    Mahalanobis distances with references fit on the sessions that are not
    previous failures (as layer_prev_fail defines them), the predictors of
    --predictor_file, a model trained on the release 7.0 overall QC of the
    sessions with the modality's post QC done, and those scored above
    --predict_pthr selected; the model, its records and the distances are
    saved in <outdir>/manual_postqc_selected_<date>_predict/."""
    import cmig_autoqc_mahal as mahal
    outdir = ctx['predict_dir']
    excl = sorted(set((r['participant_id'], r['session_id'], r['modality'])
                      for r in layer_prev_fail(args, {})
                      if r['reason'] == 'prev_fail'))
    mdir = os.path.join(outdir, 'mahal')
    with contextlib.redirect_stdout(io.StringIO()):
        infos = mahal.run_families_80(
            args.tabulated_dir, mdir, ref_exclude=excl,
            ref_source='previous failures (layer_prev_fail)',
            scanner_flag=1, n_draws=args.predict_n_draws, verbose=False)
    errs = ['{0}: {1}'.format(k, v['error']) for k, v in infos.items()
            if v.get('error')]
    if errs and args.verbose:
        print('{0}: WARNING: Mahalanobis families failed: {1}'.format(
            APPNAME, '; '.join(errs)))
    refs = {}
    for fam in infos:
        fname = mahal.reference_fname(mdir, fam)
        if os.path.isfile(fname):
            refs[fam] = mahal.load_reference(fname)
    rows = []
    for mod in MODALITIES:
        rows += predict_modality(args, ctx, mod,
                                 os.path.join(mdir, 'mahal_mdist.tsv'), refs,
                                 excl, outdir)
    return rows


# (option, function), in layer order
LAYERS = [('prev_fail_flag', layer_prev_fail), ('rawqc_flag', layer_rawqc),
          ('random_flag', layer_random), ('predict_flag', layer_predict)]

#############################################################################
# union and output
#############################################################################

def merge_rows(rows):
    """One row per (pid, sid, modality), every reason in layer order."""
    merged = OrderedDict()
    for r in rows:
        key = tuple(r[k] for k in KEY)
        m = merged.setdefault(key, {'participant_id': key[0],
                                    'session_id': key[1], 'modality': key[2],
                                    'reasons': [], 'details': []})
        if r['reason'] not in m['reasons']:
            m['reasons'].append(r['reason'])
        if r.get('detail'):
            m['details'].append(r['detail'])
        for k in ('prev_score', 'prev_notes', 'random_seed'):
            if r.get(k) not in (None, '') and m.get(k) in (None, ''):
                m[k] = r[k]
        if r.get('p_score') is not None:
            m['p_score'] = max(m.get('p_score', -INF), float(r['p_score']))
    for m in merged.values():
        m['reasons'].sort(key=list(REASONS).index)
        m['reason'] = ';'.join(m['reasons'])
        m['priority'] = min(REASONS[x] for x in m['reasons'])
        m['priority_score'] = priority_score(m)
        m['detail'] = '; '.join(m['details'])
    return list(merged.values())


def priority_score(m):
    """+inf if a layer other than predict picked the row, else its p score
    (-inf to inf; -inf if the predict layer gave none)."""
    if any(x != 'predicted' for x in m['reasons']):
        return INF
    return m.get('p_score', -INF)


def read_prev_selected(fname):
    """{(pid, sid, modality)} of a previous selection file."""
    rows = read_tsv(fname)
    if rows and 'modality' not in rows[0]:
        raise IOError('{0} lacks a modality column'.format(fname))
    return set(tuple(r[k] for k in KEY) for r in rows)


def apply_cap(args, kept, ctx):
    """(listed, over): per modality, the rows (those --filter kept) capped at
    --max_rows: first those in --prev_selected_file, then by priority_score,
    highest first, ties broken by a draw seeded with --random_seed (a
    generator of its own). ctx['cap_counts'][mod]: rows, kept from the
    previous selection, added, over."""
    prev = read_prev_selected(args.prev_selected_file) \
        if args.prev_selected_file else set()
    rng = random.Random('{0}_max_rows'.format(args.random_seed))
    listed, over = [], []
    counts = OrderedDict()
    for mod in MODALITIES:
        rows = sorted((m for m in kept if m['modality'] == mod),
                      key=lambda m: tuple(m[k] for k in KEY))
        draw = dict((tuple(m[k] for k in KEY), rng.random()) for m in rows)

        def rank(m):
            key = tuple(m[k] for k in KEY)
            return (key not in prev, -m['priority_score'], draw[key])

        rows.sort(key=rank)
        take = rows[:args.max_rows]
        n_prev = sum(1 for m in take if tuple(m[k] for k in KEY) in prev)
        counts[mod] = (len(rows), n_prev, len(take) - n_prev,
                       len(rows) - len(take))
        listed += take
        over += rows[args.max_rows:]
    ctx['cap_counts'] = counts
    return listed, over


def short_why(why):
    """why_dropped with paths and file lists left out, for the counts."""
    why = re.sub(r' in /\S+', '', why)
    why = re.sub(r'(missing recon files): .*', r'\1', why)
    why = re.sub(r'(status errcode \S+) in \S+: .*', r'\1', why)
    return re.sub(r'(no status file) \S+', r'\1', why)


def tally(rows, field):
    """'mod: {value: n, ...}' lines of rows counted by modality and field."""
    lines = []
    for mod in MODALITIES:
        c = Counter(str(r.get(field, '')) for r in rows if r['modality'] == mod)
        lines.append('  {0}: {1}'.format(mod, dict(sorted(c.items()))))
    return lines


def write_summary(fname, args, rows, kept, dropped, ctx):
    lines = ['{0} {1}'.format(APPNAME, args.date)]
    lines.append('layers: ' + ', '.join(
        '{0} {1}'.format(opt, getattr(args, opt)) for opt, _ in LAYERS))
    if args.random_flag:
        if args.random_sample_file:
            lines.append('random sample read from {0}'.format(
                args.random_sample_file))
        else:
            lines.append('random sample: seed {0}, fraction {1}, floor {2}'
                         .format(args.random_seed, args.random_fraction,
                                 args.random_floor))
            lines.append('random sample per modality (prev rows, drawn, '
                         'population):')
            for mod, size in (ctx.get('random_sizes') or {}).items():
                lines.append('  {0}: {1}'.format(mod, size))
    if args.rawqc_flag:
        lines.append('rawqc: {0}; passed series with a QU >= {1} or reviewers '
                     'split: {2}, used in proc: {3}'.format(
                         args.rawqc_series_file, RAWQC_QU_MIN,
                         *ctx['rawqc_counts']))
    if args.predict_flag:
        lines.append('predict: {0}, pthr {1}, scanner_flag {2}, max_frac_nan '
                     '{3}, min_n {4}, min_nfail {5}; per modality (sessions '
                     'scored, trained, predictors, status, rows):'.format(
                         args.predictor_file, args.predict_pthr,
                         args.predict_scanner_flag, args.predict_max_frac_nan,
                         args.predict_min_n, args.predict_min_nfail))
        for mod, n in (ctx.get('predict_counts') or {}).items():
            lines.append('  {0}: {1}'.format(mod, n))
    lines.append('layer rows per modality and reason, before drops:')
    lines += tally(rows, 'reason')
    lines.append('filter {0} (max_nrev {1}); rows kept, left out:'.format(
        args.filter, args.max_nrev))
    for mod, n in ctx['filter_counts'].items():
        lines.append('  {0}: {1}'.format(mod, n))
    if args.max_rows is not None:
        lines.append('max_rows {0} per modality (random_seed {1}){2}; rows '
                     'the filter kept with their inputs, kept from the '
                     'previous selection, added, over the cap:'.format(
                         args.max_rows, args.random_seed,
                         ', previous selection {0}'.format(
                             args.prev_selected_file)
                         if args.prev_selected_file else ''))
        for mod, n in ctx['cap_counts'].items():
            lines.append('  {0}: {1}'.format(mod, n))
    lines.append('selected (session, modality) rows: {0}, sessions {1}'.format(
        len(kept), len(set((r['participant_id'], r['session_id'])
                           for r in kept))))
    lines.append('selected per modality and reason (overlaps as a;b):')
    lines += tally(kept, 'reason')
    lines.append('selected per modality and event:')
    lines += tally(kept, 'event')
    lines.append('selected per modality and manufacturer:')
    lines += tally(kept, 'manufacturer')
    lines.append('selected per modality and nrev_done:')
    lines += tally(kept, 'nrev_done')
    lines.append('dropped rows: {0}'.format(len(dropped)))
    c = Counter((r['modality'], short_why(r['why_dropped'])) for r in dropped)
    for (mod, why), n in sorted(c.items()):
        lines.append('  {0}: {1}: {2}'.format(mod, why, n))
    with open(fname, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    return lines


def run(args):
    if args.predict_flag and not args.tabulated_dir:
        return 1, '--predict_flag 1 needs --tabulated_dir'
    if args.predict_flag and not args.predictor_file:
        try:
            args.predictor_file = os.path.join(
                cmig_sing_utils.cmig_getenv('CMIG_PIPELINE_DIR'),
                'parms', PREDICTOR_FILE)
        except RuntimeError as e:
            return 1, str(e)
    if args.predict_flag and not os.path.isfile(args.predictor_file):
        return 1, 'missing --predictor_file {0}'.format(args.predictor_file)
    if args.rawqc_flag and not args.rawqc_series_file:
        return 1, '--rawqc_flag 1 needs --rawqc_series_file'
    if args.rawqc_flag and not os.path.isfile(args.rawqc_series_file):
        return 1, 'missing --rawqc_series_file {0}'.format(args.rawqc_series_file)
    if args.max_rows is not None and args.max_rows < 0:
        return 1, '--max_rows must be 0 or more'
    if args.prev_selected_file and not os.path.isfile(args.prev_selected_file):
        return 1, 'missing --prev_selected_file {0}'.format(
            args.prev_selected_file)
    if args.tabulated_dir and not os.path.isdir(args.tabulated_dir):
        return 1, 'missing --tabulated_dir {0}'.format(args.tabulated_dir)
    if args.random_flag and args.random_sample_file and \
            args.exclude_sampled_files:
        return 1, '--random_sample_file excludes --exclude_sampled_files'
    out = dict((k, os.path.join(args.outdir, 'manual_postqc_selected_{0}{1}'
                                .format(args.date, ext)))
               for k, ext in (('sel', '.tsv'), ('drop', '_dropped.tsv'),
                              ('sum', '_summary.txt'),
                              ('rawqc', '_rawqc_detail.tsv'))
               if k != 'rawqc' or args.rawqc_flag)
    if args.predict_flag:
        out['predict'] = os.path.join(
            args.outdir, 'manual_postqc_selected_{0}_predict'.format(args.date))
    if not args.forceflag and any(os.path.exists(f) for f in out.values()):
        return 1, 'outputs exist in {0}; --forceflag 1 to overwrite'.format(
            args.outdir)

    ids = read_id_events(args.id_events_file)
    visits = set(ids.values())
    adm_info = {}
    if args.tabulated_dir:
        if args.verbose:
            print('{0}: reading {1}.tsv in {2}'.format(
                APPNAME, ADM_INFO, args.tabulated_dir))
        try:
            adm_info = read_adm_info(args.tabulated_dir)
        except IOError as e:
            return 1, str(e)
    cache = {}

    def processed(mod):
        if mod not in cache:
            cache[mod] = processed_sessions(args, mod)
        return cache[mod]

    ctx = {'visits': visits, 'adm_info': adm_info, 'processed': processed,
           'rows': []}
    if args.predict_flag:
        ctx['predict_dir'] = os.path.join(
            args.outdir, 'manual_postqc_selected_{0}_predict'.format(args.date))
        os.makedirs(ctx['predict_dir'], exist_ok=True)
    for opt, layer in LAYERS:
        if not getattr(args, opt):
            continue
        rows = layer(args, ctx)
        if args.verbose:
            print('{0}: {1}: {2} rows'.format(APPNAME, layer.__name__, len(rows)))
        # visits not in the id events file are dropped before later layers
        #   count the rows
        ctx['rows'] += [r for r in rows
                        if (r['participant_id'], r['session_id']) in visits]
        ctx.setdefault('not_visits', []).extend(
            r for r in rows
            if (r['participant_id'], r['session_id']) not in visits)

    dropped = [dict(m, why_dropped='not in {0}'.format(
        os.path.basename(args.id_events_file)))
        for m in merge_rows(ctx.get('not_visits', []))]
    kept = []
    nkeep, nout = Counter(), Counter()
    for m in merge_rows(ctx['rows']):
        pid, sid, mod = m['participant_id'], m['session_id'], m['modality']
        m['nrev_done'] = count_reviews(args.review_rootdir, pid, sid, mod)
        score = consensus_score(args.review_rootdir, pid, sid, mod) \
            if args.filter in ('pass', 'fail') and m['nrev_done'] else None
        if not filter_keep(args.filter, m['nrev_done'], args.max_nrev, score):
            nout[mod] += 1
            continue
        nkeep[mod] += 1
        why = check_inputs(args, pid, sid, mod)
        if why:
            dropped.append(dict(m, why_dropped=why))
            continue
        info = adm_info.get((pid, sid), {})
        m.update({'event': sid,
                  'manufacturer': info.get('manufacturer', 'unknown'),
                  'model': info.get('model', ''),
                  'pipeline_rootdir': pipeline_root(args, pid, sid, mod)})
        kept.append(m)
    ctx['filter_counts'] = OrderedDict((mod, (nkeep[mod], nout[mod]))
                                       for mod in MODALITIES)
    if args.verbose:
        print('{0}: filter {1}: rows left out: {2}'.format(
            APPNAME, args.filter, sum(nout.values())))
    if args.max_rows is not None:
        try:
            kept, over = apply_cap(args, kept, ctx)
        except IOError as e:
            return 1, str(e)
        dropped += [dict(m, why_dropped='over the cap (--max_rows {0})'.format(
            args.max_rows)) for m in over]
    kept.sort(key=lambda m: (m['priority'], m['participant_id'],
                             m['session_id'], MODALITIES.index(m['modality'])))
    os.makedirs(args.outdir, exist_ok=True)
    write_tsv(out['sel'], COLUMNS, kept)
    write_tsv(out['drop'], DROPPED_COLUMNS, dropped)
    if args.rawqc_flag:
        write_tsv(out['rawqc'], RAWQC_DETAIL_COLUMNS, ctx['rawqc_detail'])
    lines = write_summary(out['sum'], args, ctx['rows'], kept, dropped, ctx)
    if args.verbose:
        print('\n'.join(lines))
    return 0, ''

#############################################################################
# CLI
#############################################################################

def parse_args(argv=None):
    today = datetime.date.today().strftime('%y%m%d')
    p = argparse.ArgumentParser(
        description='Select sessions for manual post-processing QC')
    p.add_argument('--outdir', required=True,
                   help='where the selection files go')
    p.add_argument('--date', default=today,
                   help='YYMMDD in the file names (default: today)')
    for kind, what in (('proc_rootdir', 'cmig_proc'),
                       ('postqc_rootdir', 'cmig_postqc'),
                       ('analysis_rootdir', 'cmig_analysis'),
                       ('fsurf_rootdir', 'FreeSurfer')):
        p.add_argument('--' + kind, nargs='+',
                       help='{0} roots, searched in turn for each session'
                            .format(what))
    p.add_argument('--review_rootdir', required=True,
                   help='cmig_manual_postqc_bids output root, for nrev_done')
    p.add_argument('--id_events_file', required=True,
                   help='csv of the visits to select from (id_redcap, '
                        'redcap_event_name, participant_id, session_id)')
    p.add_argument('--tabulated_dir',
                   help='directory of the compiled tabulated tables '
                        '(cmig_compile_bids_tabdat.py): manufacturer and '
                        'model from {0}.tsv (default: none, manufacturer '
                        'unknown)'.format(ADM_INFO))
    for opt, _ in LAYERS:
        default = 1
        p.add_argument('--' + opt, type=int, default=default, choices=[0, 1],
                       help='run the {0} layer (default: {1})'.format(
                           opt[:-5], default))
    p.add_argument('--prev_dir', required=True,
                   help='release 7.0 {0}<modality>.tsv tables'.format(PREFIX))
    p.add_argument('--legacy_dir',
                   help='legacy {0}<modality>_desc-legacy.tsv tables '
                        '(default: none)'.format(PREFIX))
    p.add_argument('--rawqc_series_file',
                   help='compiled manual raw QC series table ({0}.tsv, '
                        'cmig_compile_bids_tabdat.py), for the rawqc layer'
                        .format(RAWQC))
    p.add_argument('--predictor_file',
                   help='predictor list of the predict layer (default: '
                        '$CMIG_PIPELINE_DIR/parms/{0})'.format(PREDICTOR_FILE))
    p.add_argument('--predict_pthr', type=float, default=0.0,
                   help='predict layer: least p score selected, uncalibrated '
                        '(default: %(default)s)')
    p.add_argument('--predict_scanner_flag', type=int, default=0,
                   choices=[0, 1],
                   help='predict layer: manufacturer dummies in the model '
                        '(default: %(default)s)')
    p.add_argument('--predict_max_frac_nan', type=float, default=0.5,
                   help='predict layer: predictors missing for more than this '
                        'fraction of the training sessions are left out '
                        '(default: %(default)s)')
    p.add_argument('--predict_min_n', type=int, default=50,
                   help='predict layer: least training sessions with a value '
                        'for a metric to be used (default: %(default)s)')
    p.add_argument('--predict_min_nfail', type=int, default=5,
                   help='predict layer: least failures among them '
                        '(default: %(default)s)')
    p.add_argument('--predict_n_draws', type=int, default=10,
                   help='predict layer: random column draws of the '
                        'Mahalanobis references (default: %(default)s)')
    p.add_argument('--random_fraction', type=float, default=0.25,
                   help='random sample size per modality, as a fraction of '
                        'the previous-failure rows (default: %(default)s)')
    p.add_argument('--random_floor', type=int, default=200,
                   help='least random sample size per modality '
                        '(default: %(default)s)')
    p.add_argument('--random_seed', type=int, default=int(today),
                   help='seed of the random sample (default: today, YYMMDD)')
    p.add_argument('--random_sample_file',
                   help='tsv of a random sample drawn before (participant_id, '
                        'session_id, modality, seed), used instead of a new '
                        'draw')
    p.add_argument('--exclude_sampled_files', nargs='+',
                   help='selection files (rows with reason random) or random '
                        'sample files of earlier rounds, whose sessions are '
                        'not drawn again')
    p.add_argument('--filter', default='new', choices=FILTERS,
                   help='rows kept by their reviews in --review_rootdir, as '
                        'autoqc_select_all: new not fully reviewed, old fully '
                        'reviewed, any every row, pass / fail fully reviewed '
                        'with consensus _score 1 / 0 (default: %(default)s)')
    p.add_argument('--max_rows', type=int,
                   help='most rows selected per modality, of those --filter '
                        'keeps (default: no cap)')
    p.add_argument('--prev_selected_file',
                   help='a previous selection file: with --max_rows, its rows '
                        'still selected and kept by --filter keep their place '
                        '(default: none)')
    p.add_argument('--max_nrev', type=int, default=2,
                   help='reviewers per session at which a row is fully '
                        'reviewed, as cmig_manual_postqc_bids --max_nrev '
                        '(default: %(default)s); 0 = no limit, no row fully '
                        'reviewed')
    p.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                   help='overwrite existing outputs (default: 0)')
    p.add_argument('--verbose', type=int, default=1, choices=[0, 1],
                   help='print progress and the summary (default: 1)')
    return p.parse_args(argv)


def main():
    errcode, errmsg = run(parse_args())
    if errcode:
        print('{0}: ERROR: {1}'.format(APPNAME, errmsg))
    sys.exit(errcode)


if __name__ == '__main__':
    main()
