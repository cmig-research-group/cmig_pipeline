#!/usr/bin/env python3
"""
Predictors and labels for the post QC prediction layer
(cmig_select_manual_postqc_bids layer_predict).

load_predictors(tabulated_dir, modality, predictor_file, ...)
    the predictor matrix of one modality from the compiled tables
    (cmig_compile_bids_tabdat.py) in tabulated_dir, with the columns the
    predictor file (parms/autoqc_predictors.tsv) lists for the modality.
    Per-task fMRI columns stay per task; a session missing a task (all of the
    task's columns empty) gets a 0/1 predictor fmri__<task>__missing_indicator
    and the task's columns filled: a raw QC column with the median of the
    same statistic over the session's other tasks, any other column (and a
    raw QC column with no other task's value) with the median of the
    reference sessions that have the task. Those medians are returned, so a
    saved model can apply the same fill to new sessions (fill=).

postqc_ids(tabulated_dir, modality, predictor_file)
    the sessions with the modality's post QC done: a value in any of its
    post QC (mr_y_qc__post__aut) columns.

screen_missing(X, rows, max_frac_nan)
    the columns missing for at most max_frac_nan of the given rows.

load_labels(rel_dir, modality)
    the selecting label of one modality: the overall QC of the release 7.0
    manual post QC table (REL/mr_y_qc__post__man__<mod>.tsv), pass 1, fail
    0, NaN where not reviewed (LABEL_SPEC gives the coding).

Values that are empty, n/a or not numbers are read as NaN. Python 3.6.9;
numpy and csv only.

Created:  11/19/2026 by Don Hagler
Prev mod: 11/19/2026 by Don Hagler
Last mod: 11/22/2026 by Don Hagler
"""

import csv
import os
import re

import numpy as np

MODALITIES = ['fsurf', 't2', 'dmri', 'fmri']
TASKS = ['rest', 'mid', 'sst', 'nback']
TYPES = ['exact', 'pattern', 'mahal']
IDS = ('participant_id', 'session_id')
POSTQC_TABLE = 'mr_y_qc__post__aut'
RAWQC_PREFIX = 'mr_y_qc__raw__'

# the selecting label per modality: column, value of a pass and of a fail
LABEL_PREFIX = 'mr_y_qc__post__man__'
LABEL_SPEC = dict((mod, {'column': LABEL_PREFIX + mod + '_score',
                         'type': 'binary', 'pass': 1, 'fail': 0})
                  for mod in MODALITIES)


def missing_name(task):
    """The name of a task's 0/1 missing-task predictor."""
    return 'fmri__{0}__missing_indicator'.format(task)


def num(val):
    """A table value as a float, NaN if empty, n/a or not a number."""
    try:
        return float((val or '').strip())
    except ValueError:
        return float('nan')


def read_predictor_file(fname, modality=None):
    """The rows (dicts: modality, table, type, column, task, note) of the
    predictor file, '#' lines skipped; those of one modality if given."""
    with open(fname, newline='') as f:
        lines = [l for l in f if not l.startswith('#')]
    rows = [dict(r) for r in csv.DictReader(lines, delimiter='\t')]
    for r in rows:
        r['task'] = (r.get('task') or '').strip()
        if r['modality'] not in MODALITIES:
            raise ValueError('{0}: unknown modality {1}'.format(
                fname, r['modality']))
        if r['type'] not in TYPES:
            raise ValueError('{0}: unknown type {1}'.format(fname, r['type']))
        if r['task'] and r['task'] not in TASKS:
            raise ValueError('{0}: unknown task {1}'.format(fname, r['task']))
    if modality:
        rows = [r for r in rows if r['modality'] == modality]
    return rows


def read_table(fname):
    """(columns, {(pid, sid): row dict}) of a tsv with participant_id and
    session_id; an id pair on more than one row is an error."""
    with open(fname, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        cols = list(reader.fieldnames or [])
        if any(c not in cols for c in IDS):
            raise IOError('{0} lacks participant_id or session_id'.format(fname))
        data = {}
        for r in reader:
            key = (r['participant_id'], r['session_id'])
            if key in data:
                raise IOError('{0}: more than one row for {1} {2}'.format(
                    fname, key[0], key[1]))
            data[key] = r
    return cols, data


def resolve_columns(prow, cols):
    """The columns of a table a predictor row names (exact: the column;
    pattern, mahal: the full matches of the regular expression)."""
    if prow['type'] == 'exact':
        return [prow['column']] if prow['column'] in cols else []
    pat = re.compile(prow['column'])
    return [c for c in cols if c not in IDS and pat.fullmatch(c)]


def load_predictors(tabulated_dir, modality, predictor_file, ids=None,
                    ref_ids=None, fill=None, mahal_file=None, strict=True):
    """(ids, X, colnames, info) for one modality.

    ids          the sessions (pid, sid) of the rows, in this order; default
                 every session with a value in any of the modality's
                 session (non-task) columns, sorted
    ref_ids      the sessions the task medians are computed on (those with
                 the task); default ids
    fill         {colname: value}: the task fill of a saved model, used in
                 place of medians of ref_ids
    mahal_file   tsv of the Mahalanobis distances (participant_id,
                 session_id, columns); without it the mahal rows are skipped
    strict       a listed column absent from its table is an error; else it
                 is all NaN and named in info['absent']

    X is a float array (sessions x columns), NaN where missing (task columns
    of a missing task are filled, and a task's Mahalanobis columns that are
    empty for a session with the task: the distance with the median, the
    NaN fraction with 1). info: fill {colname: value}, tasks {colname:
    task}, n_missing_task {task: n}, n_filled_other_task (values filled from
    the session's other tasks), n_mahal_gap {colname: n}, absent, skipped.
    """
    prows = read_predictor_file(predictor_file, modality)
    tables, sources = {}, {}
    colnames, coltask, absent, skipped = [], {}, [], []
    rawstat = {}
    for prow in prows:
        table = prow['table']
        if prow['type'] == 'mahal':
            if not mahal_file:
                skipped.append(prow['column'])
                continue
            fname = mahal_file
        else:
            fname = os.path.join(tabulated_dir, table + '.tsv')
        if fname not in tables:
            if not os.path.isfile(fname):
                raise IOError('missing {0}'.format(fname))
            tables[fname] = read_table(fname)
        cols = resolve_columns(prow, tables[fname][0])
        if not cols:
            if strict:
                raise IOError('{0}: no column {1}'.format(fname, prow['column']))
            absent.append(prow['column'])
            cols = [prow['column']]
        for c in cols:
            if c in sources:
                continue
            sources[c] = fname
            colnames.append(c)
            coltask[c] = prow['task']
            # a raw QC column's statistic, the same name in each task's table
            if table.startswith(RAWQC_PREFIX) and c.startswith(table + '__'):
                rawstat[c] = c[len(table) + 2:]

    def value(key, c):
        return num(tables[sources[c]][1].get(key, {}).get(c))

    session_cols = [c for c in colnames if not coltask[c]]
    if ids is None:
        keys = set()
        for fname in set(sources[c] for c in session_cols):
            keys.update(tables[fname][1].keys())
        ids = sorted(k for k in keys
                     if any(not np.isnan(value(k, c)) for c in session_cols))
    ids = [tuple(k) for k in ids]
    X = np.array([[value(k, c) for c in colnames] for k in ids],
                 dtype=float).reshape(len(ids), len(colnames))

    # per task: missing where every column of the task is NaN
    tasks = [t for t in TASKS if any(coltask[c] == t for c in colnames)]
    index = dict((k, i) for i, k in enumerate(ids))
    ref = [index[tuple(k)] for k in (ref_ids if ref_ids is not None else ids)
           if tuple(k) in index]
    out_fill, n_missing, extra = {}, {}, []
    miss = {}
    for t in tasks:
        jt = [j for j, c in enumerate(colnames) if coltask[c] == t]
        miss[t] = np.all(np.isnan(X[:, jt]), axis=1)
        n_missing[t] = int(miss[t].sum())
    X0 = X.copy()
    n_other = 0
    n_mahal_gap = {}
    for t in tasks:
        jt = [j for j, c in enumerate(colnames) if coltask[c] == t]
        for j in jt:
            c = colnames[j]
            if fill is not None and c in fill:
                v = float(fill[c])
            else:
                vals = X0[[i for i in ref if not miss[t][i]], j]
                vals = vals[~np.isnan(vals)]
                v = float(np.median(vals)) if vals.size else float('nan')
            out_fill[c] = v
            # same statistic in the other tasks' raw QC tables
            other = [jo for jo, co in enumerate(colnames)
                     if c in rawstat and coltask[co] not in ('', t) and
                     rawstat.get(co) == rawstat[c]]
            for i in np.nonzero(miss[t])[0]:
                vo = [X0[i, jo] for jo in other
                      if not miss[coltask[colnames[jo]]][i] and
                      not np.isnan(X0[i, jo])]
                if vo:
                    X[i, j] = float(np.median(vo))
                    n_other += 1
                else:
                    X[i, j] = v
            # a Mahalanobis family of a task the session has, but with no
            #   values of its own (e.g. the task corr, empty where no run
            #   kept enough frames): the distance takes the median, the NaN
            #   fraction 1
            if mahal_file and sources[c] == mahal_file:
                gap = np.nonzero(~miss[t] & np.isnan(X0[:, j]))[0]
                X[gap, j] = 1.0 if c.endswith('_nanfrac') else v
                n_mahal_gap[c] = len(gap)
        extra.append((missing_name(t), miss[t].astype(float)))
    if extra:
        X = np.column_stack([X] + [v for _, v in extra])
        for name, _ in extra:
            colnames.append(name)
            coltask[name] = name.split('__')[1]
    info = {'fill': out_fill, 'tasks': dict((c, t) for c, t in coltask.items()
                                            if t),
            'n_missing_task': n_missing, 'n_filled_other_task': n_other,
            'n_mahal_gap': n_mahal_gap,
            'absent': absent, 'skipped': skipped}
    return ids, X, colnames, info


def load_labels(rel_dir, modality):
    """{(pid, sid): label} from the release 7.0 manual post QC table of the
    modality: 1.0 pass, 0.0 fail, NaN not reviewed (LABEL_SPEC[modality])."""
    spec = LABEL_SPEC[modality]
    fname = os.path.join(rel_dir, LABEL_PREFIX + modality + '.tsv')
    if not os.path.isfile(fname):
        raise IOError('missing {0}'.format(fname))
    cols, data = read_table(fname)
    if spec['column'] not in cols:
        raise IOError('{0} lacks column {1}'.format(fname, spec['column']))
    labels = {}
    for key, r in data.items():
        v = num(r.get(spec['column']))
        if not np.isnan(v) and v not in (spec['pass'], spec['fail']):
            raise ValueError('{0}: {1} {2}: label {3} not {4} or {5}'.format(
                fname, key[0], key[1], v, spec['pass'], spec['fail']))
        labels[key] = v
    return labels


def postqc_ids(tabulated_dir, modality, predictor_file):
    """The sessions (pid, sid), sorted, with a value in any of the
    modality's post QC (POSTQC_TABLE) columns of the predictor file."""
    prows = [r for r in read_predictor_file(predictor_file, modality)
             if r['table'] == POSTQC_TABLE]
    fname = os.path.join(tabulated_dir, POSTQC_TABLE + '.tsv')
    if not os.path.isfile(fname):
        raise IOError('missing {0}'.format(fname))
    cols, data = read_table(fname)
    use = []
    for prow in prows:
        use += [c for c in resolve_columns(prow, cols) if c not in use]
    return sorted(k for k, r in data.items()
                  if any(not np.isnan(num(r.get(c))) for c in use))


def screen_missing(X, rows, max_frac_nan=0.5):
    """Indices of the columns of X missing (NaN) for at most max_frac_nan of
    the given rows (all columns if there are no rows)."""
    X = np.asarray(X, dtype=float)
    rows = np.asarray(rows, dtype=int)
    if rows.size == 0:
        return list(range(X.shape[1]))
    frac = np.mean(np.isnan(X[rows]), axis=0)
    return [j for j in range(X.shape[1]) if frac[j] <= max_frac_nan]
