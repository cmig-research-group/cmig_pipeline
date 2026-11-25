#!/usr/bin/env python3
"""Mahalanobis distances for post QC selection (port of MMPS autoqc_mahal.m).

The reference (column draws, scanner residualization betas, imputation means,
mean and covariance factor per draw, median scaling of the median over draws
of the log distance) is fitted on a reference set of rows by fit_reference
and applied unchanged to any rows by apply_reference.

parity_mahal reproduces autoqc_mahal.m as called by autoqc_mahal_allrois.m
(scanner_flag 1, combined_flag 1, excl_excess_nan_flag 1, normdist_flag 1):
the reference is all rows, and the scanner type is the first letter of
VisitID. read_roi_summary reads an old ROI summary csv as MATLAB readtable
does (column names made valid identifiers before the patterns are applied).

FAMILIES_80 maps the 11 MMPS families, and a task corr family per task, to
the release 8.0 compiled tables (cmig_compile_bids_tabdat.py); load_family_80 reads one family's columns,
and run_families_80 (the command line) fits the references on a set of
sessions (a list of reference sessions, or every session with values less an
exclude list), saves them, applies them to every session with values and
writes the distances (<family>_mdist, <family>_nanfrac). With scanner_flag 1
the references are residualized for the scanner type, the manufacturer in
mr_y_adm__info; with 0 they are not.

Python 3.6.9 compatible; numpy only.

Created:  11/19/2026 by Don Hagler
Prev mod: 11/21/2026 by Don Hagler
Last mod: 11/22/2026 by Don Hagler
"""

import argparse
import csv
import json
import math
import os
import re
import sys

import numpy as np

REFERENCE_VERSION = 3

# old families: (ROI summary file, column patterns), from autoqc_mahal_allrois.m
OLD_FAMILIES = [
  ('fsurf', 'MRI_all.csv', ['thick', 'area']),
  ('sMRI_T1w', 'MRI_all.csv', ['T1w']),
  ('sMRI_T2w', 'MRI_all.csv', ['T2w']),
  ('dMRI_DTI', 'DTI_full_all.csv', ['FA', 'MD']),
  ('dMRI_RSI', 'RSI_wmbc_all.csv', ['F0', 'FD', 'If']),
  ('rsfMRI_corr', 'rsBOLD_thresh2_network_corr_all.csv',
    ['DorsalAttn__', 'Default__']),
  ('rsfMRI_subcort', 'rsBOLD_thresh2_network_subcort_all.csv', ['__']),
  ('rsfMRI_var', 'rsBOLD_thresh2_var_all.csv', ['var']),
  ('tfMRI_MID', 'taskBOLD_MID_avg_beta_aseg_aparc2_all.csv', ['beta']),
  ('tfMRI_SST', 'taskBOLD_SST_avg_beta_aseg_aparc2_all.csv', ['beta']),
  ('tfMRI_nBack', 'taskBOLD_nBack_avg_beta_aseg_aparc2_all.csv', ['beta']),
]

# strings that autoqc_check_vals treats as NaN
NAN_STRINGS = set(['', 'nan', 'na', 'n/a', '_nan_', 'inf', '_inf_'])

NAMELENGTHMAX = 63


###############################################################################
# column names and patterns

def make_valid_names(names):
  """Column names as MATLAB readtable makes them (matlab.lang.makeValidName,
  then makeUniqueStrings)."""
  out = []
  for name in names:
    s = name.strip()
    # whitespace removed, the following letter capitalized
    s = re.sub(r'\s+([a-z])', lambda m: m.group(1).upper(), s)
    s = re.sub(r'\s+', '', s)
    s = re.sub(r'[^A-Za-z0-9_]', '_', s)
    if not s or not s[0].isalpha():
      s = 'x' + s
    out.append(s[:NAMELENGTHMAX])
  # later duplicates get _1, _2, ...
  seen = set(out)
  counts = {}
  first = set()
  for i, s in enumerate(out):
    if s not in first:
      first.add(s)
      continue
    k = counts.get(s, 0)
    while True:
      k += 1
      suffix = '_%d' % k
      t = s[:NAMELENGTHMAX - len(suffix)] + suffix
      if t not in seen:
        break
    counts[s] = k
    seen.add(t)
    out[i] = t
  return out


def match_patterns(names, patterns):
  """Sorted union of the indices of names matching any regexp pattern
  (as autoqc_mahal.m: regexp, then union); error if a pattern matches none."""
  if isinstance(patterns, str):
    patterns = [patterns]
  idx = set()
  counts = []
  for pat in patterns:
    rx = re.compile(pat)
    m = [i for i, n in enumerate(names) if rx.search(n)]
    if not m:
      raise ValueError('no columns with pattern %s' % pat)
    counts.append((pat, len(m)))
    idx.update(m)
  return sorted(idx), counts


###############################################################################
# reading an old ROI summary csv

def read_roi_summary(fname, patterns, id_col='VisitID', verbose=False):
  """Read the id column and the columns matching patterns from a csv.

  Returns a dict: ids (list), names (valid names of all matched columns, in
  file order), data (nrows x nmatched float array), isnum (bool per matched
  column: False for a column with a token that is not a number), ninf
  (count of Inf values per matched column), pattern_counts.
  """
  with open(fname, 'r', newline='') as f:
    reader = csv.reader(f)
    header = next(reader)
    names = make_valid_names(header)
    if id_col not in names:
      raise ValueError('%s not found in %s' % (id_col, fname))
    j_id = names.index(id_col)
    idx, pattern_counts = match_patterns(names, patterns)
    if verbose:
      for pat, n in pattern_counts:
        print('selected %d/%d columns with pattern %s' % (n, len(names), pat))
      sys.stdout.flush()
    ncol = len(idx)
    ids = []
    rows = []
    bad = np.zeros(ncol, dtype=bool)
    for row in reader:
      if not row:
        continue
      ids.append(row[j_id].strip())
      vals = [row[j] for j in idx]
      try:
        v = np.array(vals, dtype=float)
      except ValueError:
        v = np.empty(ncol)
        for k, s in enumerate(vals):
          try:
            v[k] = float(s)
          except ValueError:
            v[k] = np.nan
            if s.strip().lower() not in NAN_STRINGS:
              bad[k] = True
      rows.append(v)
  data = np.vstack(rows) if rows else np.zeros((0, ncol))
  ninf = np.sum(np.isinf(data), axis=0)
  return {
    'ids': ids,
    'names': [names[j] for j in idx],
    'data': data,
    'isnum': ~bad,
    'ninf': ninf,
    'pattern_counts': pattern_counts,
  }


###############################################################################
# helpers matching MATLAB numerics

def _matlab_rank(A):
  """rank as MATLAB: singular values above max(size(A))*eps(max(s))."""
  if A.size == 0:
    return 0
  s = np.linalg.svd(A, compute_uv=False)
  tol = max(A.shape) * np.spacing(np.max(s))
  return int(np.sum(s > tol))


def _nanstd(D):
  """nanstd (n-1 normalization) per column; NaN where fewer than 2 values."""
  n = np.sum(~np.isnan(D), axis=0)
  m = _nanmean_cols(D)
  ss = np.nansum((D - m) ** 2, axis=0)
  with np.errstate(invalid='ignore', divide='ignore'):
    sd = np.sqrt(ss / (n - 1))
  sd[n < 2] = np.nan
  return sd


def _nanmean_cols(D):
  with np.errstate(invalid='ignore'):
    n = np.sum(~np.isnan(D), axis=0)
    s = np.nansum(D, axis=0)
    m = s / np.maximum(n, 1)
  m[n == 0] = np.nan
  return m


def _mahal_factor(Dr):
  """Mean and scaled R factor, as MATLAB mahal: qr of the centered data;
  distances are sum((R'\\(y-mu)').^2)*(n-1)."""
  n, p = Dr.shape
  if n <= p:
    raise ValueError('mahal: the number of rows (%d) must exceed the number '
                     'of columns (%d)' % (n, p))
  mu = np.mean(Dr, axis=0)
  R = np.linalg.qr(Dr - mu, mode='r')
  return mu, R / math.sqrt(n - 1)


def _mahal_apply(Y, mu, Rs):
  z = np.linalg.solve(Rs.T, (Y - mu).T)
  return np.sum(z * z, axis=0)


###############################################################################
# fit and apply

def _valid_columns(X, excess_nan_thresh):
  nan_ratio = np.mean(np.isnan(X), axis=0)
  sd = _nanstd(X)
  # NaN sd counts as nonzero, as std_vals~=0 in MATLAB
  return np.where((nan_ratio < excess_nan_thresh) & ~(sd == 0))[0]


def fit_reference(X, colnames, scanner=None, ref_mask=None,
                  isnum=None, excess_nan_thresh=0.5, normdist_flag=True,
                  drop_dependent=True, dependent_tol=1e-10,
                  subsample='random', subsample_seed=0, rows_per_col=10,
                  n_draws=10, return_distances=False):
  """Fit the Mahalanobis reference on the rows in ref_mask.

  X: nrows x ncols array (NaN for missing); colnames: names of X columns.
  scanner: per-row scanner labels (None: no residualization).
  ref_mask: bool per row (None: all rows).
  isnum: bool per column (False: excluded, as a non-numeric column).
  drop_dependent: drop a draw's columns that are linear combinations of
    earlier ones (|R(i,i)| <= dependent_tol * max|R(i,i)|); with exact
    duplicates (e.g. A->B and B->A correlations) the covariance is singular
    and the distances are rounding noise. False for parity.
  subsample: 'random' fits n_draws seeded random draws (one RandomState
    from subsample_seed) of min(rank-1, floor(nrows/rows_per_col)) columns
    each, nrows the reference rows used; the distance is the median over
    draws of the log distance. 'step' fits one draw of every
    ceil(ncols/(rank-1))-th column (parity).
  The median and MAD scaling (normdist_flag) is of the combined log
  distance of the reference rows.
  Returns the reference dict (JSON-serializable), and with
  return_distances also the distances of the reference rows.
  """
  X = np.asarray(X, dtype=float)
  nrows, ncols = X.shape
  if len(colnames) != ncols:
    raise ValueError('colnames does not match X')
  if ref_mask is None:
    ref_mask = np.ones(nrows, dtype=bool)
  ref_mask = np.asarray(ref_mask, dtype=bool)
  idx_num = np.arange(ncols) if isnum is None else np.where(isnum)[0]
  Xr = X[ref_mask][:, idx_num]
  scan_r = None if scanner is None else \
    [s for s, m in zip(scanner, ref_mask) if m]

  # columns with < 50% NaN and nonzero SD, on the reference rows
  jv = _valid_columns(Xr, excess_nan_thresh)
  cols = [colnames[idx_num[j]] for j in jv]
  D = Xr[:, jv]

  # reference rows with excess NaNs are left out
  nanfrac = np.mean(np.isnan(D), axis=1) if D.shape[1] else \
    np.ones(D.shape[0])
  keep = nanfrac <= excess_nan_thresh
  D = D[keep]
  nanfrac = nanfrac[keep]

  ref = {
    'version': REFERENCE_VERSION,
    'excess_nan_thresh': excess_nan_thresh,
    'normdist_flag': bool(normdist_flag),
    'n_ref': int(np.sum(ref_mask)),
    'n_ref_used': int(D.shape[0]),
    'n_cols_input': int(ncols),
    'n_cols_numeric': int(len(idx_num)),
    'colnames': cols,
    'scanner_levels': None,
    'scanner_ref_level': None,
    'scanner_betas': None,
  }

  # residualize for scanner: per-level betas from rows without NaNs,
  #   offset of the most common level added back
  if scanner is not None:
    levels = sorted(set(scan_r))
    sv = np.array([levels.index(s) for s in scan_r])[keep]
    Xd = np.zeros((len(sv), len(levels)))
    Xd[np.arange(len(sv)), sv] = 1
    nonan = nanfrac == 0
    B = np.dot(np.linalg.pinv(Xd[nonan]), D[nonan])
    i_most = int(np.argmax(np.bincount(sv, minlength=len(levels))))
    D = D - np.dot(Xd, B) + B[i_most]
    ref['scanner_levels'] = levels
    ref['scanner_ref_level'] = levels[i_most]
    ref['scanner_betas'] = B.tolist()

  # impute NaNs with column means (0 for a column with no values)
  cmean = _nanmean_cols(D)
  cmean[np.isnan(cmean)] = 0
  D = np.where(np.isnan(D), cmean, D)
  ref['col_means'] = cmean.tolist()

  # rank once (the slow step), shared by all draws
  mmax = _matlab_rank(D) - 1
  if mmax < 1:
    raise ValueError('rank of the reference data is too low')
  ref['rank'] = int(mmax + 1)
  ref['subsample'] = subsample
  ref['subsample_step'] = None
  ref['subsample_seed'] = None
  ref['subsample_rows_per_col'] = None
  if subsample == 'step':
    # every ceil(ncol/(rank-1))-th column, one draw
    step = int(math.ceil(D.shape[1] / float(mmax)))
    subs = [np.arange(0, D.shape[1], step)]
    ref['subsample_step'] = step
  elif subsample == 'random':
    # ~rows_per_col reference rows per column, at most rank-1 columns
    m = min(mmax, D.shape[0] // rows_per_col, D.shape[1])
    if m < 1:
      raise ValueError('too few reference rows (%d) for %d rows per column' %
                       (D.shape[0], rows_per_col))
    if n_draws < 1:
      raise ValueError('n_draws must be at least 1')
    rng = np.random.RandomState(subsample_seed)
    subs = [np.sort(rng.choice(D.shape[1], m, replace=False))
            for _ in range(n_draws)]
    ref['subsample_seed'] = int(subsample_seed)
    ref['subsample_rows_per_col'] = int(rows_per_col)
  else:
    raise ValueError('subsample must be random or step')
  ref['n_draws'] = len(subs)
  ref['subsample_n_drawn'] = int(len(subs[0]))

  # per draw: column indices (into colnames), dependent columns dropped,
  #   mean, R upper triangle row by row
  ref['draws'] = []
  lds = []
  for sub in subs:
    mu, Rs = _mahal_factor(D[:, sub])
    dropped = []
    if drop_dependent:
      dg = np.abs(np.diag(Rs))
      dep = dg <= dependent_tol * np.max(dg)
      if np.any(dep):
        dropped = [cols[j] for j in sub[dep]]
        sub = sub[~dep]
        mu, Rs = _mahal_factor(D[:, sub])
    ref['draws'].append({
      'idx': sub.tolist(), 'dropped': dropped, 'mu': mu.tolist(),
      'R': [Rs[i, i:].tolist() for i in range(Rs.shape[0])]})
    lds.append(np.log(_mahal_apply(D[:, sub], mu, Rs)))
  ld = np.median(np.vstack(lds), axis=0)

  ref['median'] = 0.0
  ref['mad'] = 1.0
  if normdist_flag:
    med = float(np.median(ld))
    mad = float(np.median(np.abs(ld - med)))
    ref['median'] = med
    ref['mad'] = mad
    ld = (ld - med) / mad
  if return_distances:
    out = np.full(len(ref_mask), np.nan)
    i_ref = np.where(ref_mask)[0]
    out[i_ref[keep]] = ld
    return ref, out
  return ref


def _unpack_R(rows):
  p = len(rows)
  Rs = np.zeros((p, p))
  for i, r in enumerate(rows):
    Rs[i, i:] = r
  return Rs


def reference_draws(reference):
  """The draws of a reference; a version 1-2 reference (one subsample at the
  top level) as one draw."""
  if 'draws' in reference:
    return reference['draws']
  return [{'idx': reference['subsample_idx'],
           'dropped': reference.get('subsample_dropped', []),
           'mu': reference['mu'], 'R': reference['R']}]


def apply_reference(X, colnames, reference, scanner=None):
  """Distances of all rows of X from the reference.

  Returns (mdist, nanfrac): mdist NaN for rows with more than
  excess_nan_thresh NaNs in the reference columns; a scanner label not in
  the reference is not residualized. mdist is the median over draws of the
  log distance, scaled by the reference's median and MAD.
  """
  X = np.asarray(X, dtype=float)
  pos = dict((n, j) for j, n in enumerate(colnames))
  missing = [n for n in reference['colnames'] if n not in pos]
  if missing:
    raise ValueError('%d reference columns missing, e.g. %s' %
                     (len(missing), missing[0]))
  D = X[:, [pos[n] for n in reference['colnames']]]
  nrows = D.shape[0]
  nanfrac = np.mean(np.isnan(D), axis=1)
  keep = nanfrac <= reference['excess_nan_thresh']
  D = D[keep]

  if reference.get('scanner_levels') is not None:
    if scanner is None:
      raise ValueError('reference is residualized for scanner; '
                       'scanner labels required')
    levels = reference['scanner_levels']
    B = np.array(reference['scanner_betas'])
    i_ref = levels.index(reference['scanner_ref_level'])
    sc = [s for s, k in zip(scanner, keep) if k]
    iv = np.array([levels.index(s) if s in levels else i_ref for s in sc],
                  dtype=int)
    D = D - B[iv] + B[i_ref]

  cmean = np.array(reference['col_means'])
  D = np.where(np.isnan(D), cmean, D)
  lds = [np.log(_mahal_apply(D[:, dr['idx']], np.array(dr['mu']),
                             _unpack_R(dr['R'])))
         for dr in reference_draws(reference)]
  ld = np.median(np.vstack(lds), axis=0)
  ld = (ld - reference['median']) / reference['mad']
  mdist = np.full(nrows, np.nan)
  mdist[keep] = ld
  return mdist, nanfrac


def save_reference(reference, fname):
  with open(fname, 'w') as f:
    json.dump(reference, f)


def load_reference(fname):
  with open(fname, 'r') as f:
    return json.load(f)


###############################################################################
# parity with autoqc_mahal.m

def parity_mahal(ids, names, data, isnum, patterns, return_reference=False):
  """autoqc_mahal.m on columns read by read_roi_summary: the patterns select
  among names, the reference is all rows, the scanner type is the first
  letter of the id. Returns (mdist, nanfrac, scanner_types, info)."""
  idx, pattern_counts = match_patterns(names, patterns)
  scanner = [v[0] for v in ids]
  sel_names = [names[j] for j in idx]
  X = data[:, idx]
  ref, mdist = fit_reference(X, sel_names, scanner=scanner,
                             isnum=np.asarray(isnum)[idx],
                             drop_dependent=False, subsample='step',
                             return_distances=True)
  # nanfrac of all rows over the valid columns, as written by MATLAB
  pos = dict((n, j) for j, n in enumerate(sel_names))
  D = X[:, [pos[n] for n in ref['colnames']]]
  nanfrac = np.mean(np.isnan(D), axis=1)
  info = {
    'pattern_counts': pattern_counts,
    'n_matched': len(idx),
    'n_numeric': ref['n_cols_numeric'],
    'n_valid': len(ref['colnames']),
    'rank': ref['rank'],
    'n_subsample': len(ref['draws'][0]['idx']),
    'n_rows': len(ids),
    'n_rows_used': ref['n_ref_used'],
  }
  if return_reference:
    info['reference'] = ref
  return mdist, nanfrac, scanner, info


###############################################################################
# release 8.0 families, from the compiled tables

IDS = ('participant_id', 'session_id')
ADM_TABLE = 'mr_y_adm__info'
MANUFACTURER_COL = 'mr_y_adm__info__dev_manufact'
UNKNOWN_SCANNER = 'unknown'

# 8.0 contrast code and the old contrast it pairs with
#   (taskBOLD_<task>_avg_beta_aseg_aparc2_all.csv), by table name
TASK_CONTRASTS_80 = {
  'mid': [
    ('alvn', 'antic_loss_vs_neutral'),
    ('alrvn', 'antic_large_reward_vs_neutral'),
    ('allvn', 'antic_large_loss_vs_neutral'),
    ('alvsl', 'antic_large_vs_small_loss'),
    ('alvsr', 'antic_large_vs_small_reward'),
    ('arvn', 'antic_reward_vs_neutral'),
    ('aslvn', 'antic_small_loss_vs_neutral'),
    ('asrvn', 'antic_small_reward_vs_neutral'),
    ('lpvnf', 'loss_pos_vs_neg_feedback'),
    ('rpvnf', 'reward_pos_vs_neg_feedback'),
  ],
  'sst': [
    ('cg', 'correct_go_vs_fixation'),
    ('asvcg', 'any_stop_vs_correct_go'),
    ('csvcg', 'correct_stop_vs_correct_go'),
    ('csvis', 'correct_stop_vs_incorrect_stop'),
    ('igvcg', 'incorrect_go_vs_correct_go'),
    ('igvis', 'incorrect_go_vs_incorrect_stop'),
    ('isvcg', 'incorrect_stop_vs_correct_go'),
  ],
  'nback': [
    ('0b', '0_back'),
    ('2b', '2_back'),
    ('2bv0b', '2_back_vs_0_back'),
    ('emo', 'emotion'),
    ('emovntf', 'emotion_vs_neutface'),
    ('fvplc', 'face_vs_place'),
    ('ngfvntf', 'negface_vs_neutface'),
    ('plc', 'place'),
    ('psfvntf', 'posface_vs_neutface'),
  ],
}


def _tables(prefix, parts):
  """prefix__a__b... for every combination of the lists in parts."""
  out = [prefix]
  for p in parts:
    out = ['%s__%s' % (t, s) for t in out for s in p]
  return out


def _family(name, modality, task, tables, pattern=r'.*', drop=()):
  return {'name': name, 'modality': modality, 'task': task,
          'tables': tables, 'pattern': pattern, 'drop': list(drop)}


def _tfmri_family(name, task):
  tables = ['mr_y_tfmri__%s__%s__%s' % (task, c, a)
            for c, _ in TASK_CONTRASTS_80[task] for a in ('aseg', 'dst')]
  return _family(name, 'fmri', task, tables, pattern=r'.*_beta')


_GPNET = ['aud', 'cio', 'cip', 'def', 'doa', 'frp', 'non', 'ret', 'sal',
          'smh', 'smm', 'vea', 'vis']


def _tfmri_corr_family(name, task):
  """The resid (after the task model) network correlations and variances
  of a task; the gpnet pairs b__a (b > a) dropped, as duplicates of a__b."""
  pre = 'mr_y_tfmri__%s' % task
  tables = (['%s__corr__gpnet' % pre, '%s__corr__gpnet__aseg' % pre] +
            _tables(pre + '__var', [['aseg', 'dsk', 'dst', 'gp']]))
  pattern = r'%s__(corr__gpnet(__aseg)?|var)__resid__.*' % pre
  drop = ['%s__corr__gpnet__resid__%s__%s_mean' % (pre, b, a)
          for a in _GPNET for b in _GPNET if b > a]
  return _family(name, 'fmri', task, tables, pattern=pattern, drop=drop)


_SURF = ['dsk', 'dst']
_RSI = ['ri', 'hi', 'rd', 'hd', 'fi', 'fni']

# family: post QC modality, task, tables (tabulated/<table>.tsv), pattern
#   (full match on the data columns), columns dropped
FAMILIES_80 = [
  _family('fsurf', 'fsurf', '',
          _tables('mr_y_smri', [['thk', 'area'], ['dsk', 'dst', 'fzy']])),
  _family('sMRI_T1w', 'fsurf', '',
          _tables('mr_y_smri__t1', [['gm', 'wm', 'gwc'],
                                    ['dsk', 'dst', 'fzy']]) +
          ['mr_y_smri__t1__aseg']),
  _family('sMRI_T2w', 't2', '',
          _tables('mr_y_smri__t2', [['gm', 'wm', 'gwc'],
                                    ['dsk', 'dst', 'fzy']]) +
          ['mr_y_smri__t2__aseg']),
  _family('dMRI_DTI', 'dmri', '',
          _tables('mr_y_dti__fs', [['fa', 'md'], ['at', 'aseg']]) +
          _tables('mr_y_dti__fs', [['fa', 'md'], ['gm', 'wm', 'gwc'],
                                   _SURF])),
  _family('dMRI_RSI', 'dmri', '',
          _tables('mr_y_rsi', [_RSI, ['at', 'aseg']]) +
          _tables('mr_y_rsi', [_RSI, ['gm', 'wm', 'gwc'], _SURF])),
  # the two old networks; doa__def is def__doa, so one of the pair is dropped
  _family('rsfMRI_corr', 'fmri', 'rest', ['mr_y_rsfmri__corr__gpnet'],
          pattern=r'mr_y_rsfmri__corr__gpnet__(doa|def)__.*',
          drop=['mr_y_rsfmri__corr__gpnet__doa__def_mean']),
  _family('rsfMRI_subcort', 'fmri', 'rest',
          ['mr_y_rsfmri__corr__gpnet__aseg']),
  _family('rsfMRI_var', 'fmri', 'rest',
          _tables('mr_y_rsfmri__var', [['dsk', 'dst', 'gp', 'aseg']])),
  _tfmri_family('tfMRI_MID', 'mid'),
  _tfmri_family('tfMRI_SST', 'sst'),
  _tfmri_family('tfMRI_nBack', 'nback'),
  _tfmri_corr_family('tfMRI_MID_corr', 'mid'),
  _tfmri_corr_family('tfMRI_SST_corr', 'sst'),
  _tfmri_corr_family('tfMRI_nBack_corr', 'nback'),
]

FAMILY_NAMES_80 = [f['name'] for f in FAMILIES_80]


def _num(s):
  try:
    return float(s)
  except ValueError:
    return np.nan


def read_table_80(fname, pattern=r'.*', drop=()):
  """(keys, names, data) of a compiled tsv: keys (pid, sid) per row, the
  data columns fully matching pattern (less drop), data float (NaN for
  empty, n/a or not a number)."""
  rx = re.compile(pattern)
  drop = set(drop)
  with open(fname, 'r', newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    if any(c not in header for c in IDS):
      raise IOError('%s lacks participant_id or session_id' % fname)
    jp, js = header.index(IDS[0]), header.index(IDS[1])
    idx = [j for j, c in enumerate(header)
           if c not in IDS and c not in drop and rx.fullmatch(c)]
    keys, rows = [], []
    for row in reader:
      if not row:
        continue
      keys.append((row[jp], row[js]))
      rows.append([_num(row[j]) for j in idx])
  if len(set(keys)) != len(keys):
    raise IOError('%s: an id pair on more than one row' % fname)
  data = np.array(rows, dtype=float).reshape(len(rows), len(idx))
  return keys, [header[j] for j in idx], data


def load_family_80(tabulated_dir, family):
  """(keys, names, X) of one FAMILIES_80 entry (or name): the sessions with
  a value in any of the family's columns (sorted), the columns in table
  order, NaN where a session lacks a table or a value."""
  if isinstance(family, str):
    family = FAMILIES_80[FAMILY_NAMES_80.index(family)]
  parts = []
  for t in family['tables']:
    fname = os.path.join(tabulated_dir, t + '.tsv')
    if not os.path.isfile(fname):
      raise IOError('missing %s' % fname)
    keys, names, data = read_table_80(fname, family['pattern'],
                                      family['drop'])
    if not names:
      raise ValueError('%s: no columns match %s' % (fname, family['pattern']))
    parts.append((keys, names, data))
  allkeys = sorted(set(k for p in parts for k in p[0]))
  pos = dict((k, i) for i, k in enumerate(allkeys))
  names = [n for p in parts for n in p[1]]
  X = np.full((len(allkeys), len(names)), np.nan)
  j0 = 0
  for keys, nms, data in parts:
    X[[pos[k] for k in keys], j0:j0 + len(nms)] = data
    j0 += len(nms)
  has = np.any(~np.isnan(X), axis=1)
  return [k for k, h in zip(allkeys, has) if h], names, X[has]


def read_manufacturer(tabulated_dir):
  """{(pid, sid): manufacturer} from mr_y_adm__info ('' if n/a)."""
  fname = os.path.join(tabulated_dir, ADM_TABLE + '.tsv')
  out = {}
  with open(fname, 'r', newline='') as f:
    for r in csv.DictReader(f, delimiter='\t'):
      v = (r.get(MANUFACTURER_COL) or '').strip()
      out[(r[IDS[0]], r[IDS[1]])] = '' if v.lower() in NAN_STRINGS else v
  return out


def read_ids(fname):
  """(pid, sid) pairs of a tsv with participant_id and session_id."""
  with open(fname, 'r', newline='') as f:
    return [(r[IDS[0]], r[IDS[1]])
            for r in csv.DictReader(f, delimiter='\t')]


def read_exclude(fname):
  """(pid, sid, modality) of a tsv with participant_id, session_id and an
  optional modality column (a post QC modality: the row excludes the session
  from the references of that modality's families; '' or no column: from
  every family)."""
  with open(fname, 'r', newline='') as f:
    return [(r[IDS[0]], r[IDS[1]], (r.get('modality') or '').strip())
            for r in csv.DictReader(f, delimiter='\t')]


def reference_fname(refdir, name):
  return os.path.join(refdir, 'mahal_ref_%s.json' % name)


def run_families_80(tabulated_dir, outdir, ref_ids=None, refdir=None,
                    ref_exclude=None, ref_source=None, scanner_flag=1,
                    subsample_seed=0, rows_per_col=10, n_draws=10,
                    families=None, outfile='mahal_mdist.tsv', verbose=True):
  """Distances of the release 8.0 families.

  With ref_ids, each family's reference is fit on the ref_ids sessions with
  the family's values; with ref_exclude ((pid, sid, modality) as
  read_exclude), on every session with the family's values less those
  excluded for the family's modality. A fitted reference is saved to outdir
  (mahal_ref_<family>.json) with scanner_flag, ref_definition, ref_source
  (e.g. the file name), ref_ids and ref_excluded_ids. With refdir the
  references are loaded instead. scanner_flag 1 residualizes the fit for
  manufacturer, 0 does not; subsample_seed, rows_per_col and n_draws set the
  random column draws (fit_reference); the scanner labels are passed to the apply only
  for a residualized reference. The references are applied to every
  session with values; outdir/outfile gets participant_id, session_id,
  <family>_mdist and <family>_nanfrac (n/a without values). A family whose
  fit or apply fails gets n/a, its error in info. Returns {family: info}.
  """
  if sum(x is not None for x in (ref_ids, ref_exclude, refdir)) != 1:
    raise ValueError('give one of ref_ids, ref_exclude and refdir')
  families = families or FAMILY_NAMES_80
  manuf = read_manufacturer(tabulated_dir)
  ref_set = None if ref_ids is None else set(tuple(k) for k in ref_ids)
  excl = None
  if ref_exclude is not None:
    excl = {}
    for pid, sid, mod in ref_exclude:
      excl.setdefault((pid, sid), set()).add(mod)
  if not os.path.isdir(outdir):
    os.makedirs(outdir)
  results, infos = {}, {}
  for name in families:
    fam = FAMILIES_80[FAMILY_NAMES_80.index(name)]
    keys, names, X = load_family_80(tabulated_dir, fam)
    scanner = [manuf.get(k) or UNKNOWN_SCANNER for k in keys]
    info = {'n_cols': len(names), 'n_sessions': len(keys),
            'n_no_manufacturer': sum(s == UNKNOWN_SCANNER for s in scanner),
            'error': None}
    try:
      if refdir is None:
        if ref_set is not None:
          mask = np.array([k in ref_set for k in keys], dtype=bool)
          ref_def = 'ref_ids: the listed sessions with values'
        else:
          mask = np.array([not (excl.get(k, set()) & set(['', fam['modality']]))
                           for k in keys], dtype=bool)
          ref_def = ('all minus excluded: the sessions with values, less '
                     'those excluded for modality %s' % fam['modality'])
        info['n_ref_no_manufacturer'] = int(sum(
          s == UNKNOWN_SCANNER for s, m in zip(scanner, mask) if m))
        ref = fit_reference(X, names,
                            scanner=scanner if scanner_flag else None,
                            ref_mask=mask, subsample_seed=subsample_seed,
                            rows_per_col=rows_per_col, n_draws=n_draws)
        ref['family'] = name
        ref['modality'] = fam['modality']
        ref['tables'] = fam['tables']
        ref['scanner_flag'] = int(bool(scanner_flag))
        ref['ref_definition'] = ref_def
        ref['ref_source'] = ref_source
        ref['ref_ids'] = [list(k) for k, m in zip(keys, mask) if m]
        ref['ref_excluded_ids'] = [list(k) for k, m in zip(keys, mask)
                                   if not m] if excl is not None else None
        save_reference(ref, reference_fname(outdir, name))
      else:
        ref = load_reference(reference_fname(refdir, name))
      residualized = ref.get('scanner_levels') is not None
      draws = reference_draws(ref)
      mdist, nanfrac = apply_reference(
        X, names, ref, scanner=scanner if residualized else None)
      results[name] = (keys, mdist, nanfrac)
      info.update({
        'n_ref': ref['n_ref'], 'n_ref_used': ref['n_ref_used'],
        'n_valid': len(ref['colnames']), 'rank': ref['rank'],
        'subsample': ref.get('subsample', 'step'),
        'subsample_seed': ref.get('subsample_seed'),
        'n_draws': len(draws),
        'n_kept': [len(d['idx']) for d in draws],
        'n_dropped_dependent': [len(d['dropped']) for d in draws],
        'median': ref['median'], 'mad': ref['mad'],
        'scanner_flag': int(residualized),
        'scanner_levels': ref['scanner_levels'],
        'n_scored': int(np.sum(~np.isnan(mdist))),
        'n_nan': int(np.sum(np.isnan(mdist))),
      })
    except (ValueError, np.linalg.LinAlgError) as e:
      info['error'] = str(e)
    infos[name] = info
    if verbose:
      print('%s: %s' % (name, json.dumps(info)))
      sys.stdout.flush()

  allkeys = sorted(set(k for keys, _, _ in results.values() for k in keys))
  pos = dict((k, i) for i, k in enumerate(allkeys))
  cols = []
  M = np.full((len(allkeys), 2 * len(families)), np.nan)
  for j, name in enumerate(families):
    cols += ['%s_mdist' % name, '%s_nanfrac' % name]
    if name in results:
      keys, mdist, nanfrac = results[name]
      i = [pos[k] for k in keys]
      M[i, 2 * j] = mdist
      M[i, 2 * j + 1] = nanfrac
  with open(os.path.join(outdir, outfile), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(list(IDS) + cols)
    for k, row in zip(allkeys, M):
      w.writerow(list(k) + ['n/a' if np.isnan(v) else '%.10g' % v
                            for v in row])
  return infos


def parse_args(argv=None):
  p = argparse.ArgumentParser(
    description='Mahalanobis distances of the release 8.0 post QC families')
  p.add_argument('--tabulated_dir', required=True,
                 help='the compiled tables (cmig_compile_bids_tabdat.py)')
  p.add_argument('--outdir', required=True,
                 help='output dir: distances tsv and fitted references')
  p.add_argument('--ref_file', default=None,
                 help='tsv (participant_id, session_id) of the reference '
                      'sessions; the references are fit and saved')
  p.add_argument('--ref_exclude_file', default=None,
                 help='tsv (participant_id, session_id, optional modality) '
                      'of the sessions left out of the references, which '
                      'are fit on every other session with values and '
                      'saved; a row with a post QC modality excludes the '
                      'session from that modality\'s families only')
  p.add_argument('--refdir', default=None,
                 help='dir of saved references (mahal_ref_<family>.json), '
                      'applied without a fit')
  p.add_argument('--scanner_flag', type=int, default=1, choices=[0, 1],
                 help='residualize the fitted references for manufacturer '
                      '(default: %(default)s)')
  p.add_argument('--subsample_seed', type=int, default=0,
                 help='seed of the random column draws (default: %(default)s)')
  p.add_argument('--rows_per_col', type=int, default=10,
                 help='reference rows per kept column (default: %(default)s)')
  p.add_argument('--n_draws', type=int, default=10,
                 help='random column draws; the distance is the median over '
                      'draws of the log distance (default: %(default)s)')
  p.add_argument('--families', nargs='+', default=None,
                 choices=FAMILY_NAMES_80, help='default all 11')
  p.add_argument('--outfile', default='mahal_mdist.tsv')
  p.add_argument('--verbose', type=int, default=1)
  return p.parse_args(argv)


def main(argv=None):
  args = parse_args(argv)
  ref_ids = read_ids(args.ref_file) if args.ref_file else None
  ref_exclude = read_exclude(args.ref_exclude_file) \
    if args.ref_exclude_file else None
  ref_source = args.ref_file or args.ref_exclude_file
  infos = run_families_80(args.tabulated_dir, args.outdir, ref_ids=ref_ids,
                          refdir=args.refdir, ref_exclude=ref_exclude,
                          ref_source=ref_source,
                          scanner_flag=args.scanner_flag,
                          subsample_seed=args.subsample_seed,
                          rows_per_col=args.rows_per_col,
                          n_draws=args.n_draws,
                          families=args.families, outfile=args.outfile,
                          verbose=bool(args.verbose))
  return int(any(i['error'] for i in infos.values()))


if __name__ == '__main__':
  sys.exit(main())
