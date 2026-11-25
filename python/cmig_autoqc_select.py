#!/usr/bin/env python3
"""Multivariate selection of sessions for manual post QC review.

Port of the multivariate path (multivar_flag 1) of MMPS autoqc_select_all.m and
the functions it calls (autoqc_calc_xcorr, autoqc_calc_xauc, autoqc_classify_score,
autoqc_calc_roc, autoqc_rocp, autoqc_train/test_multivar_prediction,
autoqc_select_cases with calibration, mmil_wtd_mad).  numpy only.
Changed from the MATLAB: no histogram gate in classify_score, NaN x keeps a NaN
score, weighted median/MAD over non-missing values, labels outside idx_train
masked, a separated fit rejected as 'separated'.

Library use:
    res = train(X, xnames, Y, labels, ids=ids, parms=parms)
        screen metrics per label (rthr, athr, min_n, min_nfail), z-score, SVD,
        logistic fit; res['model'] is a dict for <mod>_model.json, res['pscore']
        holds each label's score for every row, res['metrics_rows'] and
        res['labels_rows'] the screening and per-label records
    pscore = score(model, X, xnames)
        apply a saved model to new rows
    sel = select(pscore, Y, model, parms, idx_allowed=...)
        per-label pthr/nselect, weighting, one row per case, nall cap

Labels are typed: {'name': ..., 'type': 'pass1'} for a binary QC with pass = 1,
fail = 0, or {'name': ..., 'type': 'score'} for a QU score that fails above ythr.
label_mode 'parity' instead applies the MATLAB rules by name (a 0/1 label not
named QU_* is inverted; the AUC screen inverts any exact {0,1} label).

Not ported: univariate path, plots, nrand, scanner_flag 2-3, qcnew and
train_new_flag.
"""

import json
import math
import datetime

import numpy as np

EPS = np.finfo(float).eps

DEFAULT_PARMS = {
    'rthr': 0.2,
    'max_auc_flag': False,
    'athr': 0.5,
    'min_n': 100,
    'min_nfail': 2,
    'max_frac_nan': 0.1,
    'max_accept_frac': 0.5,
    'max_fapf': math.inf,
    'min_boost': 0.0,
    'rthr_margin': 0.0,
    'max_marginal_fapf': None,
    'max_shift_ratio': math.inf,
    'max_nvar': 10,
    'sthr': 0.05,
    'monoflag': True,
    'scanner_flag': 0,
    'nselect': math.inf,
    'nall': math.inf,
    'weighted_flag': True,
    'filter': 'new',
    'pthr': 0.0,
    'ythr': 2.0,
    'alpha': 0.05,
    'label_mode': 'typed',
}

SCREEN_COLS = ['yname', 'xname', 'n', 'n_fail', 'n_pass', 'frac_nan', 'R', 'AUC',
               'cohens_d', 'selected', 'rthr', 'athr', 'min_n', 'min_nfail',
               'rthr_margin']

LABEL_COLS = ['yname', 'xname', 'status', 'n_metrics', 'n', 'n_fail', 'base_rate',
              'n_sel', 'n_notsel', 'precision', 'boost', 'false_omission_rate',
              'fa_per_fail', 'false_discovery_rate', 'recall', 'roc_auc', 'pos_rate',
              'weight', 'max_abs_R', 'near_rthr', 'calib_thr', 'ythr', 'max_fapf',
              'min_boost', 'rthr', 'rthr_margin', 'max_marginal_fapf',
              'n_sel_fail', 'n_new', 'n_new_sel', 'shift_ratio',
              # not in the MATLAB csv
              'n_unscored', 'glm_iter', 'iter_limit', 'separated', 'score_method',
              'accept_frac']


def check_parms(parms=None):
    p = dict(DEFAULT_PARMS)
    if parms:
        unknown = set(parms) - set(p)
        if unknown:
            raise ValueError('unknown parms: %s' % ', '.join(sorted(unknown)))
        p.update(parms)
    if p['scanner_flag'] not in (0, 1):
        raise ValueError('scanner_flag must be 0 or 1')
    if p['filter'] not in ('new', 'old', 'any', 'pass', 'fail'):
        raise ValueError('unrecognized filter %s' % p['filter'])
    if p['label_mode'] not in ('typed', 'parity'):
        raise ValueError('label_mode must be typed or parity')
    return p


###########################################################################
# MATLAB numerics

def round_half_away(x):
    """MATLAB round: half away from zero (python round is half to even);
    +-Inf and NaN as they are."""
    if math.isinf(x) or math.isnan(x):
        return x
    return math.copysign(math.floor(abs(x) + 0.5), x)


def nanmean(x):
    x = np.asarray(x, dtype=float)
    v = x[~np.isnan(x)]
    return float(np.mean(v)) if v.size else math.nan


def nanstd(x):
    """MATLAB nanstd: N-1 normalization, 0 for one value, NaN for none."""
    x = np.asarray(x, dtype=float)
    v = x[~np.isnan(x)]
    if v.size == 0:
        return math.nan
    if v.size == 1:
        return 0.0
    return float(np.std(v, ddof=1))


def normpdf(x, mu, sigma):
    if not (sigma > 0):
        sigma = math.nan
    with np.errstate(invalid='ignore', over='ignore', under='ignore'):
        return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (math.sqrt(2 * math.pi) * sigma)


def colon_matlab(a, d, b):
    """a:d:b as MATLAB's colon (end point within tolerance, symmetric fill)."""
    if d == 0 or (d > 0 and a > b) or (d < 0 and a < b):
        return np.zeros(0)
    tol = 2.0 * EPS * max(abs(a), abs(b))
    sig = 1.0 if d > 0 else -1.0
    n = int(math.floor((b - a) / d))
    if sig * (a + (n + 1) * d - b) <= tol:
        n += 1
    c = a + n * d
    if sig * (c - b) > -tol:
        c = b
    v = np.zeros(n + 1)
    k = np.arange(0, n // 2 + 1)
    v[k] = a + k * d
    v[n - k] = c - k * d
    if n % 2 == 0:
        v[n // 2] = (a + c) / 2
    return v


def wtd_median(x):
    """mmil_wtd_median with unit weights, over the non-missing values."""
    x = np.asarray(x, dtype=float).ravel()
    x = x[~np.isnan(x)]
    n = x.size
    if n == 0:
        return math.nan
    sortx = np.sort(x, kind='mergesort')
    csumw = np.arange(1, n + 1, dtype=float)
    midpoint = n / 2.0
    j = np.nonzero(csumw <= midpoint)[0]
    if j.size == 0:
        return float(sortx[0])
    j = int(j[-1])
    if j == n - 1:
        return float(sortx[j])
    if csumw[j] == midpoint:
        return float((sortx[j] + sortx[j + 1]) / 2)
    return float(sortx[j + 1])


def wtd_mad(x):
    """mmil_wtd_mad with unit weights: (mad / 0.6745, median)."""
    med = wtd_median(x)
    return wtd_median(np.abs(np.asarray(x, dtype=float) - med)) / 0.6745, med


def corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xc = x - np.mean(x)
    yc = y - np.mean(y)
    with np.errstate(invalid='ignore', divide='ignore'):
        return float(np.sum(xc * yc) / math.sqrt(np.sum(xc * xc) * np.sum(yc * yc)))


def _stable_desc(v):
    """Indices sorting v descending, stable (MATLAB sort 'descend')."""
    return np.argsort(-np.asarray(v, dtype=float), kind='mergesort')


###########################################################################
# labels

def _label_spec(label):
    if isinstance(label, str):
        return {'name': label, 'type': None}
    spec = dict(label)
    spec.setdefault('type', None)
    return spec


def _parity_inverts(y, name):
    """check_ythr: 0/1 values and a name not starting with QU_."""
    v = np.unique(y[~np.isnan(y)])
    return bool(np.all(np.isin(v, [0, 1]))) and not name.startswith('QU_')


def _xauc_inverts(y):
    """autoqc_calc_xauc: values exactly {0, 1}."""
    v = np.unique(y[~np.isnan(y)])
    return v.size == 2 and v[0] == 0 and v[1] == 1


def prepare_label(y, spec, parms, idx_train):
    """Fail-oriented y and threshold for one label: (y, ythr, yname, inverted)."""
    y = np.asarray(y, dtype=float).copy()
    name = spec['name']
    if parms['label_mode'] == 'parity':
        invert = _parity_inverts(y, name)
    else:
        if spec['type'] not in ('pass1', 'score'):
            raise ValueError('label %s: type must be pass1 or score' % name)
        invert = spec['type'] == 'pass1'
    ythr = spec.get('ythr', parms['ythr'])
    if invert:
        ok = ~np.isnan(y)
        y[ok] = 1.0 * (y[ok] == 0)
        ythr = 0.5
        yname = '~' + name if parms['label_mode'] == 'parity' else name
    else:
        yname = name
    return y, ythr, yname, invert


def _auc_label(y_raw, spec, parms, y_fail, ythr_fail):
    """y and ythr the AUC screen uses (parity: autoqc_calc_xauc's own rule)."""
    if parms['label_mode'] != 'parity':
        return y_fail, ythr_fail
    y = np.asarray(y_raw, dtype=float).copy()
    if _xauc_inverts(y):
        ok = ~np.isnan(y)
        y[ok] = 1.0 * (y[ok] == 0)
        return y, 0.5
    return y, parms['ythr']


###########################################################################
# score, ROC

def classify_score(x, y, ythr=2.0, monoflag=False):
    """autoqc_classify_score; returns (pscore, params) for apply_score."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = x.size
    params = {'method': 'none'}
    idxo = np.nonzero(~np.isnan(x) & ~np.isnan(y))[0]
    with np.errstate(invalid='ignore'):
        if ythr > 0:
            pmask = y > abs(ythr)
            fmask = y <= abs(ythr)
        else:
            pmask = y <= abs(ythr)
            fmask = y > abs(ythr)
    x1 = x[pmask & ~np.isnan(x)]
    x2 = x[fmask & ~np.isnan(x)]

    # Gaussian log likelihood ratio, fail (1) vs pass (2): needs each group
    #   with >= 2 values and SD > 0
    pscore = np.full(n, math.nan)
    if x1.size >= 2 and x2.size >= 2:
        mu1 = nanmean(x1)
        sig1 = nanstd(x1)
        mu2 = nanmean(x2)
        sig2 = nanstd(x2)
        if sig1 > 0 and sig2 > 0:
            params.update({'method': 'gauss', 'mu1': mu1, 'sig1': sig1,
                           'mu2': mu2, 'sig2': sig2, 'monoflag': bool(monoflag)})
            with np.errstate(divide='ignore', invalid='ignore'):
                pscore = np.log(normpdf(x, mu1, sig1)) - np.log(normpdf(x, mu2, sig2))
            if monoflag:
                pscore, mono = _mono(x, pscore, mu1, sig1, mu2, sig2)
                params['mono'] = mono

    if not np.any(~np.isnan(pscore)):
        # robust z of x, signed by the correlation with y
        mad_x, med_x = wtd_mad(x)
        # MAD 0 (most values equal, e.g. a 0/1 indicator): scale by the SD;
        #   SD 0 too (constant): no score
        scale, scale_type = mad_x, 'mad'
        if not mad_x > 0:
            scale, scale_type = nanstd(x), 'sd'
        sgn = math.nan  # no labeled rows: no sign
        if idxo.size:
            r = corr(x[idxo], y[idxo])
            sgn = float(np.sign(r)) if not math.isnan(r) else math.nan
        if scale > 0:
            pscore = (x - med_x) / scale * sgn
        else:
            pscore = np.full(n, math.nan)
        params = {'method': 'fallback_mad', 'med': med_x, 'mad': mad_x,
                  'scale': scale, 'scale_type': scale_type, 'sign': sgn}
    return pscore, params


def _mono(x, pscore, mu1, sig1, mu2, sig2):
    """monoflag: hold the tails at the min / max, as autoqc_classify_score."""
    mono = {'applied': False}
    order = np.argsort(x, kind='mergesort')  # NaN last, as MATLAB sort
    ps = pscore[order]
    xs = x[order]
    ok = ~np.isnan(ps)
    if not ok.any():
        return pscore, mono
    idx_min = int(np.nanargmin(ps))
    idx_max = int(np.nanargmax(ps))
    minp = float(ps[idx_min])
    maxp = float(ps[idx_max])
    if sig1 > sig2:
        keep_high = mu1 > mu2
    else:
        keep_high = idx_max > idx_min
    nvalid = int(np.sum(~np.isnan(xs)))  # NaN x sorts last and stays NaN
    if keep_high:
        ps[:idx_min + 1] = minp
    else:
        ps[idx_min:nvalid] = minp
    hold_max = None
    if sig1 <= sig2:
        if mu1 > mu2:
            if idx_max < nvalid - 1:
                ps[idx_max:nvalid] = maxp
            hold_max = 'high'
        else:
            ps[:idx_max + 1] = maxp
            hold_max = 'low'
    out = np.empty_like(pscore)
    out[order] = ps
    mono = {'applied': True, 'keep_high': bool(keep_high),
            'x_min': float(xs[idx_min]), 'minp': minp,
            'x_max': float(xs[idx_max]), 'maxp': maxp, 'hold_max': hold_max}
    return out, mono


def apply_score(params, x):
    """pscore of new values x from classify_score params."""
    x = np.asarray(x, dtype=float).ravel()
    method = params.get('method')
    if method == 'gauss':
        with np.errstate(divide='ignore', invalid='ignore'):
            ps = (np.log(normpdf(x, params['mu1'], params['sig1'])) -
                  np.log(normpdf(x, params['mu2'], params['sig2'])))
        mono = params.get('mono')
        if mono and mono.get('applied'):
            valid = ~np.isnan(x)
            # order matters where the two holds overlap: min first, then max
            if mono['keep_high']:
                ps[valid & (x <= mono['x_min'])] = mono['minp']
            else:
                ps[valid & (x >= mono['x_min'])] = mono['minp']
            if mono['hold_max'] == 'high':
                ps[valid & (x >= mono['x_max'])] = mono['maxp']
            elif mono['hold_max'] == 'low':
                ps[valid & (x <= mono['x_max'])] = mono['maxp']
        return ps
    if method == 'fallback_mad':
        with np.errstate(divide='ignore', invalid='ignore'):
            if not params['scale'] > 0:
                return np.full(x.size, math.nan)
            return (x - params['med']) / params['scale'] * params['sign']
    return np.full(x.size, math.nan)


LAMBDA = np.arange(7201, dtype=float) * 0.5 - 150.0  # -150:0.5:3450


def calc_roc(x, y, pscore, ythr):
    """autoqc_calc_roc: AUC over the fixed pscore threshold grid, via rocp."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    pscore = np.asarray(pscore, dtype=float)
    out = {'auc': math.nan, 'fp0': math.nan, 'tp0': math.nan}
    if not np.any(np.isfinite(pscore)):
        return out
    idx = ~np.isnan(x) & ~np.isnan(y)
    ps = pscore[idx]
    if ythr > 0:
        truth = y[idx] > abs(ythr)
    else:
        truth = y[idx] < abs(ythr)
    npos = int(truth.sum())
    nneg = int((~truth).sum())
    # count pscore >= lambda (NaN never counts)
    sp = np.sort(ps[truth & ~np.isnan(ps)])
    sn = np.sort(ps[~truth & ~np.isnan(ps)])
    tpn = sp.size - np.searchsorted(sp, LAMBDA, side='left')
    fpn = sn.size - np.searchsorted(sn, LAMBDA, side='left')
    with np.errstate(invalid='ignore', divide='ignore'):
        tp = tpn / float(npos) if npos else np.full(LAMBDA.size, math.nan)
        fp = fpn / float(nneg) if nneg else np.full(LAMBDA.size, math.nan)
    tp = np.where(np.isnan(tp), 0.0, tp)
    fp = np.where(np.isnan(fp), 0.0, fp)
    out['auc'] = rocp(fp, tp, 8)
    i0 = int(np.argmin(np.abs(LAMBDA)))
    out['fp0'] = float(fp[i0])
    out['tp0'] = float(tp[i0])
    return out


def rocp(fp, tp, N=8):
    """autoqc_rocp: interpolated ROC closed by polyarea."""
    fp = np.asarray(fp, dtype=float).ravel()
    tp = np.asarray(tp, dtype=float).ravel()
    s = tp * 100 + fp
    o = np.argsort(s, kind='mergesort')
    tp = tp[o]
    fp = fp[o]
    fp = np.concatenate([[0.0], fp, [1.0]])
    tp = np.concatenate([[0.0], tp, [1.0]])
    s = tp * 100 + fp
    o = np.argsort(s, kind='mergesort')
    tp = tp[o]
    fp = fp[o]
    s = s[o]
    _, iu = np.unique(s, return_index=True)
    tp = tp[iu]
    fp = fp[iu]
    s = fp * 100 + tp
    o = np.argsort(s, kind='mergesort')
    tp = list(tp[o])
    fp = list(fp[o])
    ddy = 1e-4
    ddx = 1e-4
    xs = []
    ys = []
    for i in range(len(s) - 1):
        if tp[i + 1] == tp[i]:
            tp[i + 1] = tp[i] + ddy
        if fp[i + 1] == fp[i]:
            fp[i + 1] = fp[i] + ddx
        dx = (fp[i + 1] - fp[i]) / N
        dy = (tp[i + 1] - tp[i]) / N
        xs.append(colon_matlab(fp[i], dx, fp[i + 1]))
        ys.append(colon_matlab(tp[i], dy, tp[i + 1]))
    xx = np.concatenate(xs + [np.array([1.0, 0.0])])
    yy = np.concatenate(ys + [np.array([0.0, 0.0])])
    if xx.size != yy.size:
        return math.nan
    return float(abs(np.sum((np.roll(xx, -1) - xx) * (np.roll(yy, -1) + yy)) / 2))


###########################################################################
# screening

def calc_xcorr(X, Y):
    """autoqc_calc_xcorr: Pearson r of each x with each y, pairwise complete."""
    nx = X.shape[1]
    ny = Y.shape[1]
    R = np.full((nx, ny), math.nan)
    for i in range(nx):
        x = X[:, i]
        for j in range(ny):
            y = Y[:, j]
            ok = ~np.isnan(x) & ~np.isnan(y)
            if ok.any():
                R[i, j] = corr(x[ok], y[ok])
    return R


def calc_auc(x, y, ythr, monoflag=False):
    """One cell of autoqc_calc_xauc (y already fail-oriented with ythr)."""
    ok = ~np.isnan(x) & ~np.isnan(y)
    if not ok.any():
        return math.nan
    pscore, _ = classify_score(x, y, ythr, monoflag)
    return calc_roc(x, y, pscore, ythr)['auc']


def calc_metric_stats(x, y_fail, ythr):
    """n, n_fail, n_pass, frac_nan, cohens_d over the training rows."""
    with np.errstate(invalid='ignore'):
        fail = y_fail > ythr
        pas = y_fail <= ythr
    n = int(np.sum(~np.isnan(x) & ~np.isnan(y_fail)))
    frac_nan = float(np.sum(np.isnan(x))) / max(1, x.size)
    nf = int(np.sum(~np.isnan(x[fail])))
    npass = int(np.sum(~np.isnan(x[pas])))
    mf = nanmean(x[fail]) if fail.any() else math.nan
    sf = nanstd(x[fail]) if fail.any() else math.nan
    mp = nanmean(x[pas]) if pas.any() else math.nan
    sp = nanstd(x[pas]) if pas.any() else math.nan
    d = math.nan
    with np.errstate(invalid='ignore'):
        sd_pool = math.sqrt(((nf - 1) * sf ** 2 + (npass - 1) * sp ** 2) /
                            max(1, nf + npass - 2)) \
            if not (math.isnan(sf) or math.isnan(sp)) else math.nan
    if sd_pool > 0:
        d = (mf - mp) / sd_pool
    return n, nf, npass, frac_nan, d


def select_metrics(Xt, xnames, R_col, y_train_fail, ythr, y_auc, ythr_auc, parms, log):
    """select_metrics of autoqc_select_all: rthr, min_n/min_nfail, athr."""
    absR = np.abs(R_col)
    with np.errstate(invalid='ignore'):
        idx_rthr = np.nonzero(absR >= parms['rthr'])[0]
    stats = [calc_metric_stats(Xt[:, i], y_train_fail, ythr) for i in idx_rthr]
    ok = [k for k, st in enumerate(stats)
          if st[0] >= parms['min_n'] and st[1] >= parms['min_nfail']]
    for k, st in enumerate(stats):
        if k not in ok:
            log('excluding %s: n = %d, n_fail = %d (min_n = %d, min_nfail = %d)' %
                (xnames[idx_rthr[k]], st[0], st[1], parms['min_n'], parms['min_nfail']))
    if parms['max_auc_flag'] or parms['athr'] > 0:
        aucs = np.array([calc_auc(Xt[:, i], y_auc, ythr_auc, False)
                         for i in idx_rthr])
        with np.errstate(invalid='ignore'):
            keep = sorted(set(np.nonzero(aucs >= parms['athr'])[0]) & set(ok))
        if parms['max_auc_flag'] and len(keep) > 1:
            keep = [keep[int(np.argmax(aucs[keep]))]]
        idx_sel = idx_rthr[keep] if keep else np.zeros(0, dtype=int)
    else:
        aucs = np.full(idx_rthr.size, math.nan)
        idx_sel = idx_rthr[ok] if ok else np.zeros(0, dtype=int)
    rows = []
    for k, i in enumerate(idx_rthr):
        n, nf, npass, frac_nan, d = stats[k]
        if i in idx_sel and frac_nan > parms['max_frac_nan']:
            log('WARNING: %s is missing for %0.1f%% of the training set (n = %d)' %
                (xnames[i], 100 * frac_nan, n))
        rows.append({'xname': xnames[i], 'n': n, 'n_fail': nf, 'n_pass': npass,
                     'frac_nan': frac_nan, 'R': float(R_col[i]),
                     'AUC': float(aucs[k]), 'cohens_d': d,
                     'selected': 1.0 * (i in idx_sel)})
    if parms['rthr_margin'] > 0:
        with np.errstate(invalid='ignore'):
            near = np.nonzero((absR >= parms['rthr'] - parms['rthr_margin']) &
                              (absR < parms['rthr']))[0]
        for i in near:
            n, nf, npass, frac_nan, d = calc_metric_stats(Xt[:, i], y_train_fail, ythr)
            rows.append({'xname': xnames[i], 'n': n, 'n_fail': nf, 'n_pass': npass,
                         'frac_nan': frac_nan, 'R': float(R_col[i]),
                         'AUC': calc_auc(Xt[:, i], y_auc, ythr_auc, False),
                         'cohens_d': d, 'selected': 0.0})
    return [int(i) for i in idx_sel], rows


###########################################################################
# logistic regression (glmfit binomial logit, R2023a)

def _back_sub(R, b):
    n = R.shape[0]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(R[i, i + 1:], x[i + 1:])) / R[i, i]
    return x


def _qr_pivot(A):
    """Column-pivoted Householder QR: (rank as glmfit counts it, pivot order)."""
    A = np.array(A, dtype=float)
    m, n = A.shape
    k = min(m, n)
    piv = np.arange(n)
    norms = np.sum(A * A, axis=0)
    for j in range(k):
        jm = j + int(np.argmax(norms[j:]))
        if jm != j:
            A[:, [j, jm]] = A[:, [jm, j]]
            norms[[j, jm]] = norms[[jm, j]]
            piv[[j, jm]] = piv[[jm, j]]
        v = A[j:, j].copy()
        nv = np.linalg.norm(v)
        if nv == 0:
            continue
        alpha = -nv if v[0] >= 0 else nv
        v[0] -= alpha
        v /= np.linalg.norm(v)
        A[j:, j:] -= 2.0 * np.outer(v, v @ A[j:, j:])
        norms[j + 1:] = np.sum(A[j + 1:, j + 1:] ** 2, axis=0)
    d = np.abs(np.diag(A[:k, :k]))
    if d.size == 0:
        return 0, piv
    return int(np.sum(d > d[0] * max(m, n) * EPS)), piv


def _diagnose_separation(eta, y):
    """glmfit's diagnoseSeparation, as a flag."""
    o = np.argsort(eta, kind='mergesort')
    x = eta[o]
    p = y[o]
    if np.all(p == p[0]) or x[0] == x[-1]:
        return False
    no_front = 0 < p[0] < 1
    no_end = 0 < p[-1] < 1
    if no_front and no_end:
        return False
    dx = 100 * max(np.spacing(abs(x[0])), np.spacing(abs(x[-1])))
    n = p.size
    if no_front:
        A = 0
    else:
        A = int(np.nonzero(p != p[0])[0][0])
        cutoff = x[A] - dx
        A = int(np.sum(x[:A] < cutoff))
    if no_end:
        B = n + 1
    else:
        B = int(np.nonzero(p != p[-1])[0][-1]) + 2
        cutoff = x[B - 2] + dx
        B = (n + 1) - int(np.sum(x[B - 1:] > cutoff))
    if A + 1 < B - 1:
        if x[B - 2] - x[A] > dx:
            return False
    return True


def glmfit_logit(X, y, iter_lim=100, convcrit=1e-6):
    """glmfit(X, [y 1], 'binomial', 'link', 'logit'); returns a dict."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n = X.shape[0]
    x = np.column_stack([np.ones(n), X])
    ncolx = x.shape[1]
    rankx, piv = _qr_pivot(x)
    perm = np.arange(ncolx)
    rank_deficient = rankx < ncolx
    if rank_deficient:
        # keep the first rankx pivoted columns, as glmfit does
        perm = piv[:rankx]
        x = x[:, perm]
    p = x.shape[1]
    mu = (y + 0.5) / 2.0
    eta = np.log(mu / (1 - mu))
    b = np.zeros(p)
    seps = math.sqrt(EPS)
    lo = math.log(EPS)
    hi = -lo
    it = 0
    Rw = None
    sqrtirls = None
    while it <= iter_lim:
        it += 1
        deta = 1.0 / (mu * (1 - mu))
        z = eta + (y - mu) * deta
        sqrtirls = np.abs(deta) * np.sqrt(mu) * np.sqrt(1 - mu)
        sqrtw = 1.0 / sqrtirls
        wtol = np.max(sqrtw) * EPS ** (2.0 / 3.0)
        t = (sqrtw < wtol) & (sqrtw != 0)
        sqrtw[t] = wtol
        b_old = b
        Q, Rw = np.linalg.qr(x * sqrtw[:, None])
        b = _back_sub(Rw, Q.T @ (z * sqrtw))
        eta = x @ b
        mu = 1.0 / (1.0 + np.exp(-np.clip(eta, lo, hi)))
        mu = np.clip(mu, EPS, 1 - EPS)
        if not np.any(np.abs(b - b_old) > convcrit * np.maximum(seps, np.abs(b_old))):
            break
    iter_limit = it > iter_lim
    separated = _diagnose_separation(eta, y) if iter_limit else False
    RI = np.vstack([_back_sub(Rw, e) for e in np.eye(p)]).T
    C = RI @ RI.T
    se = np.sqrt(np.diag(C))
    bb = np.zeros(ncolx)
    bb[perm] = b
    tt = np.full(ncolx, math.nan)
    tt[perm] = b / se
    pv = np.array([math.erfc(abs(v) / math.sqrt(2)) if not math.isnan(v) else math.nan
                   for v in tt])
    return {'beta': bb, 'p': pv, 'iter': it, 'iter_limit': bool(iter_limit),
            'separated': bool(separated), 'rank_deficient': bool(rank_deficient)}


def train_multivar(X, y, ythr=2.0, max_nvar=10, sthr=0.05):
    """autoqc_train_multivar_prediction: z-score, SVD, logistic fit."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    nx = X.shape[1]
    out = {'xmean': None, 'xstd': None, 'coef': None, 'eigvals': None, 'nvar': nx}
    if nx > 1:
        ok = ~np.isnan(np.sum(X, axis=1))
        Xt = X[ok]
        xmean = np.mean(Xt, axis=0)
        xstd = np.std(Xt, axis=0, ddof=1)
        with np.errstate(invalid='ignore', divide='ignore'):
            Xn = (Xt - xmean) / xstd
        Xn[:, xstd == 0] = 0
        _, S, Vt = np.linalg.svd(Xn, full_matrices=False)
        eig = S ** 2
        eig = eig / np.max(eig)
        nvar = int(np.sum(eig > sthr))
        nvar = min(max_nvar, min(nx, nvar))
        coef = Vt.T[:, :nvar]
        Z = np.full((X.shape[0], nvar), math.nan)
        Z[ok] = Xn @ coef
        X = Z
        out.update({'xmean': xmean, 'xstd': xstd, 'coef': coef, 'eigvals': eig,
                    'nvar': nvar})
    with np.errstate(invalid='ignore'):
        Y = 1.0 * (y > ythr)
    Y[np.isnan(y)] = math.nan
    ok = ~np.isnan(np.sum(X, axis=1)) & ~np.isnan(Y)
    fit = glmfit_logit(X[ok], Y[ok])
    out.update(fit)
    out['n_fit'] = int(ok.sum())
    out['yp'] = np.column_stack([np.ones(X.shape[0]), X]) @ fit['beta']
    return out


def test_multivar(X, xmean, xstd, coef, beta):
    """autoqc_test_multivar_prediction."""
    X = np.asarray(X, dtype=float)
    if X.shape[1] > 1:
        xstd = np.asarray(xstd, dtype=float)
        with np.errstate(invalid='ignore', divide='ignore'):
            X = (X - np.asarray(xmean, dtype=float)) / xstd
        X[:, xstd == 0] = 0
        X = X @ np.asarray(coef, dtype=float)
    return np.column_stack([np.ones(X.shape[0]), X]) @ np.asarray(beta, dtype=float)


###########################################################################
# calibration (autoqc_select_cases calib_pscore)

def fit_calibration(pscore, y, ythr, idx_calib):
    """Isotonic fit of failure on pscore over idx_calib; None if not fittable."""
    p = np.clip(np.asarray(pscore, dtype=float), -1e6, 1e6)
    idx = np.asarray(idx_calib, dtype=int)
    idx = idx[~np.isnan(y[idx]) & ~np.isnan(p[idx])]
    pc = p[idx]
    f = 1.0 * (y[idx] > ythr)
    if not f.any() or f.all():
        return None
    up, ic = np.unique(pc, return_inverse=True)
    bn = list(np.bincount(ic).astype(float))
    bsum = list(np.bincount(ic, weights=f))
    bx = list(up * np.array(bn))
    k = 0
    while k < len(bsum) - 1:
        if bsum[k] / bn[k] > bsum[k + 1] / bn[k + 1]:
            bsum[k] += bsum[k + 1]
            bn[k] += bn[k + 1]
            bx[k] += bx[k + 1]
            del bsum[k + 1], bn[k + 1], bx[k + 1]
            if k > 0:
                k -= 1
        else:
            k += 1
    rate = np.array(bsum) / np.array(bn)
    xk = np.array(bx) / np.array(bn)
    xk, iu = np.unique(xk, return_index=True)
    rate = rate[iu]
    rmin = 0.5 / idx.size
    rate = np.minimum(np.maximum(rate, rmin), 1 - rmin)
    return {'xk': xk, 'rate': rate}


def apply_calibration(calib, pscore, calib_k):
    """logit(rate) + log(calib_k); all NaN when calib is None."""
    p = np.clip(np.asarray(pscore, dtype=float), -1e6, 1e6)
    out = np.full(p.size, math.nan)
    if calib is None:
        return out
    xk = np.asarray(calib['xk'], dtype=float)
    rate = np.asarray(calib['rate'], dtype=float)
    ival = ~np.isnan(p)
    pv = p[ival]
    if xk.size > 1:
        r = np.interp(pv, xk, rate)
        r[pv < xk[0]] = rate[0]
        r[pv > xk[-1]] = rate[-1]
    else:
        r = np.full(pv.size, rate[0])
    out[ival] = np.log(r / (1 - r)) + math.log(calib_k)
    return out


###########################################################################
# per-label selection (autoqc_select_cases) and measures

def select_label_cases(pscore, y, ythr, pthr, filt='new', nselect=math.inf,
                       idx_allowed=None):
    """Rows selected for one label, in descending pscore (stable)."""
    ps = np.asarray(pscore, dtype=float).copy()
    with np.errstate(invalid='ignore'):
        ps[ps <= pthr] = math.nan
    if idx_allowed is not None:
        mask = np.ones(ps.size, dtype=bool)
        mask[np.asarray(idx_allowed, dtype=int)] = False
        ps[mask] = math.nan
    o = _stable_desc(np.where(np.isnan(ps), -math.inf, ps))
    o = o[~np.isnan(ps[o])]
    yo = y[o]
    with np.errstate(invalid='ignore'):
        if filt == 'new':
            keep = np.isnan(yo)
        elif filt == 'old':
            keep = ~np.isnan(yo)
        elif filt == 'any':
            keep = np.ones(o.size, dtype=bool)
        elif filt == 'pass':
            keep = yo <= ythr
        elif filt == 'fail':
            keep = yo > ythr
        else:
            raise ValueError('unrecognized filter')
    cases = o[keep]
    if nselect < cases.size:
        cases = cases[:int(nselect)]
    return cases, ps[cases]


def _safe_ratio(a, b):
    return a / b if b > 0 else math.nan


def label_measures(y, pscore, ythr, idx_train, idx_allowed):
    """calc_label_measures of autoqc_select_all (in-sample, at pscore > 0)."""
    m = {}
    yt = y[idx_train]
    pt = pscore[idx_train]
    ok = ~np.isnan(yt)
    with np.errstate(invalid='ignore'):
        fail = ok & (yt > ythr)
        sel = ok & (pt > 0)
    n = int(ok.sum())
    nf = int(fail.sum())
    ns = int(sel.sum())
    nsf = int((sel & fail).sum())
    m['n'] = n
    m['n_fail'] = nf
    m['n_sel'] = ns
    m['n_notsel'] = n - ns
    m['base_rate'] = _safe_ratio(nf, n)
    m['precision'] = _safe_ratio(nsf, ns)
    m['boost'] = _safe_ratio(m['precision'], m['base_rate'])
    m['false_omission_rate'] = _safe_ratio(nf - nsf, n - ns)
    if ns > 0 and nsf == 0:
        m['fa_per_fail'] = math.inf
    else:
        m['fa_per_fail'] = _safe_ratio(ns - nsf, nsf)
    m['false_discovery_rate'] = 1 - m['precision']
    m['recall'] = _safe_ratio(nsf, nf)
    m['n_sel_fail'] = nsf
    idx_new = np.setdiff1d(np.asarray(idx_allowed), np.asarray(idx_train))
    m['n_new'] = int(idx_new.size)
    with np.errstate(invalid='ignore'):
        m['n_new_sel'] = int(np.sum(pscore[idx_new] > 0))
    m['shift_ratio'] = _safe_ratio(_safe_ratio(m['n_new_sel'], m['n_new']),
                                   _safe_ratio(ns, n))
    return m


def label_status(m, parms, log, yname):
    fails = []
    if m['fa_per_fail'] > parms['max_fapf']:
        fails.append('max_fapf')
    if m['boost'] < parms['min_boost']:
        fails.append('min_boost')
    if parms['max_marginal_fapf'] is not None and m['n_sel_fail'] < parms['min_nfail']:
        fails.append('min_nfail')
    if fails:
        log('WARNING: excluding multivar vs %s (%s)' % (yname, '+'.join(fails)))
        return '+'.join(fails)
    if m['shift_ratio'] > parms['max_shift_ratio'] and m['n_new_sel'] >= 10:
        log('WARNING: multivar vs %s selects %d of %d cases outside the training set'
            % (yname, m['n_new_sel'], m['n_new']))
    return 'selected'


###########################################################################
# scanner regressors (scanner_flag 1)

def scanner_dummies(manufacturer, levels=None):
    """Dummy columns for manufacturer, dropping the most common level.

    levels: (all_levels, kept_levels) from a saved model, to code new rows alike.
    """
    vals = np.array(['NaN' if (v is None or (isinstance(v, float) and math.isnan(v)))
                     else str(v) for v in manufacturer])
    if levels is None:
        all_levels = sorted(set(vals) - {'NaN'})
        counts = [int(np.sum(vals == lv)) for lv in all_levels]
        imax = int(np.argmax(counts))
        kept = [lv for i, lv in enumerate(all_levels) if i != imax and counts[i] > 0]
    else:
        all_levels, kept = levels
    D = np.column_stack([1.0 * (vals == lv) for lv in kept]) if kept \
        else np.zeros((vals.size, 0))
    D[vals == 'NaN', :] = math.nan
    return D, (list(all_levels), list(kept))


###########################################################################
# train / score / select

def _drop_uninformative(X, names):
    """get_datamat: drop all-NaN and constant columns (over every row)."""
    keep = []
    for i in range(X.shape[1]):
        d = X[:, i]
        v = d[~np.isnan(d)]
        if v.size == 0 or np.unique(v).size < 2:
            continue
        keep.append(i)
    return X[:, keep], [names[i] for i in keep]


def _jsonable(v):
    if isinstance(v, np.ndarray):
        return [_jsonable(a) for a in v.tolist()] if v.ndim else _jsonable(v.item())
    if isinstance(v, (list, tuple)):
        return [_jsonable(a) for a in v]
    if isinstance(v, dict):
        return {k: _jsonable(a) for k, a in v.items()}
    if isinstance(v, (np.floating, float)):
        f = float(v)
        if math.isnan(f):
            return None
        if math.isinf(f):
            return 'Inf' if f > 0 else '-Inf'
        return f
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.bool_):
        return bool(v)
    return v


def _nan_label_row(yname, status, nmetrics, y_train, ythr, max_abs_R, parms):
    n = int(np.sum(~np.isnan(y_train)))
    with np.errstate(invalid='ignore'):
        nf = int(np.sum(y_train > ythr))
    m = {'n': n, 'n_fail': nf, 'base_rate': _safe_ratio(nf, n)}
    return _label_row(yname, status, nmetrics, m, ythr, max_abs_R, parms)


def _label_row(yname, status, nmetrics, m, ythr, max_abs_R, parms):
    row = {c: math.nan for c in LABEL_COLS}
    row.update({k: v for k, v in m.items() if k in row})
    mmf = parms['max_marginal_fapf']
    row.update({
        'yname': yname, 'xname': 'multivar' if status not in
        ('no_fails', 'all_fails', 'no_metrics') else '',
        'status': status, 'n_metrics': nmetrics,
        'max_abs_R': max_abs_R,
        'near_rthr': 1.0 * (max_abs_R < parms['rthr'] + parms['rthr_margin'])
        if not math.isnan(max_abs_R) else 0.0,
        'ythr': ythr, 'max_fapf': parms['max_fapf'], 'min_boost': parms['min_boost'],
        'rthr': parms['rthr'], 'rthr_margin': parms['rthr_margin'],
        'max_marginal_fapf': math.nan if mmf is None else mmf})
    return row


def train(X, xnames, Y, labels, ids=None, parms=None, idx_train=None,
          idx_allowed=None, manufacturer=None, log=None):
    """Screen, fit and score each label; see the module docstring.

    X: rows x metrics; Y: rows x labels (same order as labels); ids: row ids,
    saved as the training ids.  idx_train defaults to the rows with every
    label (as autoqc_select_all); idx_allowed to every row.
    """
    parms = check_parms(parms)
    log = log or (lambda s: None)
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    nrows = X.shape[0]
    xnames = list(xnames)
    specs = [_label_spec(lb) for lb in labels]
    if Y.ndim == 1:
        Y = Y[:, None]
    if idx_train is None:
        idx_train = np.nonzero(~np.isnan(np.sum(Y, axis=1)))[0]
    idx_train = np.asarray(idx_train, dtype=int)
    if idx_allowed is None:
        idx_allowed = np.arange(nrows)
    idx_allowed = np.asarray(idx_allowed, dtype=int)

    R = calc_xcorr(X[idx_train], Y[idx_train])

    model = {
        'version': 1,
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
        'parms': {k: parms[k] for k in sorted(parms)},
        'xnames': xnames,
        'n_rows': int(nrows),
        'training_ids': [str(ids[i]) for i in idx_train] if ids is not None else None,
        'mahal_reference': None,
        'labels': [],
    }
    screen_rows = []
    label_rows = []
    pscores = {}
    yps = {}
    for j, spec in enumerate(specs):
        y_raw = Y[:, j]
        y, ythr, yname, invert = prepare_label(y_raw, spec, parms, idx_train)
        # labels outside idx_train do not reach the score, ROC or pos_rate
        y_out = np.ones(nrows, dtype=bool)
        y_out[idx_train] = False
        y[y_out] = math.nan
        lm = {'name': spec['name'], 'yname': yname, 'type': spec['type'],
              'inverted': invert, 'ythr': ythr, 'status': None}
        model['labels'].append(lm)
        y_train = y[idx_train]
        with np.errstate(invalid='ignore'):
            nfail = int(np.sum(y_train > ythr))
        with np.errstate(invalid='ignore'):
            all_fail = bool(np.all(y_train > ythr))
        if nfail == 0 or all_fail:
            status = 'no_fails' if nfail == 0 else 'all_fails'
            log('%s for %s in training set' % (status, yname))
            lm['status'] = status
            label_rows.append(_nan_label_row(yname, status, 0, y_train, ythr,
                                             math.nan, parms))
            continue

        # screening uses parms ythr and the fail definition of check_ythr
        y_auc, ythr_auc = _auc_label(y_raw[idx_train], spec, parms,
                                     y_train, ythr)
        idx_xs, rows = select_metrics(X[idx_train], xnames, R[:, j], y_train,
                                      ythr, y_auc, ythr_auc, parms, log)
        for r in rows:
            r.update({'yname': spec['name'], 'rthr': parms['rthr'],
                      'athr': parms['athr'], 'min_n': parms['min_n'],
                      'min_nfail': parms['min_nfail'],
                      'rthr_margin': parms['rthr_margin']})
        screen_rows.extend(rows)
        lm['screen'] = rows
        if not idx_xs:
            log('selected 0 metrics for %s' % yname)
            lm['status'] = 'no_metrics'
            label_rows.append(_nan_label_row(yname, 'no_metrics', 0, y_train, ythr,
                                             math.nan, parms))
            continue
        max_abs_R = float(np.max(np.abs(R[idx_xs, j])))
        log('selected %d metrics for %s: %s' %
            (len(idx_xs), yname, ' '.join(xnames[i] for i in idx_xs)))

        Xa, xkept = _drop_uninformative(X[:, idx_xs], [xnames[i] for i in idx_xs])
        scanner_levels = None
        if parms['scanner_flag'] == 1:
            if manufacturer is None:
                raise ValueError('scanner_flag 1 needs manufacturer')
            Xs, scanner_levels = scanner_dummies(manufacturer)
            Xa = np.column_stack([Xs, Xa])
        if Xa.shape[1] == 0 or not np.any(~np.isnan(np.sum(Xa[idx_train], axis=1))):
            log('WARNING: no complete training rows for %s' % yname)
            lm['status'] = 'fit_failed'
            label_rows.append(_nan_label_row(yname, 'fit_failed', len(idx_xs), y_train,
                                             ythr, max_abs_R, parms))
            continue
        fit = train_multivar(Xa[idx_train], y_train, ythr, parms['max_nvar'],
                             parms['sthr'])
        lm.update({'metrics': xkept, 'scanner_levels': scanner_levels,
                   'xmean': fit['xmean'], 'xstd': fit['xstd'], 'coef': fit['coef'],
                   'eigvals': fit['eigvals'], 'nvar': fit['nvar'],
                   'beta': fit['beta'], 'p': fit['p'], 'glm_iter': fit['iter'],
                   'iter_limit': fit['iter_limit'], 'separated': fit['separated'],
                   'rank_deficient': fit['rank_deficient'], 'n_fit': fit['n_fit']})
        if fit['iter_limit']:
            log('WARNING: glmfit iteration limit for %s (separated %d)' %
                (yname, fit['separated']))
        pv = fit['p'][1:]
        if np.all(pv > parms['alpha']):
            if fit['separated']:
                status = 'separated'
                log('WARNING: perfectly separated logistic fit for %s (%d fails), '
                    'rejected' % (yname, nfail))
            else:
                status = 'fit_failed'
                log('WARNING: multivariate prediction failed for %s' % yname)
            lm['status'] = status
            row = _nan_label_row(yname, status, len(idx_xs), y_train, ythr,
                                 max_abs_R, parms)
            row.update({'xname': 'multivar', 'glm_iter': fit['iter'],
                        'iter_limit': 1.0 * fit['iter_limit'],
                        'separated': 1.0 * fit['separated']})
            label_rows.append(row)
            continue
        if idx_train.size < nrows:
            yp = test_multivar(Xa, fit['xmean'], fit['xstd'], fit['coef'], fit['beta'])
        else:
            yp = fit['yp']

        pscore, sparams = classify_score(yp, y, ythr, parms['monoflag'])
        roc = calc_roc(yp, y, pscore, ythr)
        nonan = ~np.isnan(y)
        with np.errstate(invalid='ignore'):
            pos_rate = float(np.sum(y[nonan] > ythr)) / max(1, int(nonan.sum()))
        calib = None
        calib_thr = math.nan
        if parms['max_marginal_fapf'] is not None:
            calib = fit_calibration(pscore, y, ythr, idx_train)
            praw = pscore
            pscore = apply_calibration(calib, pscore, parms['max_marginal_fapf'])
            with np.errstate(invalid='ignore'):
                isel = pscore >= 0
            if isel.any():
                calib_thr = float(np.min(praw[isel]))
        m = label_measures(y, pscore, ythr, idx_train, idx_allowed)
        m.update({'roc_auc': roc['auc'], 'pos_rate': pos_rate, 'calib_thr': calib_thr})
        status = label_status(m, parms, log, yname)
        lm.update({'status': status, 'score': sparams, 'calibration': calib,
                   'roc_auc': roc['auc'], 'pos_rate': pos_rate,
                   'calib_thr': calib_thr})
        row = _label_row(yname, status, len(idx_xs), m, ythr, max_abs_R, parms)
        n_unscored = int(np.sum(np.isnan(yp[idx_allowed])))
        with np.errstate(invalid='ignore'):
            accept = len(np.unique(idx_allowed[pscore[idx_allowed] >= 0]))
        row.update({'n_unscored': n_unscored, 'glm_iter': fit['iter'],
                    'iter_limit': 1.0 * fit['iter_limit'],
                    'separated': 1.0 * fit['separated'],
                    'score_method': sparams['method'],
                    'accept_frac': _safe_ratio(accept, idx_allowed.size)})
        label_rows.append(row)
        pscores[spec['name']] = pscore
        yps[spec['name']] = yp
    return {'model': model, 'pscore': pscores, 'yp': yps, 'metrics_rows': screen_rows,
            'labels_rows': label_rows, 'R': R, 'idx_train': idx_train,
            'idx_allowed': idx_allowed}


def score(model, X, xnames, manufacturer=None, return_yp=False):
    """Apply a saved model to rows of X: {label name: pscore} (and yp)."""
    X = np.asarray(X, dtype=float)
    col = {n: i for i, n in enumerate(xnames)}
    out = {}
    yps = {}
    for lm in model['labels']:
        if 'beta' not in lm or lm.get('status') in (None, 'fit_failed', 'separated'):
            continue
        missing = [n for n in lm['metrics'] if n not in col]
        if missing:
            raise ValueError('label %s: missing metrics %s' %
                             (lm['name'], ', '.join(missing)))
        Xa = X[:, [col[n] for n in lm['metrics']]]
        if lm.get('scanner_levels'):
            Xs, _ = scanner_dummies(manufacturer, lm['scanner_levels'])
            Xa = np.column_stack([Xs, Xa])
        yp = test_multivar(Xa, lm['xmean'], lm['xstd'], lm['coef'], lm['beta'])
        ps = apply_score(lm['score'], yp)
        mmf = model['parms'].get('max_marginal_fapf')
        if mmf is not None:
            ps = apply_calibration(lm['calibration'], ps, mmf)
        out[lm['name']] = ps
        yps[lm['name']] = yp
    return (out, yps) if return_yp else out


def select(pscore, Y, model, parms=None, idx_allowed=None, label_names=None, yp=None):
    """Selected rows over labels: list of dicts (row, p, xname, x, yname, y), by p.

    yp: {label name: linear predictor}, reported as x (as autoqc_select_all).

    Per label: pthr, idx_allowed, filter and nselect (autoqc_select_cases);
    then the pos_rate weighting, one row per case at its highest p, and nall.
    """
    parms = check_parms(parms if parms is not None else
                        {k: v for k, v in model['parms'].items()
                         if k in DEFAULT_PARMS})
    Y = np.asarray(Y, dtype=float)
    if Y.ndim == 1:
        Y = Y[:, None]
    names = label_names or [lm['name'] for lm in model['labels']]
    per_label = []
    for lm in model['labels']:
        if lm.get('status') != 'selected' or lm['name'] not in pscore:
            continue
        j = names.index(lm['name'])
        y = Y[:, j].copy()
        if lm.get('inverted'):
            ok = ~np.isnan(y)
            y[ok] = 1.0 * (y[ok] == 0)
        cases, ps = select_label_cases(pscore[lm['name']], y, lm['ythr'],
                                       parms['pthr'], parms['filter'],
                                       parms['nselect'], idx_allowed)
        if cases.size:
            per_label.append((lm, cases, ps, y))
    max_pos_rate = EPS
    for lm, _, _, _ in per_label:
        if lm['pos_rate'] > max_pos_rate:
            max_pos_rate = lm['pos_rate']
    rows = []
    weights = {}
    for lm, cases, ps, y in per_label:
        if parms['weighted_flag']:
            w = lm['pos_rate'] / max_pos_rate
            weights[lm['name']] = w
            nsel = int(min(cases.size, round_half_away(w * parms['nselect'])))
            cases = cases[:nsel]
            ps = ps[:nsel]
        for c, p in zip(cases, ps):
            xv = float(yp[lm['name']][c]) if yp and lm['name'] in yp else math.nan
            rows.append({'row': int(c), 'p': float(p), 'xname': 'multivar', 'x': xv,
                         'yname': lm['yname'], 'y': float(y[c])})
    # one row per case at its highest p, then by p; stable, as sortrows/unique
    rows.sort(key=lambda r: (r['row'], -r['p']))
    seen = set()
    uniq = []
    for r in rows:
        if r['row'] in seen:
            continue
        seen.add(r['row'])
        uniq.append(r)
    uniq.sort(key=lambda r: -r['p'])
    if parms['nall'] < len(uniq):
        uniq = uniq[:int(parms['nall'])]
    return {'rows': uniq, 'weights': weights, 'max_pos_rate': max_pos_rate}


def model_to_json(model, fname):
    with open(fname, 'w') as f:
        json.dump(_jsonable(model), f, indent=1)


def _unjson(v):
    if isinstance(v, list):
        return [_unjson(a) for a in v]
    if isinstance(v, dict):
        return {k: _unjson(a) for k, a in v.items()}
    if v == 'Inf':
        return math.inf
    if v == '-Inf':
        return -math.inf
    return v


def model_from_json(fname):
    """Load a saved model; arrays come back as numpy (None -> NaN)."""
    with open(fname) as f:
        model = _unjson(json.load(f))
    for lm in model['labels']:
        for k in ('xmean', 'xstd', 'coef', 'eigvals', 'beta', 'p'):
            if lm.get(k) is not None:
                lm[k] = np.array(lm[k], dtype=float)
        sp = lm.get('score')
        if sp:
            for k in ('mu1', 'sig1', 'mu2', 'sig2', 'med', 'mad', 'sign'):
                if k in sp and sp[k] is None:
                    sp[k] = math.nan
        if lm.get('calibration'):
            lm['calibration'] = {k: np.array(v, dtype=float)
                                 for k, v in lm['calibration'].items()}
    return model


def write_rows(rows, cols, fname):
    """Tab-separated table of dict rows (NaN as n/a)."""
    def fmt(v):
        if isinstance(v, str):
            return v
        if v is None:
            return 'n/a'
        if isinstance(v, (bool, np.bool_)):
            return str(int(v))
        f = float(v)
        if math.isnan(f):
            return 'n/a'
        if math.isinf(f):
            return 'Inf' if f > 0 else '-Inf'
        if f == int(f) and abs(f) < 1e15:
            return str(int(f))
        return repr(f)
    with open(fname, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(fmt(r.get(c, math.nan)) for c in cols) + '\n')
