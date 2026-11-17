#!/usr/bin/env python3
"""
Set the ABCD imaging inclusion flags (mr_y_qc__incl) from the compiled
tabulated data, one row per session.

Run across sessions after cmig_compile_bids_tabdat.py (with its release
filter, --id_events_file) and abcd_write_legacy_tabdat.py. The rows are those
of mr_y_adm__info. Reimplements abcd_imgincl.py (release 7.0 imgincl) with the
release 8.0 variable names; the criteria are in CRITERIA, one place.

Public function
---------------
set_imgincl_abcd(tabdir, outdir, legacy_dir=None, report_flag=0,
                 forceflag=0, verbose=1)

  tabdir      directory of the compiled tables, <table>.tsv
  legacy_dir  directory of the frozen _desc-legacy tables (E-Prime timing
              mismatch override, mr_y_qc__raw__tfmri__<task>_desc-legacy.tsv);
              default tabdir

  Returns (fname, errcode, errmsg), fname the mr_y_qc__incl.tsv written.
  A missing input table or column, the manual post QC tables included, is
  an error (errcode 1, errmsg naming each one).

Criteria
--------
A flag is 1 when the session passes every criterion of the flag. A missing
input value fails its test, except: a missing manual post QC score passes
(most sessions have no manual post QC), and a missing sst coding error
indicator passes, as in release 7.0. The GE start delay test fails only a GE
session with a start delay min or max of 6.4 (a missing delay passes); the
E-Prime test passes on a match or on the legacy ignore flag. The
series_count and eprime run valid_count criteria are implied by others (the
flags are unchanged by them); they give the cause of an exclusion.

  t1_series, <type>_series   mr_y_qc__raw__<type>__series_count > 0
  t1_ok, <type>_ok           mr_y_qc__raw__<type>__ok_count > 0
  fsurf_man, t2_man,         mr_y_qc__post__man__<fsurf|t2|dmri|fmri>_score
    dmri_man, fmri_man         != 0 (missing passes)
  dmri_reps                  mr_y_qc__raw__dmri__frame__ok_mean >= 103, or
                               Philips with dmri ok_count >= 2 and 51 reps
  dmri_b0unwarp, fmri_b0unwarp
                             mr_y_qc__post__aut__<dmri|fmri>__b0__unwarp_indicator == 1
  rsfmri_frames              mr_y_qc__mot__rsfmri__frame__censor_count > 375
  <task>_eprime_match        mr_y_qc__raw__tfmri__<task>__eprime__match_indicator
                               == 1 or (legacy) ...__eprime__tdiff__ign_indicator == 1
  <task>_eprime_valid        mr_y_qc__raw__tfmri__<task>__eprime__scan__run__valid_count > 0
  <task>_ge_delay            not (GE and mr_y_tfmri__<task>__beh__start__delay_min
                               or _max == 6.4)
  sst_coderr                 mr_y_tfmri__sst__beh__coderr_indicator == 0 (missing passes)
  <task>_beh_qc              mr_y_tfmri__<task>__beh__qc_indicator == 1
  <task>_dof                 mr_y_qc__mot__tfmri__<task>_dof > 200
  mid_trials                 mr_y_tfmri__mid__beh__trial_count == 100
  <flag>_roi                 ROI value not missing: mr_y_smri__t1__aseg__cwm__lh_mean,
                               mr_y_smri__t2__aseg__cwm__lh_mean, mr_y_rsi__rnd__at_wmean,
                               mr_y_rsfmri__corr__gpnet__def__sal_mean,
                               mr_y_tfmri__mid__arvn__aseg__cwm__lh_beta,
                               mr_y_tfmri__nback__0b__aseg__cbwm__lh_beta,
                               mr_y_tfmri__sst__cg__aseg__cwm__lh_beta

Not criteria (release 7.0 criteria omitted for now): the T1w registration
rigid scores (fMRI, dMRI, T2) and the dMRI and fMRI FOV cutoffs.

Outputs (in outdir)
-------
  mr_y_qc__incl.tsv               participant_id, session_id and the 7
                                  mr_y_qc__incl__*_indicator, 0/1
  mr_y_qc__incl_desc-criteria.tsv participant_id, session_id and, per flag,
                                  one 0/1 column per criterion (1 = passes),
                                  <flag>__<criterion>, e.g.
                                  tfmri__mid__mid_eprime_match, then the
                                  flag, <flag>__incl
  imgincl_report.txt              with report_flag: per flag and criterion,
                                  in order, the number of sessions failing
                                  it, failing it in sequence (after passing
                                  every criterion before it), failing only
                                  it, with it missing; and n included

Created:  11/09/2026 by Don Hagler
Last mod: 11/09/2026 by Don Hagler
"""

import argparse
import os
import sys
from collections import OrderedDict

import pandas as pd

APPNAME = 'cmig_set_imgincl_abcd'

KEYS = ['participant_id', 'session_id']
ROW_TABLE = 'mr_y_adm__info'
OUT_TABLE = 'mr_y_qc__incl'
CRITERIA_SUFFIX = '_desc-criteria'
REPORT_NAME = 'imgincl_report.txt'
LEGACY_SUFFIX = '_desc-legacy'
NA_VALUES = ['', 'n/a', 'NA', 'NaN', 'nan']
TASKS = ['mid', 'nback', 'sst']
MANUFACT = 'mr_y_adm__info__dev_manufact'

#############################################################################
# input columns and their tables
#############################################################################

# ROI value of each flag (not-missing test): (table, column)
ROI_COLUMNS = OrderedDict([
    ('smri__t1', ('mr_y_smri__t1__aseg', 'mr_y_smri__t1__aseg__cwm__lh_mean')),
    ('smri__t2', ('mr_y_smri__t2__aseg', 'mr_y_smri__t2__aseg__cwm__lh_mean')),
    ('dmri', ('mr_y_rsi__rnd__at', 'mr_y_rsi__rnd__at_wmean')),
    ('rsfmri', ('mr_y_rsfmri__corr__gpnet',
                'mr_y_rsfmri__corr__gpnet__def__sal_mean')),
    ('tfmri__mid', ('mr_y_tfmri__mid__arvn__aseg',
                    'mr_y_tfmri__mid__arvn__aseg__cwm__lh_beta')),
    ('tfmri__nback', ('mr_y_tfmri__nback__0b__aseg',
                      'mr_y_tfmri__nback__0b__aseg__cbwm__lh_beta')),
    ('tfmri__sst', ('mr_y_tfmri__sst__cg__aseg',
                    'mr_y_tfmri__sst__cg__aseg__cwm__lh_beta')),
])

def _input_tables():
    """Return OrderedDict column -> table; a table ending in LEGACY_SUFFIX
    is read from legacy_dir."""
    tabs = OrderedDict()
    tabs[MANUFACT] = ROW_TABLE
    raw = [('t1', 'smri__t1'), ('t2', 'smri__t2'), ('dmri', 'dmri'),
           ('rsfmri', 'rsfmri')] + [(t, 'tfmri__' + t) for t in TASKS]
    for _, rtype in raw:
        table = 'mr_y_qc__raw__' + rtype
        for suffix in ['series_count', 'ok_count']:
            tabs['{}__{}'.format(table, suffix)] = table
    tabs['mr_y_qc__raw__dmri__frame__ok_mean'] = 'mr_y_qc__raw__dmri'
    for mtype in ['fsurf', 't2', 'dmri', 'fmri']:
        table = 'mr_y_qc__post__man__' + mtype
        tabs[table + '_score'] = table
    for mtype in ['dmri', 'fmri']:
        tabs['mr_y_qc__post__aut__{}__b0__unwarp_indicator'.format(mtype)] = \
            'mr_y_qc__post__aut'
    tabs['mr_y_qc__mot__rsfmri__frame__censor_count'] = 'mr_y_qc__mot'
    for task in TASKS:
        table = 'mr_y_qc__raw__tfmri__' + task
        tabs[table + '__eprime__match_indicator'] = table
        tabs[table + '__eprime__scan__run__valid_count'] = table
        tabs[table + '__eprime__tdiff__ign_indicator'] = table + LEGACY_SUFFIX
        table = 'mr_y_tfmri__{}__beh'.format(task)
        tabs[table + '__qc_indicator'] = table
        tabs[table + '__start__delay_min'] = table
        tabs[table + '__start__delay_max'] = table
        tabs['mr_y_qc__mot__tfmri__{}_dof'.format(task)] = 'mr_y_qc__mot'
    tabs['mr_y_tfmri__mid__beh__trial_count'] = 'mr_y_tfmri__mid__beh'
    tabs['mr_y_tfmri__sst__beh__coderr_indicator'] = 'mr_y_tfmri__sst__beh'
    for table, col in ROI_COLUMNS.values():
        tabs[col] = table
    return tabs

#############################################################################
# criteria
#############################################################################

# tests on the input columns (pandas Series, numeric; manufacturer text)
def gt0(v):
    return v > 0

def eq1(v):
    return v == 1

def ne0(v):
    return v != 0

def notna(v):
    return v.notna()

def dmri_reps(nreps, nok, manufact):
    return (nreps >= 103) | ((manufact == 'Philips') & (nok >= 2) &
                             (nreps == 51))

def eprime_match(match, ign):
    return (match == 1) | (ign == 1)

def ge_delay(manufact, dmin, dmax):
    # NaN != 6.4, so a missing delay passes, as in release 7.0
    return (manufact != 'GE') | ((dmin != 6.4) & (dmax != 6.4))

def _crit(name, cols, test, miss_pass=False):
    """A criterion: name, input columns, test(*columns) -> bool Series;
    miss_pass: a missing value (single input) passes, else fails."""
    return OrderedDict([('name', name), ('cols', cols), ('test', test),
                        ('miss_pass', miss_pass)])


def _raw(rtype, suffix):
    return 'mr_y_qc__raw__{}__{}'.format(rtype, suffix)


def _common_t1(stype, rtype):
    return [_crit(stype + '_series', [_raw(rtype, 'series_count')], gt0),
            _crit(stype + '_ok', [_raw(rtype, 'ok_count')], gt0),
            _crit('t1_series', [_raw('smri__t1', 'series_count')], gt0),
            _crit('t1_ok', [_raw('smri__t1', 'ok_count')], gt0)]


def _man(mtype):
    return _crit(mtype + '_man', ['mr_y_qc__post__man__{}_score'.format(mtype)],
                 ne0, miss_pass=True)


def _roi(flag, stype):
    return _crit(stype + '_roi', [ROI_COLUMNS[flag][1]], notna)


def _tfmri(task):
    rtype = 'tfmri__' + task
    beh = 'mr_y_tfmri__{}__beh__'.format(task)
    crits = _common_t1(task, rtype) + [
        _crit('fmri_b0unwarp',
              ['mr_y_qc__post__aut__fmri__b0__unwarp_indicator'], eq1),
        _crit(task + '_eprime_match',
              [_raw(rtype, 'eprime__match_indicator'),
               _raw(rtype, 'eprime__tdiff__ign_indicator')], eprime_match),
        _crit(task + '_eprime_valid',
              [_raw(rtype, 'eprime__scan__run__valid_count')], gt0),
        _crit(task + '_ge_delay', [MANUFACT, beh + 'start__delay_min',
                                   beh + 'start__delay_max'], ge_delay)]
    if task == 'sst':
        crits.append(_crit('sst_coderr', [beh + 'coderr_indicator'],
                           lambda v: v == 0, miss_pass=True))
    crits += [_crit(task + '_beh_qc', [beh + 'qc_indicator'], eq1),
              _crit(task + '_dof',
                    ['mr_y_qc__mot__tfmri__{}_dof'.format(task)],
                    lambda v: v > 200)]
    if task == 'mid':
        crits.append(_crit('mid_trials', [beh + 'trial_count'],
                           lambda v: v == 100))
    crits += [_man('fsurf'), _man('fmri'), _roi('tfmri__' + task, task)]
    return crits


# per flag, the criteria in the order of the exclusion sequence (the order of
# the release 7.0 per-flag breakdown scripts)
CRITERIA = OrderedDict([
    ('smri__t1', [
        _crit('t1_series', [_raw('smri__t1', 'series_count')], gt0),
        _crit('t1_ok', [_raw('smri__t1', 'ok_count')], gt0),
        _man('fsurf'),
        _roi('smri__t1', 't1')]),
    ('smri__t2', _common_t1('t2', 'smri__t2') + [
        _man('fsurf'),
        _man('t2'),
        _roi('smri__t2', 't2')]),
    ('dmri', _common_t1('dmri', 'dmri') + [
        _crit('dmri_b0unwarp',
              ['mr_y_qc__post__aut__dmri__b0__unwarp_indicator'], eq1),
        _crit('dmri_reps', ['mr_y_qc__raw__dmri__frame__ok_mean',
                            _raw('dmri', 'ok_count'), MANUFACT], dmri_reps),
        _man('fsurf'),
        _man('dmri'),
        _roi('dmri', 'dmri')]),
    ('rsfmri', _common_t1('rsfmri', 'rsfmri') + [
        _crit('fmri_b0unwarp',
              ['mr_y_qc__post__aut__fmri__b0__unwarp_indicator'], eq1),
        _crit('rsfmri_frames', ['mr_y_qc__mot__rsfmri__frame__censor_count'],
              lambda v: v > 375),
        _man('fsurf'),
        _man('fmri'),
        _roi('rsfmri', 'rsfmri')]),
    ('tfmri__mid', _tfmri('mid')),
    ('tfmri__nback', _tfmri('nback')),
    ('tfmri__sst', _tfmri('sst')),
])

INPUT_TABLES = _input_tables()

#############################################################################
# helpers
#############################################################################

def _msg(verbose, text):
    if verbose:
        print('{}: {}'.format(APPNAME, text))


def flag_column(flag):
    return 'mr_y_qc__incl__{}_indicator'.format(flag)


def read_inputs(tabdir, legacy_dir):
    """Return (df, errmsgs): the input columns merged onto the ROW_TABLE
    rows; errmsgs lists each missing table or column."""
    cols = []
    for crits in CRITERIA.values():
        for crit in crits:
            cols += [c for c in crit['cols'] if c not in cols]
    bytable = OrderedDict([(ROW_TABLE, [])])
    errmsgs = []
    for col in cols:
        if col not in INPUT_TABLES:
            errmsgs.append('no input table for column {}'.format(col))
            continue
        bytable.setdefault(INPUT_TABLES[col], []).append(col)
    df = None
    for table, tcols in bytable.items():
        tdir = legacy_dir if table.endswith(LEGACY_SUFFIX) else tabdir
        fname = os.path.join(tdir, table + '.tsv')
        if not os.path.isfile(fname):
            errmsgs.append('missing table {}'.format(fname))
            continue
        header = pd.read_csv(fname, sep='\t', nrows=0).columns
        missing = [c for c in KEYS + tcols if c not in header]
        if missing:
            errmsgs += ['missing column {} in {}'.format(c, fname)
                        for c in missing]
            continue
        tdf = pd.read_csv(fname, sep='\t', usecols=KEYS + tcols, dtype=str,
                          na_values=NA_VALUES, keep_default_na=False)
        ndup = tdf.duplicated(KEYS).sum()
        if ndup:
            errmsgs.append('{} duplicate sessions in {}'.format(ndup, fname))
            continue
        for col in tcols:
            if col != MANUFACT:
                tdf[col] = pd.to_numeric(tdf[col], errors='coerce')
        if df is None:
            df = tdf
        else:
            df = df.merge(tdf, on=KEYS, how='left')
    return df, errmsgs


def apply_criteria(df):
    """Return the criteria table: KEYS, then per flag one bool column per
    criterion, <flag>__<criterion>, and the flag, <flag>__incl."""
    crit_df = df[KEYS].copy()
    for flag, crits in CRITERIA.items():
        incl = pd.Series(True, index=df.index)
        for crit in crits:
            vals = [df[c] for c in crit['cols']]
            res = crit['test'](*vals).fillna(False).astype(bool)
            if len(vals) == 1:
                res[vals[0].isna()] = crit['miss_pass']
            crit_df['{}__{}'.format(flag, crit['name'])] = res
            incl &= res
        crit_df['{}__incl'.format(flag)] = incl
    return crit_df


def write_report(fname, df, crit_df):
    with open(fname, 'w') as fid:
        fid.write('{} report: {} sessions\n'.format(APPNAME, len(crit_df)))
        fid.write('per criterion, sessions failing it: fail (any), seq (in '
                  'sequence, after passing every criterion above), only (no '
                  'other), miss (input missing)\n')
        for flag, crits in CRITERIA.items():
            fid.write('\n{}\n'.format(flag_column(flag)))
            fid.write('  {:<22s} {:>7s} {:>7s} {:>7s} {:>7s}  {}\n'.format(
                'criterion', 'fail', 'seq', 'only', 'miss', 'inputs'))
            names = ['{}__{}'.format(flag, c['name']) for c in crits]
            fails = ~crit_df[names]
            nfail = fails.sum(axis=1)
            remain = pd.Series(True, index=crit_df.index)
            for crit, name in zip(crits, names):
                seq = (remain & fails[name]).sum()
                remain &= ~fails[name]
                only = (fails[name] & (nfail == 1)).sum()
                miss = df[crit['cols']].isna().any(axis=1).sum()
                fid.write('  {:<22s} {:>7d} {:>7d} {:>7d} {:>7d}  {}\n'.format(
                    crit['name'], fails[name].sum(), seq, only, miss,
                    ' '.join(crit['cols'])))
            fid.write('  {:<22s} {:>7d}\n'.format(
                'included', crit_df['{}__incl'.format(flag)].sum()))

#############################################################################
# main function
#############################################################################

def set_imgincl_abcd(tabdir, outdir, legacy_dir=None, report_flag=0,
                     forceflag=0, verbose=1):
    """Write <outdir>/mr_y_qc__incl.tsv and its criteria table; with
    report_flag, imgincl_report.txt. Returns (fname, errcode, errmsg)."""
    fname = os.path.join(outdir, OUT_TABLE + '.tsv')
    fname_crit = os.path.join(outdir, OUT_TABLE + CRITERIA_SUFFIX + '.tsv')
    fname_report = os.path.join(outdir, REPORT_NAME)
    if legacy_dir is None:
        legacy_dir = tabdir
    if (not forceflag and os.path.isfile(fname) and os.path.isfile(fname_crit)
            and (not report_flag or os.path.isfile(fname_report))):
        _msg(verbose, 'output exists: {}'.format(fname))
        return fname, 0, ''
    df, errmsgs = read_inputs(tabdir, legacy_dir)
    if errmsgs:
        return None, 1, '; '.join(errmsgs)
    df = df.sort_values(KEYS).reset_index(drop=True)
    crit_df = apply_criteria(df)
    out_df = crit_df[KEYS].copy()
    for flag in CRITERIA:
        out_df[flag_column(flag)] = crit_df['{}__incl'.format(flag)]
    os.makedirs(outdir, exist_ok=True)
    out_df.astype({c: int for c in out_df.columns[2:]}).to_csv(
        fname, sep='\t', index=False)
    crit_df.astype({c: int for c in crit_df.columns[2:]}).to_csv(
        fname_crit, sep='\t', index=False)
    _msg(verbose, 'wrote {} ({} sessions)'.format(fname, len(out_df)))
    if report_flag:
        write_report(fname_report, df, crit_df)
        _msg(verbose, 'wrote {}'.format(fname_report))
    return fname, 0, ''

#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Set the ABCD imaging inclusion flags (mr_y_qc__incl) '
                    'from the compiled tabulated data.')
    parser.add_argument('--tabdir', required=True,
                        help='Directory of the compiled tables '
                             '(cmig_compile_bids_tabdat.py --outdir)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory')
    parser.add_argument('--legacy_dir', default=None,
                        help='Directory of the _desc-legacy tables '
                             '(abcd_write_legacy_tabdat.py; default: --tabdir)')
    parser.add_argument('--report_flag', type=int, default=0, choices=[0, 1],
                        help='Write imgincl_report.txt, the exclusion counts '
                             'per criterion (default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing output (default: 0)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1],
                        help='Display status messages (default: 1)')
    return parser.parse_args()


def run(args):
    fname, errcode, errmsg = set_imgincl_abcd(
        args.tabdir, args.outdir, legacy_dir=args.legacy_dir,
        report_flag=args.report_flag, forceflag=args.forceflag,
        verbose=args.verbose)
    if errcode:
        print('{}: ERROR: {}'.format(APPNAME, errmsg), file=sys.stderr)
    return errcode, errmsg


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
