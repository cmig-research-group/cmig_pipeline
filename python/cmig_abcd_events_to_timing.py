#!/usr/bin/env python
"""
cmig_abcd_events_to_timing.py

Public functions
----------------
events_to_timing(fname_events_tsv, TR, numTRs=500, minfrac=0.5,
                 outdir='.', outstem=None, skipTRs=None, block_gap=3.0,
                 task=None, forceflag=False, verbose=True)

  Reads a single BIDS _events.tsv and writes per-condition FSL/AFNI timing
  files for that run. Does not access the original E-Prime file or JSON
  sidecar. Returns (conditions_written, errcode, errmsg). The per-run timing
  QC columns of the events files are summarized by events_to_behavioral() in
  cmig_abcd_events_to_behavioral.py.

  outstem is the full per-run stem (including any BIDS run-NN entity, e.g.
  sub-XXXXXXXX_ses-XXA_task-<task>_run-NN); no separate scan number is
  appended. Per condition, six file types are written (content matches
  abcd_extract_eprime_*.m's own per-file semantics exactly, including two
  file-vs-name pairings that read as swapped at first glance -- see below):
    <outstem>_<cond>.1D             TR-aligned binary mask (AFNI)
    <outstem>_<cond>.txt            every active TR's onset time in seconds,
                                     TR-snapped, one per registered TR --
                                     NOT grouped into contiguous blocks,
                                     despite being separate from _block.txt
    <outstem>_<cond>_block.txt      raw per-event onset time (sub-TR
                                     precision), one per real event -- NOT
                                     deduplicated/grouped despite the name
    <outstem>_<cond>_block_dur.txt  mean of the raw per-event durations
                                     (not derived from TR-mask blocks)
    <outstem>_<cond>_fsl.txt        onset duration 1 rows (FSL 3-column)
    <outstem>_<cond>_rt.txt         response times in seconds

  A condition in the task's fixed list (_TASK_CONDS) with zero events
  this run still gets all six files, as MATLAB-style placeholders
  (needed for 3dDeconvolve): see _write_placeholder_files.

  NBACK exception: every condition except 'cue' collapses _block.txt/
  _block_dur.txt to blocks instead of the per-event list/mean above -- a
  block runs from an event's onset through every following event that starts
  within block_gap seconds of the previous one's offset; _block.txt has each
  block's first onset, _block_dur.txt the mean block span (first onset to
  last offset). With one block per condition per run, as the task presents
  it, this matches abcd_extract_eprime_nback.m's write_files() exactly; a
  condition split by a longer gap gives two blocks where MATLAB gives one
  spanning the gap. 'cue' and all MID/SST conditions use the generic
  per-event behavior.

  The task (MID, SST or NBACK) is read from the file name (task-<label>);
  task, when given, overrides it.

  SST: trial_type ssd_stop, the name in events files written before the
  08/14/2026 rename, is read as early_stop.

  The _events.tsv onsets are from the first stored volume, as BIDS
  requires, dummy volumes included; the dummy row's duration gives skipTRs.
  skipTRs, when given, overrides it, so the timing files match BOLD data
  with that many volumes dropped.
  .1D/.txt are masked on that full grid, then lose the first skipTRs TRs;
  _block.txt and _fsl.txt give onsets from the first non-dummy TR (the
  dummy duration subtracted), matching abcd_extract_eprime_*.m's
  rel_onset/rel_offset, for BOLD data with its skipTRs dummy volumes
  dropped; _block_dur.txt and _rt.txt do not depend on the zero point.

Created:  05/18/2026 by Don Hagler
Prev mod: 06/25/2026 by Don Hagler
Prev mod: 06/29/2026 by Don Hagler
Prev mod: 07/07/2026 by Don Hagler
Prev mod: 07/08/2026 by Don Hagler
Prev mod: 07/09/2026 by Don Hagler
Prev mod: 07/13/2026 by Don Hagler
Prev mod: 07/27/2026 by Don Hagler
Prev mod: 07/28/2026 by Don Hagler
Prev mod: 09/18/2026 by Don Hagler
Prev mod: 10/01/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/21/2026 by Don Hagler
Last mod: 11/01/2026 by Don Hagler

Usage
-----
  cmig_abcd_events_to_timing.py --fname_events run-01_events.tsv --TR 0.8 \
      [--numTRs 500] [--minfrac 0.5] [--skipTRs N] [--block_gap 3] [--task T] \
      [--outdir .] [--outstem STEM] [--forceflag 0|1] [--verbose 0|1]
"""
from __future__ import print_function
import argparse
import math
import os
import re
import sys

import numpy as np
import pandas as pd

from cmig_select_bids_func_runs import has_responses

# SST-specific mindur overrides (abcd_extract_eprime_sst.m, present since the
# function's original 2017 creation): go-type conditions are RT-terminated
# and can be shorter than the generic minfrac*TR margin tolerates; stop-type
# conditions are a fixed ~20ms marker, far too brief to ever clear that
# margin at all. Both would otherwise silently fail to register in any TR.
# early_stop is grouped with the go-type conditions since it is likewise
# response-defined and typically short. MID/NBACK have no such override --
# both always use the generic minfrac*TR.
# zero_ssd_response takes the go-type mindur, as early_stop does in
# abcd_extract_eprime_sst.m; its _rt.txt is from response_time (its ssd_rt)
_SST_MINDUR_GO = 0.08
_SST_MINDUR_STOP = 0.0
_SST_GO_CONDS = {'correct_go', 'correctlate_go', 'incorrect_go',
                'incorrectlate_go', 'noresp_go', 'early_stop'}
_SST_STOP_CONDS = {'correct_stop', 'incorrect_stop'}
_SST_ZR_COND = 'zero_ssd_response'

# per-task fixed condition lists, for writing MATLAB-style placeholder
# files (see _write_placeholder_files) when a condition has zero events in
# a run -- 3dDeconvolve needs one timing file set per condition across
# every run, not just the conditions that happened to occur. Matches each
# task's standard (default-flag) condition set.
_MID_CONDS = tuple('{}_{}'.format(c, p) for c in
                   ('small_reward', 'large_reward', 'small_loss', 'large_loss', 'neutral')
                   for p in ('antic', 'pos_feedback', 'neg_feedback'))
_SST_CONDS = tuple(sorted(_SST_GO_CONDS | _SST_STOP_CONDS | {_SST_ZR_COND}))
_NBACK_CONDS = tuple('{}_{}'.format(b, s) for b in ('2_back', '0_back')
                     for s in ('posface', 'neutface', 'negface', 'place')) + ('cue',)
_TASK_CONDS = {'MID': _MID_CONDS, 'SST': _SST_CONDS, 'NBACK': _NBACK_CONDS}


# ── TR-grid utilities (self-contained; mirror cmig_abcd_eprime_to_events) ──

def build_tr_grid(numTRs, TR):
    """Return (TR_onset, TR_offset) arrays for numTRs analysis TRs."""
    TR_onset = np.arange(numTRs, dtype=float) * TR
    TR_offset = TR_onset + TR
    return TR_onset, TR_offset


def compute_tr_mask(rel_onset, rel_offset, TR_onset, TR_offset, mindur):
    """Boolean array: True for each TR where the event registers.

    Event registers in TR k when:
      rel_onset < TR_offset[k] - mindur  AND  rel_offset > TR_onset[k] + mindur
    """
    return (rel_onset < TR_offset - mindur) & (rel_offset > TR_onset + mindur)


# ── private helpers ────────────────────────────────────────────────────────

def _to_float(val):
    """Convert a TSV cell to float; return NaN for 'n/a', None, or bad strings."""
    if val is None:
        return np.nan
    s = str(val).strip().lower()
    if s in ('n/a', 'na', 'nan', 'none', ''):
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def _derive_outstem(fname_tsv):
    """Derive output stem from a TSV filename, stripping only the _events.tsv
    (or .tsv) suffix -- the run identity (_run-NN, _desc-invalid, etc.) is
    part of the input filename's stem and is kept as-is, since output timing
    filenames are built as <outstem>_<cond>... with no separate scan number
    reappended.
    """
    base = os.path.basename(fname_tsv)
    base = re.sub(r'_events\.tsv$', '', base)
    base = re.sub(r'\.tsv$', '', base)
    return base


def _find_blocks(tr_mask, TR):
    """Extract contiguous active blocks from a binary TR mask.

    Returns list of (onset_s, duration_s) for each run of consecutive 1s.
    onset_s is the start time of the first active TR in the block.
    """
    blocks = []
    n = len(tr_mask)
    k = 0
    while k < n:
        if tr_mask[k]:
            start = k
            while k < n and tr_mask[k]:
                k += 1
            blocks.append((float(start) * TR, float(k - start) * TR))
        else:
            k += 1
    return blocks


def _safe_write(fname, content, forceflag, verbose):
    """Write text file; skip if exists and forceflag is False."""
    if not forceflag and os.path.exists(fname):
        if verbose:
            print('  exists, skipping: {}'.format(fname))
        return
    with open(fname, 'w') as fh:
        fh.write(content)
    if verbose:
        print('  wrote {}'.format(fname))


# ── per-condition timing writers ───────────────────────────────────────────

def _write_1d(outstem, cond, tr_mask, forceflag, verbose):
    """Write TR-aligned binary mask (.1D), one value per line."""
    fname = '{}_{}.1D'.format(outstem, cond)
    content = '\n'.join(str(int(v)) for v in tr_mask) + '\n'
    _safe_write(fname, content, forceflag, verbose)


def _write_txt(outstem, cond, active_tr_times, forceflag, verbose):
    """Write AFNI local_times file (.txt): every active TR's onset time in
    seconds (TR-snapped, one per registered TR, not grouped into blocks --
    matches abcd_extract_eprime_sst.m's event_TRs*TR, derived from the TR
    mask directly), all on one line, '*' if none.
    """
    fname = '{}_{}.txt'.format(outstem, cond)
    if len(active_tr_times) == 0:
        content = '*\n'
    else:
        content = ' '.join('{:.3f}'.format(t) for t in active_tr_times) + '\n'
    _safe_write(fname, content, forceflag, verbose)


def _write_block_txt(outstem, cond, onsets_s, forceflag, verbose):
    """Write block onset times (_block.txt): raw per-event onset (rel_onset),
    one per real event -- matches MATLAB's own _block.txt despite the name
    (not grouped/deduplicated into contiguous blocks; every event gets its
    own line, sub-TR precision)."""
    fname = '{}_{}_{}'.format(outstem, cond, 'block.txt')
    if len(onsets_s) == 0:
        content = '*\n'
    else:
        content = '\n'.join('{:.3f}'.format(t) for t in onsets_s) + '\n'
    _safe_write(fname, content, forceflag, verbose)


def _write_block_dur_txt(outstem, cond, durations_s, forceflag, verbose):
    """Write mean event duration (_block_dur.txt): mean of the raw per-event
    durations (rel_offset - rel_onset) -- matches MATLAB's
    mean(rel_offset-rel_onset), not derived from contiguous TR-mask blocks."""
    fname = '{}_{}_{}'.format(outstem, cond, 'block_dur.txt')
    finite = [d for d in durations_s if not math.isnan(d)]
    if not finite:
        content = '0\n'
    else:
        content = '{:.3f}\n'.format(float(np.mean(finite)))
    _safe_write(fname, content, forceflag, verbose)


def _write_fsl_txt(outstem, cond, onsets_s, durations_s, forceflag, verbose):
    """Write FSL 3-column timing file (_fsl.txt): onset duration 1."""
    fname = '{}_{}_{}'.format(outstem, cond, 'fsl.txt')
    lines = []
    for o, d in zip(onsets_s, durations_s):
        dur = d if not math.isnan(d) else 0.0
        lines.append('{:.3f} {:.3f} 1'.format(o, dur))
    content = '\n'.join(lines) + '\n' if lines else '*\n'
    _safe_write(fname, content, forceflag, verbose)


def _write_rt_txt(outstem, cond, rts_s, forceflag, verbose):
    """Write response times (_rt.txt): one value per row; 'n/a' for missing."""
    fname = '{}_{}_{}'.format(outstem, cond, 'rt.txt')
    lines = ['{:.3f}'.format(rt) if not math.isnan(rt) else 'n/a' for rt in rts_s]
    content = '\n'.join(lines) + '\n' if lines else '*\n'
    _safe_write(fname, content, forceflag, verbose)


def _write_placeholder_files(outstem, cond, numTRs, forceflag, verbose):
    """Write all six files for a condition with zero events this run,
    matching abcd_extract_eprime_*.m's write_files() empty-condition
    branch: numTRs zeros (full grid, not the post-skip analysis length)
    for .1D, '*' for .txt/_block.txt, '0' for _block_dur.txt, and a
    truly empty file for _fsl.txt/_rt.txt."""
    _safe_write('{}_{}.1D'.format(outstem, cond), '0\n' * numTRs, forceflag, verbose)
    _safe_write('{}_{}.txt'.format(outstem, cond), '*\n', forceflag, verbose)
    _safe_write('{}_{}_{}'.format(outstem, cond, 'block.txt'), '*\n', forceflag, verbose)
    _safe_write('{}_{}_{}'.format(outstem, cond, 'block_dur.txt'), '0\n', forceflag, verbose)
    _safe_write('{}_{}_{}'.format(outstem, cond, 'fsl.txt'), '', forceflag, verbose)
    _safe_write('{}_{}_{}'.format(outstem, cond, 'rt.txt'), '', forceflag, verbose)


# ── public function ────────────────────────────────────────────────────────

def events_to_timing(fname_events_tsv, TR, numTRs=500, minfrac=0.5,
                     outdir='.', outstem=None, skipTRs=None, block_gap=3.0,
                     task=None, forceflag=False, verbose=True):
    """Read a BIDS events TSV and write per-condition FSL/AFNI timing files.

    Parameters
    ----------
    fname_events_tsv : str
        Path to a single _events.tsv file (one scan).
    TR : float
        Repetition time in seconds.
    numTRs : int
        Total number of TRs in the run including any leading dummies (default 500).
    minfrac : float
        Minimum fraction of a TR that an event must span to register (default 0.5).
        For SST specifically (detected from the file name's task-<label>),
        this is overridden per condition to match abcd_extract_eprime_sst.m's
        long-standing (since 2017) mindur convention: 0.0 for stop-type
        conditions (correct_stop/incorrect_stop -- a fixed ~20ms
        marker, far too brief to clear any nonzero margin) and 0.08s for
        go-type conditions (correct_go/correctlate_go/incorrect_go/
        incorrectlate_go/noresp_go/early_stop -- RT-terminated, and often
        shorter than the generic minfrac*TR margin tolerates). MID/NBACK
        have no such override and always use minfrac*TR.

        _rt.txt has a similar SST override: go-type conditions read
        response_time_go_adjusted, stop-type read
        response_time_stop_adjusted, instead of response_time, matching
        abcd_extract_eprime_sst.m.
    outdir : str
        Output directory (created if absent).
    outstem : str or None
        Base output stem; derived from fname_events_tsv if None.
    skipTRs : int or None
        Dummy volumes dropped from the BOLD data; None takes them from the
        dummy row's duration (default None).
    block_gap : float
        NBACK: seconds from one event's offset to the next event's onset
        beyond which a condition starts a new block (default 3.0).
    task : str or None
        MID, SST or NBACK (any case); None reads the file name's task-<label>
        (default None).
    forceflag : bool
        Overwrite existing output files (default False).
    verbose : bool
        Print progress (default True).

    Returns
    -------
    conditions_written : list of str
    errcode : int    0=success 1=error
    errmsg : str
    """
    if not os.path.exists(fname_events_tsv):
        return [], 1, 'events_to_timing: file not found: {}'.format(fname_events_tsv)

    os.makedirs(outdir, exist_ok=True)

    if outstem is None:
        outstem = _derive_outstem(fname_events_tsv)

    outstem_path = os.path.join(outdir, outstem)

    # ── read TSV ───────────────────────────────────────────────────────────
    try:
        df = pd.read_csv(fname_events_tsv, sep='\t', dtype=str, low_memory=False)
    except Exception as exc:
        return [], 1, 'events_to_timing: cannot read {}: {}'.format(fname_events_tsv, exc)

    required = {'onset', 'duration', 'trial_type'}
    missing = required - set(df.columns)
    if missing:
        return [], 1, 'events_to_timing: missing columns {} in {}'.format(
            sorted(missing), fname_events_tsv)

    task = (task or '').strip().upper()
    if not task:
        m = re.search(r'_task-([A-Za-z0-9]+)_', os.path.basename(fname_events_tsv))
        if m:
            task = m.group(1).upper()

    # events files written before the 08/14/2026 rename say ssd_stop
    if task == 'SST':
        df.loc[df['trial_type'].str.strip() == 'ssd_stop', 'trial_type'] = 'early_stop'

    # ── infer skip TRs from dummy row ──────────────────────────────────────
    dummy_rows = df[df['trial_type'].str.strip().str.lower() == 'dummy']
    if skipTRs is not None:
        skip_time_s = int(skipTRs) * TR
    elif len(dummy_rows) > 0:
        skip_time_s = _to_float(dummy_rows['duration'].iloc[0])
        if math.isnan(skip_time_s):
            skip_time_s = 0.0
    else:
        skip_time_s = 0.0

    nskipTRs = int(round(skip_time_s / TR)) if TR > 0 else 0
    n_analysis_TRs = max(1, numTRs - nskipTRs)

    default_mindur = minfrac * TR
    TR_onset, TR_offset = build_tr_grid(n_analysis_TRs, TR)
    # full, pre-skip grid: masking happens against raw (unshifted) onset here,
    # matching MATLAB's write_files() bit-for-bit -- shifting onset instead
    # (subtracting skip_time_s before comparing) adds a rounding step MATLAB
    # never performs, which can flip an exact-tie mindur boundary comparison.
    TR_onset_full, TR_offset_full = build_tr_grid(numTRs, TR)

    if verbose:
        print('events_to_timing: outstem={} skipTRs={} analysisTRs={} TR={} minfrac={}'.format(
            outstem, nskipTRs, n_analysis_TRs, TR, minfrac))

    # ── collect conditions (all trial_type values except special rows) ─────
    _skip_types = {'dummy', 'scanner_trigger', 'task_init'}
    all_types = df['trial_type'].str.strip().dropna().unique().tolist()
    observed_conditions = sorted(t for t in all_types
                                 if t.lower() not in _skip_types and t != '')

    if not observed_conditions:
        return [], 1, 'events_to_timing: no non-dummy trial types found in {}'.format(
            fname_events_tsv)

    # union in this task's fixed condition list (_TASK_CONDS) so a
    # condition with zero events this run still gets a placeholder file
    # set, matching MATLAB
    fixed_conds = _TASK_CONDS.get(task, ())
    conditions = sorted(set(observed_conditions) | set(fixed_conds))

    has_rt = 'response_time' in df.columns

    # a run with zero recorded responses has degenerate trial-type
    # classifications (e.g. SST's correct_go/incorrect_go/noresp_go are
    # themselves defined by response outcome), so timing files built from
    # it would misrepresent a non-responsive session as real engagement.
    # Skip them for this run; not an error -- F1/F3 still produce output.
    # SST writes no response as rt 0, so for SST the trial types decide
    #   (has_responses(), as the run selection uses)
    if (has_rt or task == 'SST') and not has_responses(
            df.to_dict('records'), task):
        print('events_to_timing: WARNING: no responses recorded in {} '
              '-- timing files not generated'.format(fname_events_tsv))
        return [], 0, ''

    conditions_written = []

    for cond in conditions:
        if task == 'SST' and cond in _SST_STOP_CONDS:
            mindur = _SST_MINDUR_STOP
        elif task == 'SST' and (cond in _SST_GO_CONDS or cond == _SST_ZR_COND):
            mindur = _SST_MINDUR_GO
        else:
            mindur = default_mindur

        cond_rows = df[df['trial_type'].str.strip() == cond].copy()

        onsets_raw = [_to_float(v) for v in cond_rows['onset']]
        durs_raw = [_to_float(v) for v in cond_rows['duration']]
        # SST rt.txt uses response_time_go_adjusted/response_time_stop_adjusted,
        # not raw response_time -- matches MATLAB's write_files()
        if task == 'SST' and cond in _SST_GO_CONDS and 'response_time_go_adjusted' in df.columns:
            rt_col = 'response_time_go_adjusted'
        elif task == 'SST' and cond in _SST_STOP_CONDS and 'response_time_stop_adjusted' in df.columns:
            rt_col = 'response_time_stop_adjusted'
        elif has_rt:
            rt_col = 'response_time'
        else:
            rt_col = None
        rts_raw = ([_to_float(v) for v in cond_rows[rt_col]]
                   if rt_col else [np.nan] * len(cond_rows))

        # keep only rows with valid onset
        valid_rows = [(o, d, r) for o, d, r in zip(onsets_raw, durs_raw, rts_raw)
                      if not math.isnan(o)]

        cond_safe = re.sub(r'[^A-Za-z0-9_\-]', '_', cond)

        if not valid_rows:
            if cond in fixed_conds:
                _write_placeholder_files(outstem_path, cond_safe, numTRs, forceflag, verbose)
                conditions_written.append(cond)
            elif verbose:
                print('events_to_timing: no valid onsets for condition {!r} — skipping'.format(cond))
            continue

        # onsets in the events.tsv are from the first stored volume, as BIDS
        # requires, dummy volumes included (old-style files too). The .1D and
        # .txt are masked on that full grid and then lose the first nskipTRs
        # TRs; _block.txt and _fsl.txt instead give onsets from the first
        # non-dummy TR, as MATLAB's rel_onset in write_files() does, for BOLD
        # data with its skipTRs dummy volumes dropped, so the dummy duration
        # is subtracted for them (in integer ms, as below)
        onsets_s = np.array([row[0] for row in valid_rows])
        skip_ms = round(skip_time_s * 1000)
        rel_onsets_s = np.array([(round(o * 1000) - skip_ms) / 1000.0
                                 for o in onsets_s])
        durations_s = np.array([row[1] for row in valid_rows])
        rts_s = np.array([row[2] for row in valid_rows])

        # ── build TR binary mask ───────────────────────────────────────────
        # Masked against the full, pre-skip grid using raw onset (matching
        # MATLAB's write_files() exactly), then sliced -- not shifted onset
        # against a small grid, which adds a rounding step MATLAB doesn't
        # have (see TR_onset_full comment above).
        tr_mask_full = np.zeros(numTRs, dtype=int)
        for onset, dur in zip(onsets_s, durations_s):
            if math.isnan(dur):
                # unknown duration (not a real zero) -- no better estimate
                dur = TR
            elif dur < 0:
                dur = 0.0
            # Sum in integer milliseconds, then divide once -- matching
            # MATLAB's rel_offset = (offset_ms - event_ref_time_ms)/1000,
            # a single division of already-summed integer milliseconds.
            # onset/dur are always exact multiples of 0.001s (rounded to
            # milliseconds when Function 1 wrote them), so round(x*1000)
            # exactly recovers the original integer. A plain onset+dur
            # float addition sums two independently-rounded seconds values
            # and can differ from MATLAB's result in the last bit -- enough
            # to flip a mindur=0 exact-tie comparison (e.g. an SST stop
            # event landing exactly on a TR boundary).
            offset = (round(onset * 1000) + round(dur * 1000)) / 1000.0
            event_mask = compute_tr_mask(onset, offset, TR_onset_full, TR_offset_full, mindur)
            tr_mask_full |= event_mask.astype(int)

        tr_mask = tr_mask_full[nskipTRs:]

        # every active TR's own onset time, ungrouped -- matches MATLAB's
        # .txt (event_TRs*TR, derived directly from the TR mask, not
        # collapsed into contiguous blocks).
        active_tr_times = TR_onset[tr_mask.astype(bool)]

        # blocks are computed only for the verbose diagnostic below (not
        # used by any output file -- MATLAB's own _block.txt/.txt do not
        # group into contiguous blocks either, despite the name).
        blocks = _find_blocks(tr_mask, TR)

        if verbose:
            print('events_to_timing: condition {!r}: {} events, {} active TRs, {} blocks, mindur={}'.format(
                cond, len(onsets_s), int(tr_mask.sum()), len(blocks), mindur))

        # NBACK-specific _block.txt/_block_dur.txt collapsing (Function 2):
        # abcd_extract_eprime_nback.m's write_files() collapses every non-cue
        # stim-type condition to a single block value -- _block.txt = the
        # first event's own onset only; _block_dur.txt = the last event's
        # own offset minus the first event's own onset (the full block
        # span), not a per-event list/mean. Confirmed against real output
        # (0_back_negface_block.txt='31.72', _block_dur.txt='24.56'). 'cue'
        # is NOT collapsed this way (MATLAB writes the full per-event onset
        # list and the mean per-event duration for it, same as the generic
        # per-event default below), and MID/SST are unaffected.
        # A condition split by a gap over block_gap gives one block per
        # part, rather than one block spanning the gap.
        if task == 'NBACK' and cond != 'cue':
            order = np.argsort(onsets_s)
            ord_onsets = onsets_s[order]
            ord_rel = rel_onsets_s[order]
            ord_offsets = ord_onsets + np.array(
                [d if not math.isnan(d) else 0.0 for d in durations_s[order]])
            starts = [0] + [i for i in range(1, len(ord_onsets))
                            if ord_onsets[i] - ord_offsets[i - 1] > block_gap]
            ends = starts[1:] + [len(ord_onsets)]
            block_onsets_s = np.array([ord_rel[i] for i in starts])
            block_durations_s = np.array([ord_offsets[j - 1] - ord_onsets[i]
                                          for i, j in zip(starts, ends)])
        else:
            block_onsets_s = rel_onsets_s
            block_durations_s = durations_s

        # ── write six timing files ─────────────────────────────────────────
        _write_1d(outstem_path, cond_safe, tr_mask, forceflag, verbose)
        _write_txt(outstem_path, cond_safe, active_tr_times, forceflag, verbose)
        _write_block_txt(outstem_path, cond_safe, block_onsets_s, forceflag, verbose)
        _write_block_dur_txt(outstem_path, cond_safe, block_durations_s, forceflag, verbose)
        _write_fsl_txt(outstem_path, cond_safe, rel_onsets_s, durations_s, forceflag, verbose)
        _write_rt_txt(outstem_path, cond_safe, rts_s, forceflag, verbose)

        conditions_written.append(cond)

    return conditions_written, 0, ''


def parse_args():
    parser = argparse.ArgumentParser(
        description='Per-condition FSL/AFNI timing files from one BIDS events TSV.')
    parser.add_argument('--fname_events', required=True)
    parser.add_argument('--TR', type=float, required=True)
    parser.add_argument('--numTRs', type=int, default=500,
                        help='TRs in the run, dummy volumes included')
    parser.add_argument('--minfrac', type=float, default=0.5)
    parser.add_argument('--skipTRs', type=int, default=None,
                        help='dummy volumes dropped (default: from the dummy row)')
    parser.add_argument('--block_gap', type=float, default=3.0)
    parser.add_argument('--task', default=None,
                        help='MID, SST or NBACK (default: task-<label> of the file name)')
    parser.add_argument('--outdir', default='.')
    parser.add_argument('--outstem', default=None)
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1])
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1])
    return parser.parse_args()


def main():
    args = parse_args()
    _, errcode, errmsg = events_to_timing(
        args.fname_events, args.TR, numTRs=args.numTRs, minfrac=args.minfrac,
        outdir=args.outdir, outstem=args.outstem, skipTRs=args.skipTRs,
        block_gap=args.block_gap, task=args.task, forceflag=bool(args.forceflag),
        verbose=bool(args.verbose))
    if errcode:
        print('cmig_abcd_events_to_timing: ERROR: {}'.format(errmsg),
              file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
