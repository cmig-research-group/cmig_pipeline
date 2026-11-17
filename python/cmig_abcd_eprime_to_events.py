#!/usr/bin/env python
"""
cmig_abcd_eprime_to_events.py

Reads an ABCD E-Prime behavioral file and writes BIDS-compatible output
under outdir (the caller-supplied func/ or beh/ directory):
  <outstem>_run-<NN>_events.tsv            (one per BIDS run number, caller-
                                            assigned via run_numbers; REC:
                                            <outstem>_events.tsv, no run
                                            entity, always one administration)
  <outstem>_run-<NN>_desc-invalid_events.tsv (a run with the wrong trial
                                            count or no responses, mapped
                                            to a real run-NN)
  outdir/eprime/<outstem>_events_switched.csv (non-standard; only when a
                                            button-swap was kept)

The column descriptions are not written per session: write_root_sidecars()
writes one task-<task>_events.json per task at the rawdata root, which applies
to every events file of that task by BIDS inheritance.

Public functions
----------------
eprime_to_events(fname_eprime, outdir='.', outstem=None,
                 task=None, TR=0.8, skipTRs=None,
                 Manufacturer=None, SoftwareVersion=None,
                 run_offset=1, switch_thresh=0.3, error_flag=False,
                 extra_flag=False, ste_recode_flag=None, fname_wm=None,
                 run_numbers=None, run_times=None, run_times_tolerance=4.0,
                 forceflag=False, verbose=True)
Returns (n_runs, runs_ok, errcode, errmsg).
write_root_sidecars(rootdir, tasks=None, error_flag=False, extra_flag=False,
                    forceflag=False, verbose=True)
Returns (errcode, errmsg).

Command line: write the root sidecars
  cmig_abcd_eprime_to_events.py --rootdir DIR [--task ...] [--forceflag 1]

Created:  05/15/2026 by Don Hagler
Prev mod: 06/25/2026 by Don Hagler
Prev mod: 06/29/2026 by Don Hagler
Prev mod: 07/07/2026 by Don Hagler
Prev mod: 07/08/2026 by Don Hagler
Prev mod: 07/09/2026 by Don Hagler
Prev mod: 07/10/2026 by Don Hagler
Prev mod: 07/13/2026 by Don Hagler
Prev mod: 07/24/2026 by Don Hagler
Prev mod: 07/27/2026 by Don Hagler
Prev mod: 07/28/2026 by Don Hagler
Prev mod: 07/29/2026 by Don Hagler
Prev mod: 07/30/2026 by Don Hagler
Prev mod: 07/31/2026 by Don Hagler
Prev mod: 08/03/2026 by Don Hagler
Prev mod: 08/14/2026 by Don Hagler
Prev mod: 10/01/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/21/2026 by Don Hagler
Prev mod: 10/31/2026 by Don Hagler
Last mod: 11/01/2026 by Don Hagler
"""
from __future__ import print_function
import codecs
import datetime
import argparse
import json
import math
import os
import re
import sys

import dateutil.parser
import numpy as np
import pandas as pd
from io import StringIO
from scipy.stats import norm as _norm

from cmig_set_abcd_skipTRs import set_abcd_skipTRs
from cmig_select_bids_func_runs import has_responses

# ── column name mappings per task ──────────────────────────────────────────

_MID_COLS = [
    'ExperimentName', 'ExperimentVersion', 'RuntimeVersion', 'Condition',
    'Cue.OnsetTime', 'Cue.Duration', 'prbrt', 'prbacc', 'prbresp',
    'Probe.OnsetTime', 'Probe.Duration', 'Probe.OffsetTime', 'Probe.RT', 'Probe.ACC',
    'Probe.RESP', 'ResponseCheck',
    'TextDisplay1.RTTime', 'TextDisplay1.RESP', 'Feedback.RTTime', 'Feedback.RESP',
    'Feedback.OnsetTime', 'Feedback.OffsetTime',
    'Anticipation.OnsetTime', 'Anticipation.Duration', 'RunMoney',
    'Trigger.RTTime', 'GetReady.RTTime', 'PrepTime.OnsetTime', 'Block',
    'SessionDate', 'SessionTime',
]
_MID_FIELDS = [
    'experiment_name', 'experiment_version', 'runtime_version', 'block',
    'cue_onset', 'cue_dur', 'cue_rt', 'cue_acc', 'cue_resp',
    'probe_onset', 'probe_dur', 'probe_offset', 'probe_rt', 'probe_acc',
    'probe_resp', 'response_check',
    'text_display_rttime', 'text_display_resp', 'feedback_rttime', 'feedback_resp',
    'feedback_onset', 'feedback_offset',
    'fix2_onset', 'fix2_dur', 'money',
    'trig_time', 'getready_rttime', 'initfix_onset', 'raw_block',
    'session_date', 'session_time',
]

_NBACK_COLS = [
    'ExperimentName', 'ExperimentVersion', 'RuntimeVersion', 'Procedure[Block]',
    'BlockType', 'StimType', 'TargetType', 'Stimulus',
    'Stim.OnsetTime', 'Stim.OffsetTime', 'Stim.Duration', 'Stim.RT', 'Stim.ACC', 'Stim.RESP',
    'CorrectResponse',
    'CueFix.OnsetTime', 'CueFix.Duration',
    'Cue2Back.OffsetTime', 'CueTarget.OffsetTime',
    'Trigger.RTTime', 'GetReady.RTTime', 'GetReady2.RTTime',
    'CueFix.StartTime', 'SessionDate', 'SessionTime',
]
_NBACK_FIELDS = [
    'experiment_name', 'experiment_version', 'runtime_version', 'block',
    'block_type', 'stim_type', 'target_type', 'stim',
    'stim_onset', 'stim_offset', 'stim_dur', 'stim_rt', 'stim_acc', 'stim_resp',
    'correct_response',
    'cuefix_onset', 'cuefix_dur',
    'cue2back_offset', 'cuetarget_offset',
    'trig_time', 'getready_rttime', 'getready2_rttime',
    'cuefix_start', 'session_date', 'session_time',
]

_SST_COLS = [
    'ExperimentName', 'ExperimentVersion', 'RuntimeVersion', 'Block', 'Trial', 'TrialCode',
    'Go.OnsetTime', 'Go.OffsetTime', 'Go.Duration', 'Go.RT', 'Go.RESP', 'Go.CRESP',
    'SSDDur', 'SSD.RT', 'SSD.RESP',
    'Fix.OnsetTime', 'Fix.Duration', 'Fix.RT', 'Fix.RESP',
    'StopSignal.OnsetTime', 'StopSignal.Duration', 'StopSignal.RT', 'StopSignal.RESP',
    'SSD.OnsetTime', 'StopSignal.StartTime',
    'GetReady.RTTime', 'Wait4Scanner.RTTime', 'Wait4ScannerB.RTTime',
    'BeginFix.OnsetTime', 'Trigger.RTTime', 'Stimulus',
    'BeginFix.StartTime', 'SessionDate', 'SessionTime',
]
_SST_FIELDS = [
    'experiment_name', 'experiment_version', 'runtime_version', 'block', 'trial_run', 'trial_code',
    'go_onset', 'go_offset', 'go_dur', 'go_rt', 'go_resp', 'go_cresp',
    'ssd_dur', 'ssd_rt', 'ssd_resp',
    'fix_onset', 'fix_dur', 'fix_rt', 'fix_resp',
    'stop_onset', 'stop_dur', 'stop_rt', 'stop_resp',
    'ssd_onset', 'stop_start',
    'getready_rttime', 'wait1_rttime', 'wait2_rttime',
    'beginfix_onset', 'trig_time', 'stim_file',
    'beginfix_start', 'session_date', 'session_time',
]

_WM_COLS = _NBACK_COLS
_WM_FIELDS = _NBACK_FIELDS

_REC_COLS = [
    'ExperimentName', 'ExperimentVersion', 'RuntimeVersion', 'Block',
    'Stim.ACC[Block]', 'Stim.CRESP[Block]', 'Stim.RESP[Block]',
    'Stim.RT[Block]', 'StimType[Block]', 'Stimulus[Block]',
]
_REC_FIELDS = [
    'experiment_name', 'experiment_version', 'runtime_version', 'block',
    'stim_acc_block', 'stim_cresp_block', 'stim_resp_block',
    'stim_rt_block', 'stim_type_block', 'stim_block',
]


# ── internal: read E-Prime spreadsheet ────────────────────────────────────

def read_eprime(fname, outstem, cols, fieldnames, verbose=True):
    """Read an E-Prime exported spreadsheet into a pandas DataFrame.

    Handles UTF-8/UTF-16 encoding, per-row quotes, Windows blank lines,
    tab/comma separator, and header row offsets 0..4.

    Parameters
    ----------
    fname : str
    outstem : str    path stem (unused; kept for call-site symmetry with
        process_mid/sst/nback/rec, which do write outstem-based files)
    cols : list of str    E-Prime column names to extract
    fieldnames : list of str    Python field names (same length as cols)
    verbose : bool

    Returns
    -------
    df : pd.DataFrame or None
    errcode : int    0=success 1=error
    errmsg : str
    available_fields : set of str
        Fieldnames (Python-side names, matching `fieldnames`) whose raw
        E-Prime column was genuinely found in the file's header, as opposed
        to backfilled with NaN because the raw column was absent -- distinct
        from checking the resulting column's content, which can't tell
        "column absent" from "column present but empty for this subject".
    """
    rows_n_min = 20
    rows_head_n_max = 5
    cols_n_min = 10

    # read raw text (try encodings in order)
    text = None
    for enc in ('utf-8-sig', 'utf-16', 'latin-1'):
        try:
            with codecs.open(fname, 'r', encoding=enc, errors='strict') as fh:
                text = fh.read()
            break
        except (UnicodeDecodeError, UnicodeError):
            continue

    if text is None:
        return None, 1, 'read_eprime: cannot decode {}'.format(fname), set()

    lines = text.splitlines()
    n_lines = len(lines)

    if n_lines < rows_n_min:
        return None, 1, 'read_eprime: too few lines ({}) in {}'.format(n_lines, fname), set()

    # detect quoted rows: check a sample around the midpoint
    mid_idx = min(int(math.floor(0.5 * (rows_n_min + n_lines))), n_lines - 1)
    check_idxs = range(max(0, mid_idx - 2), min(n_lines, mid_idx + 4))
    n_quoted = sum(1 for i in check_idxs
                   if lines[i].lstrip().startswith('"'))
    quoted_rows = n_quoted > max(1, len(list(check_idxs)) // 2)

    if quoted_rows:
        cleaned = []
        for line in lines:
            s = line.strip()
            if s.startswith('"') and s.endswith('"'):
                s = s[1:-1].replace('""', '"')
            cleaned.append(s)
        lines = cleaned

    # detect separator via midpoint line
    mid_line = lines[mid_idx]
    sep = '\t' if mid_line.count('\t') > cols_n_min else ','

    # detect typical separator count in a window (for blank-line removal)
    win = lines[max(0, mid_idx - 2): min(n_lines, mid_idx + 5)]
    sep_counts = [l.count(sep) for l in win if l.strip()]
    typic_seps = min(sep_counts) if sep_counts else 0

    # remove blank/sparse lines
    if typic_seps > 0:
        lines = [l for l in lines if l.count(sep) >= typic_seps // 2]

    content = '\n'.join(lines)

    # detect header row offset 0..4
    df_raw = None
    for offset in range(rows_head_n_max + 1):
        try:
            candidate = pd.read_csv(StringIO(content), sep=sep,
                                    skiprows=offset, dtype=str,
                                    low_memory=False)
            col_list = list(candidate.columns)
            n_unnamed = sum(1 for c in col_list
                            if str(c).startswith('Unnamed') or str(c).strip() == '')
            if 'ExperimentName' in col_list and n_unnamed < 5:
                df_raw = candidate
                break
        except Exception:
            continue

    if df_raw is None:
        return None, 1, 'read_eprime: no valid header found in {}'.format(fname), set()

    # (checked.txt no longer written)

    # extract requested columns (skip those absent in file)
    available = [c for c in cols if c in df_raw.columns]
    rename = {c: fieldnames[cols.index(c)] for c in available}
    available_fields = set(rename.values())
    df_out = df_raw[available].rename(columns=rename).copy()

    # fill in missing columns as NaN
    for fn in fieldnames:
        if fn not in df_out.columns:
            df_out[fn] = np.nan

    df_out = df_out[fieldnames].reset_index(drop=True)

    if verbose:
        print('read_eprime: {} rows from {}'.format(len(df_out), fname))

    return df_out, 0, '', available_fields


# ── internal: task auto-detection ─────────────────────────────────────────

def detect_task(experiment_name):
    """Detect ABCD task from ExperimentName string.

    Returns (task, behavioral_only, scanner_type).
    task: 'MID' | 'SST' | 'NBACK' | 'REC' | None
    scanner_type: 'GE' | 'behavioral' | 'Siemens'
    """
    name = str(experiment_name).lower()
    behavioral_only = 'behavioral' in name
    if '_ge' in name:
        scanner_type = 'GE'
    elif behavioral_only:
        scanner_type = 'behavioral'
    else:
        scanner_type = 'Siemens'

    # 'rec' checked before 'nback'/'wm': nBack_Rec files go to REC
    if 'rec' in name:
        task = 'REC'
    elif 'mid' in name:
        task = 'MID'
    elif 'sst' in name:
        task = 'SST'
    elif 'nback' in name or 'wm' in name:
        task = 'NBACK'
    else:
        task = None

    return task, behavioral_only, scanner_type


# ── shared utilities ───────────────────────────────────────────────────────

def clean_fstem(s):
    """Strip problematic characters from a filename stem."""
    s = os.path.basename(s)
    s = os.path.splitext(s)[0]
    s = re.sub(r'[^A-Za-z0-9_\-]', '_', s)
    return s


def set_ref(onsets, start_times):
    """Match each onset to the nearest preceding start_time.

    Returns integer array of indices into start_times (0-based run index).
    """
    onsets = np.asarray(onsets, dtype=float)
    start_times = np.asarray(start_times, dtype=float)
    run_idx = np.zeros(len(onsets), dtype=int)
    for i, onset in enumerate(onsets):
        valid = np.where(start_times <= onset)[0]
        if len(valid) > 0:
            run_idx[i] = valid[-1]
    return run_idx


def check_offsets(df, onset_col, offset_col, verbose=True):
    """Drop rows with missing onset or offset; warn if any are dropped."""
    o = pd.to_numeric(df[onset_col], errors='coerce')
    f = pd.to_numeric(df[offset_col], errors='coerce')
    mask = o.notna() & f.notna()
    n_drop = (~mask).sum()
    if n_drop > 0 and verbose:
        print('check_offsets: dropping {} rows with missing onset/offset'.format(n_drop))
    return df[mask].copy().reset_index(drop=True)


def _derive_switch_codes(cresp):
    """Derive the two response codes to flip from cresp's own unique values.

    Device/site button-box wiring varies (e.g. some Philips sites use 3/4
    instead of 1/2), so the codes to swap must be read from the file's own
    correct-response values rather than assumed -- same philosophy as
    _build_resp_name_map's data-derived resp<->name pairing.
    Returns (code_a, code_b) when cresp has exactly 2 unique non-NaN
    values; None when indeterminate (0, 1, or more than 2 values).
    """
    uniq = sorted(pd.to_numeric(cresp, errors='coerce').dropna().unique().tolist())
    if len(uniq) == 2:
        return uniq[0], uniq[1]
    return None


def button_switch(df, resp_col, cresp_col, switch_thresh=0.3,
                  extra_resp_cols=None, min_valid_trials=10, verbose=True):
    """Detect and apply a response-code swap if accuracy is below threshold.

    The two codes to flip are derived from cresp_col's own unique values
    (see _derive_switch_codes), not assumed to be 1/2. Reverts if swapped
    accuracy is not better than original. With fewer than min_valid_trials
    valid (resp & cresp both present) trials, the accuracy estimate is too
    noisy to act on, so no switch is attempted at all; if cresp_col does
    not yield exactly 2 unique values, the code pair is indeterminate and
    no switch is attempted either.
    extra_resp_cols: additional response columns to swap alongside resp_col
    (e.g. fix_resp, stop_resp, ssd_resp for SST).
    Returns (df_out, switched).
    """
    def _swap_col(series, code_a, code_b):
        s = pd.to_numeric(series, errors='coerce')
        sw = s.copy()
        sw[s == code_a] = code_b
        sw[s == code_b] = code_a
        # preserve non-numeric values where original was NaN
        return sw.where(s.notna(), series)

    df = df.copy()
    resp = pd.to_numeric(df[resp_col], errors='coerce')
    cresp = pd.to_numeric(df[cresp_col], errors='coerce')

    valid = resp.notna() & cresp.notna()
    if valid.sum() < min_valid_trials:
        return df, False

    codes = _derive_switch_codes(cresp)
    if codes is None:
        if verbose:
            print('button_switch: {} does not yield a unique response-code '
                  'pair, skipping switch'.format(cresp_col))
        return df, False
    code_a, code_b = codes

    orig_acc = (resp[valid] == cresp[valid]).mean()
    if orig_acc >= switch_thresh:
        return df, False

    resp_sw = _swap_col(df[resp_col], code_a, code_b)
    sw_acc = (pd.to_numeric(resp_sw, errors='coerce')[valid] == cresp[valid]).mean()

    if sw_acc <= orig_acc:
        return df, False

    df[resp_col] = resp_sw
    if extra_resp_cols:
        for col in extra_resp_cols:
            if col in df.columns:
                df[col] = _swap_col(df[col], code_a, code_b)
    if verbose:
        msg = 'button_switch: switched {} (orig_acc={:.3f} sw_acc={:.3f})'.format(
            resp_col, orig_acc, sw_acc)
        if (code_a, code_b) != (1.0, 2.0):
            msg += ' using derived codes {}<->{}'.format(code_a, code_b)
        print(msg)
    return df, True


# ── private helpers ────────────────────────────────────────────────────────

def _num(val):
    """Convert val to float; return NaN on failure."""
    try:
        return float(val)
    except (TypeError, ValueError):
        return np.nan


def _isnan(x):
    try:
        return math.isnan(float(x))
    except (TypeError, ValueError):
        return True


def _na_or(val, fmt='{:.3f}'):
    """Return formatted value or 'n/a' when NaN/None."""
    if _isnan(val):
        return 'n/a'
    return fmt.format(float(val))


def _hhmmss_to_seconds(val):
    """Parse an HHMMSS-style value (e.g. 134938 or '112137.254') to seconds
    since midnight; NaN for unparseable/missing values."""
    try:
        s = str(val).strip()
        if not s or s.lower() in ('nan', 'n/a', 'na', 'none', ''):
            return np.nan
        f = float(s)
    except (TypeError, ValueError):
        return np.nan
    hh = int(f // 10000)
    mm = int((f % 10000) // 100)
    ss = f % 100
    return hh * 3600 + mm * 60 + ss


def _parse_button_response(val):
    """Parse E-Prime response string to a button code float (NaN on failure).

    Standard ABCD uses buttons 1 and 2; some sites (e.g. Philips) use 3 and 4.
    Accepts any leading digit 1-9 to handle all button-box configurations.
    Also handles patterns like '1{LEFTARROW}', 'LEFTARROW'→2, 'RIGHTARROW'→1.
    """
    s = str(val).strip()
    if not s or s.lower() in ('nan', 'none', ''):
        return np.nan
    # plain leading digit (1-9): accept any button code
    m = re.match(r'^([1-9])', s)
    if m:
        return float(m.group(1))
    # arrow names
    sl = s.upper()
    if 'LEFTARROW' in sl:
        return 2.0
    if 'RIGHTARROW' in sl:
        return 1.0
    return np.nan


def _parse_cresp(val):
    """Parse E-Prime CorrectResponse; strip trailing text like ',{LEFTARROW}'."""
    s = str(val).strip()
    if not s or s.lower() in ('nan', 'none', ''):
        return np.nan
    # remove trailing ',{...}' or '{...}'
    s2 = re.sub(r'[,\s]*\{[^}]*\}.*', '', s).strip()
    try:
        return float(s2)
    except ValueError:
        pass
    return _parse_button_response(val)


def _sst_resp_is_plain_numeric(val):
    s = str(val).strip()
    return bool(re.match(r'^-?\d+(\.\d+)?$', s))


_BARE_ARROW_MAP = {'LEFTARROW': 1.0, 'RIGHTARROW': 2.0}


def _build_resp_name_map(raw_values, comma_required, letter_fallback,
                         bare_arrow_map=_BARE_ARROW_MAP):
    """Derive {NAME: numeric code} from unique raw response-field string
    values, shared by SST, NBACK, and REC (each with its own parameters).

    Site/device button-box wiring varies, so the resp<->number
    correspondence must be read from the file's own declared pairing
    rather than assumed from a fixed table.

    Handles, in order, per unique non-empty non-numeric raw value (unique
    values sorted, matching MATLAB's unique()):
      - every unique value is a key of bare_arrow_map -> that direct
        mapping (pass bare_arrow_map=None to disable this fallback)
      - 'N,{NAME}' (comma required) or 'N{NAME}'/'N,{NAME}' (comma
        optional, per comma_required) -> {NAME: N}
      - single lowercase r/g/b/y character, if letter_fallback -> numbered
        by 1-based sorted position

    Returns (name_map, errcode, errmsg). name_map is {} (errcode=0) when no
    raw value needed string parsing (already plain numeric/empty).
    name_map is None (errcode=1) when a value matches no enabled pattern.
    """
    seen = set()
    for v in raw_values:
        s = str(v).strip()
        if not s or s.lower() in ('nan', 'none') or _sst_resp_is_plain_numeric(s):
            continue
        # unwrap a value entirely wrapped in '{...}' (matches _apply_resp_map's
        # unwrap step) so e.g. '{LEFTARROW}' is treated the same as 'LEFTARROW'
        m = re.match(r'^\{(\w+)\}$', s)
        if m:
            s = m.group(1)
        seen.add(s)
    uniq = sorted(seen)

    if not uniq:
        return {}, 0, ''

    if bare_arrow_map and all(u in bare_arrow_map for u in uniq):
        return {u: bare_arrow_map[u] for u in uniq}, 0, ''

    comma = ',' if comma_required else ',?'
    pattern = r'(\d+)' + comma + r'\{(\w+)\}'
    mapping = {}
    for i, u in enumerate(uniq, start=1):
        m = re.search(pattern, u)
        if m:
            mapping[m.group(2)] = float(m.group(1))
            continue
        if letter_fallback:
            m2 = re.match(r'^[rgby]$', u)
            if m2:
                mapping[u] = float(i)
                continue
        return None, 1, 'string response value with unexpected pattern ({!r})'.format(u)
    return mapping, 0, ''


def _apply_resp_map(val, name_map):
    """Map a raw response value to numeric via name_map (exact match,
    shared by SST, NBACK, and REC), or coerce directly if already numeric.

    Before the name-map lookup, tries extracting an embedded digit
    directly (handles a self-referential value already in
    'N,{NAME}'/'N{NAME}' format, e.g. SST's go_cresp when Go.RESP is
    plain numeric but Go.CRESP is itself 'N,{NAME}' -- matches MATLAB's
    resolve_resp_value(), sst.m). Otherwise, a value entirely wrapped in
    '{...}' is unwrapped to its bare name -- mirrors mmil_parsecells's
    iscellarr(), which unwraps a value ENTIRELY {wrapped} from position 1
    into the bare name, but leaves a digit-prefixed value like
    '2,{RIGHTARROW}' untouched (handled by the embedded-digit step above).
    """
    s = str(val).strip()
    if not s or s.lower() in ('nan', 'none'):
        return np.nan
    if _sst_resp_is_plain_numeric(s):
        return float(s)
    m0 = re.search(r'(\d+),?\{(\w+)\}', s)
    if m0:
        return float(m0.group(1))
    m = re.match(r'^\{(\w+)\}$', s)
    if m:
        s = m.group(1)
    return name_map.get(s, np.nan)


def _get_experiment_info(df):
    """Extract experiment_name, experiment_version, and runtime_version from DataFrame."""
    exp_name = ''
    exp_ver = ''
    runtime_ver = ''
    if 'experiment_name' in df.columns:
        vals = df['experiment_name'].dropna()
        if len(vals) > 0:
            exp_name = str(vals.iloc[0])
    if 'experiment_version' in df.columns:
        vals = df['experiment_version'].dropna()
        if len(vals) > 0:
            exp_ver = str(vals.iloc[0])
    if 'runtime_version' in df.columns:
        vals = df['runtime_version'].dropna()
        if len(vals) > 0:
            runtime_ver = str(vals.iloc[0])
    return exp_name, exp_ver, runtime_ver


def _find_run_times_mid_sst(df, scanner_type, behavioral_only,
                             trig_col, getready_col,
                             start_delay_Siemens=6.4,
                             expected_start_delay_GE=12.0,
                             TR=0.8, skipTRs=0):
    """Compute (trig_times, start_times) per run for MID and SST.

    Groups trigger-wait rows by ROW-INDEX CONTIGUITY, matching MATLAB's
    get_event_info() (abcd_extract_eprime_mid.m/_sst.m): ind_trig = rows
    where getready_col is non-empty; a new run starts only where there is
    a gap in ROW INDEX (diff(ind_trig)>1), i.e. some intervening row(s)
    lack a value (real task trials in between). The magnitude of the time
    gap between two trigger rows is never inspected -- a value-gap
    threshold (the previous approach here) can spuriously split one true
    run's own contiguous trigger-wait rows into two if the raw file has an
    atypical timing gap within that sequence.

    Returns lists in milliseconds.
    """
    # collect trigger candidate rows, keeping their DataFrame row index
    grrt = pd.to_numeric(df.get(getready_col, pd.Series(dtype=float)),
                         errors='coerce')
    tt = pd.to_numeric(df.get(trig_col, pd.Series(dtype=float)),
                       errors='coerce')

    # prefer getready_rttime; fall back to trig_time
    candidates = grrt if grrt.notna().sum() >= tt.notna().sum() else tt
    valid = candidates[candidates > 0].dropna()

    if len(valid) == 0:
        return [], []

    idx = valid.index.to_numpy()
    group_breaks = np.where(np.diff(idx) > 1)[0] + 1
    groups = np.split(valid.to_numpy(), group_breaks)

    trig_times = []
    start_times = []

    for grp in groups:
        trig_ms = float(grp[0])
        if scanner_type == 'GE':
            start_ms = float(grp[-1])
            # check for missed GE triggers: matches MATLAB's
            # start_time(idx_miss) = start_time(idx_miss) - start_delay_diff(idx_miss)
            # -- subtracts the excess delay (not adds it) to correct for
            # 1-3 missed scanner pulses.
            delta = start_ms - trig_ms
            expected_delta = expected_start_delay_GE * 1000.0
            diff = delta - expected_delta
            n_missed = int(round(diff / (TR * 1000.0)))
            if 1 <= n_missed <= 3 and abs(diff - n_missed * TR * 1000.0) < 0.01 * TR * 1000.0:
                start_ms -= diff
        elif behavioral_only:
            start_ms = trig_ms
        else:
            # single-trigger files: use skipTRs*TR instead of the fixed
            # start_delay_Siemens constant (numerically the same for
            # genuine Siemens/Philips hardware)
            start_ms = trig_ms + skipTRs * TR * 1000.0

        trig_times.append(trig_ms)
        start_times.append(start_ms)

    return trig_times, start_times


def _get_init_times(df, init_col, start_times):
    """Return init_time per run, one value per run in file row order --
    matches MATLAB's algorithm (abcd_extract_eprime_mid.m/_sst.m): collect
    every non-empty, non-zero value of init_col in row order and assign
    them positionally to runs (1st value -> run 1, 2nd -> run 2, ...).

    This is NOT a value-threshold search against start_time: a
    "first value > start_time" search can never select a value earlier
    than start_time, which silently masks a real negative init_delay
    (falls back to init_time=start_time, i.e. init_delay=0) when the true
    init event genuinely precedes the nominal start_time.
    """
    col = pd.to_numeric(df.get(init_col, pd.Series(dtype=float)),
                        errors='coerce')
    nonzero = col[col.notna() & (col != 0)]
    init_times = []
    for i, st in enumerate(start_times):
        init_times.append(float(nonzero.iloc[i]) if i < len(nonzero) else st)
    return init_times


def _next_row_value(series, idx, fallback):
    """Return series at position idx+1 (the row immediately after idx,
    matching MATLAB's ind_init = ind_start + 1 positional lookup) if that
    row exists and is non-NaN, else fallback."""
    next_idx = idx + 1
    if next_idx in series.index:
        v = series.loc[next_idx]
        if not _isnan(v):
            return float(v)
    return fallback


# ── internal: timing validation ────────────────────────────────────────────

def check_eprime_timing(trig_times, start_times, init_times, scanner_type,
                        behavioral_only=False,
                        expected_start_delay_GE=12.0,
                        expected_start_delay_Siemens=6.4,
                        expected_init_delays=None,
                        max_start_diff=0.5, max_init_diff=5.0,
                        verbose=True):
    """Validate E-Prime timing per run; return list of per-run info dicts.

    expected_init_delays may be a float (same for all runs) or a list.
    Does not write any files; caller writes the result into the events columns (start_delay etc.).
    """
    n_runs = len(trig_times)

    if behavioral_only:
        exp_start = 0.0
    elif scanner_type == 'GE':
        exp_start = expected_start_delay_GE
    else:
        exp_start = expected_start_delay_Siemens

    if expected_init_delays is None:
        expected_init_delays = [0.0] * n_runs
    elif not hasattr(expected_init_delays, '__len__'):
        expected_init_delays = [expected_init_delays] * n_runs

    runs_info = []
    for r in range(n_runs):
        tt_s = trig_times[r] / 1000.0
        st_s = start_times[r] / 1000.0
        it_s = init_times[r] / 1000.0
        exp_init = float(expected_init_delays[r])

        start_delay = st_s - tt_s
        init_delay = it_s - st_s
        start_diff = start_delay - exp_start
        init_diff = init_delay - exp_init

        warns = []
        if not behavioral_only:
            if abs(start_diff) > max_start_diff:
                warns.append(
                    'run {}: start_delay={:.3f}s expected={:.3f}s diff={:.3f}s'.format(
                        r + 1, start_delay, exp_start, start_diff))
            if abs(init_diff) > max_init_diff:
                warns.append(
                    'run {}: init_delay={:.3f}s expected={:.3f}s diff={:.3f}s'.format(
                        r + 1, init_delay, exp_init, init_diff))

        if warns and verbose:
            for w in warns:
                print('check_eprime_timing: WARNING: ' + w)

        runs_info.append({
            'run': r + 1,
            'trig_time_s': round(tt_s, 3),
            'start_time_s': round(st_s, 3),
            'start_delay_s': round(start_delay, 3),
            'expected_start_delay_s': exp_start,
            'start_delay_diff_s': round(start_diff, 3),
            'init_time_s': round(it_s, 3),
            'init_delay_s': round(init_delay, 3),
            'expected_init_delay_s': exp_init,
            'init_delay_diff_s': round(init_diff, 3),
            'timing_ok': len(warns) == 0,
            'timing_warnings': warns,
        })

    return runs_info


# ── task processor: MID ────────────────────────────────────────────────────

def process_mid(df, outstem, scanner_type, behavioral_only,
                TR=0.8, skipTRs=0, run_offset=1, available_fields=None,
                verbose=True):
    """Extract MID events from parsed E-Prime DataFrame.

    available_fields: set of fieldnames whose raw E-Prime column was
    genuinely present in the file's header (from read_eprime); used for
    feedback_window/fbflag detection. None/absent assumes new format.

    Returns (events_df, timing_info, meta, errcode, errmsg).
    events_df columns: onset, duration, response_time, trial_type,
                       run_num, response_correct, money
    """
    _TYPENAMES = ['SmallReward', 'LgReward', 'SmallPun', 'LgPun', 'Triangle']
    _CONDNAMES = ['small_reward', 'large_reward', 'small_loss', 'large_loss', 'neutral']
    expected_numtrials = 50
    start_delay_Siemens = 6.4
    expected_start_delay_GE = 12.0

    df = df.copy()
    for col in ['cue_onset', 'cue_dur', 'cue_rt', 'cue_acc', 'cue_resp',
                'probe_onset', 'probe_dur', 'probe_rt', 'probe_acc',
                'fix2_onset', 'fix2_dur', 'money',
                'trig_time', 'getready_rttime', 'initfix_onset', 'raw_block']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    exp_name, exp_ver, runtime_ver = _get_experiment_info(df)

    # feedback_window: True for newer MID versions with extra response
    # windows in the TextDisplay/Feedback periods. fbflag is true if the
    # raw file's header has a feedback_resp column at all, independent of
    # whether any row has real data in it.
    if available_fields is not None:
        fbflag = 'feedback_resp' in available_fields
    else:
        fbflag = True

    trig_times, start_times = _find_run_times_mid_sst(
        df, scanner_type, behavioral_only,
        trig_col='trig_time', getready_col='getready_rttime',
        start_delay_Siemens=start_delay_Siemens,
        expected_start_delay_GE=expected_start_delay_GE,
        TR=TR, skipTRs=skipTRs)

    if not start_times:
        return None, [], {}, 1, 'process_mid: no trigger times found in {}'.format(outstem)

    n_runs = len(start_times)
    init_times = _get_init_times(df, 'initfix_onset', start_times)

    timing_info = check_eprime_timing(
        trig_times, start_times, init_times, scanner_type,
        behavioral_only=behavioral_only,
        expected_start_delay_GE=expected_start_delay_GE,
        expected_start_delay_Siemens=start_delay_Siemens,
        verbose=verbose)

    # trial rows: any row with a real Condition value, not just rows with a
    # valid cue_onset -- a trial can have a genuine Condition but a missing
    # cue_onset when the raw file was archived before the E-Prime export
    # finished writing. Such a trial still counts toward totals/accuracy;
    # check_offsets drops it from timing files only.
    trial_mask = df['block'].isin(_TYPENAMES)
    df_trials = df[trial_mask].copy().reset_index(drop=True)

    if len(df_trials) == 0:
        return None, timing_info, {}, 1, 'process_mid: no trial rows found'

    # run assignment: use the raw Block column directly (E-Prime stamps
    # every row of a run with the same Block number) rather than matching
    # each trial's timestamp against a nominal start_time -- a run's actual
    # trigger-to-first-trial gap can be shorter than the nominal
    # start_delay_Siemens, which would misassign its first trial.
    # (set_ref/start_times above are only used for relative onset time.)
    run_idx = df_trials['raw_block'].values - 1
    # a missing Block value (truncated file, see trial_mask above) gets
    # assigned to the last run.
    missing_block = df_trials['raw_block'].isna().values
    if missing_block.any():
        run_idx = np.where(missing_block, n_runs - 1, run_idx)
    df_trials['_run_idx'] = run_idx.astype(int)

    events = []
    runs_ok = []
    runs_invalid = []

    for r in range(n_runs):
        run_df = df_trials[df_trials['_run_idx'] == r]
        start_ms = start_times[r]

        n_valid = run_df['money'].notna().sum()
        if n_valid == 0:
            if verbose:
                print('process_mid: run {} has 0 valid trials — skipping'.format(r + run_offset))
            continue

        # the E-Prime run, not a count of the runs kept, so a skipped run
        #   leaves the numbering of the rest, by which run_start_times,
        #   timing_info and run_numbers are indexed, in line
        run_num = r + 1

        if n_valid != expected_numtrials:
            if verbose:
                print('process_mid: run {}: {} valid trials (expected {}), marking invalid (truncated)'.format(
                    run_num, n_valid, expected_numtrials))
            runs_invalid.append(run_num)
        else:
            runs_ok.append(run_num)

        for _, row in run_df.iterrows():
            block_str = str(row.get('block', ''))
            condname = None
            for tn, cn in zip(_TYPENAMES, _CONDNAMES):
                if tn.lower() in block_str.lower():
                    condname = cn
                    break
            if condname is None:
                continue

            cue_ms = row['cue_onset']
            cue_dur_ms = row.get('cue_dur', np.nan)
            probe_ms = _num(row.get('probe_onset', np.nan))
            probe_dur_ms = _num(row.get('probe_dur', np.nan))
            probe_offset_ms = _num(row.get('probe_offset', np.nan))
            cue_acc = row.get('cue_acc', np.nan)
            money_val = row.get('money', np.nan)

            # cue_ms missing (truncated file, see trial_mask comment above):
            # still counted as a trial, just can't be timed.
            rel_onset = np.nan if _isnan(cue_ms) else (cue_ms - start_ms) / 1000.0

            # Duration: cue onset to probe offset. Fall back to
            # probe_onset + probe_dur if offset unavailable.
            if not _isnan(probe_offset_ms) and probe_offset_ms != 0:
                antic_dur_s = (probe_offset_ms - cue_ms) / 1000.0
            elif not _isnan(probe_ms) and not _isnan(probe_dur_ms):
                antic_dur_s = (probe_ms + probe_dur_ms - cue_ms) / 1000.0
            else:
                antic_dur_s = cue_dur_ms / 1000.0 if not _isnan(cue_dur_ms) else np.nan

            # Response time fallback order: premature -> n/a; correct -> prbrt;
            # incorrect+pressed -> tdrt-probe_onset; else fbrt-probe_onset;
            # no response -> n/a. RT==0 is never genuine (E-Prime's
            # not-recorded sentinel, matches get_rt(...,1)); tdrt/fbrt
            # already guard on !=0.
            response_check = str(row.get('response_check', ''))
            tdrt = _num(row.get('text_display_rttime', np.nan))
            fbrt = _num(row.get('feedback_rttime', np.nan))
            prbrt_val = _num(row.get('cue_rt', np.nan))
            if 'too soon' in response_check.lower():
                antic_rt_s = np.nan
            elif not _isnan(cue_acc) and int(_num(cue_acc)) == 1:
                antic_rt_s = (prbrt_val / 1000.0
                             if not _isnan(prbrt_val) and prbrt_val != 0
                             else np.nan)
            elif not _isnan(tdrt) and tdrt != 0 and not _isnan(probe_ms):
                antic_rt_s = (tdrt - probe_ms) / 1000.0
            elif not _isnan(fbrt) and fbrt != 0 and not _isnan(probe_ms):
                antic_rt_s = (fbrt - probe_ms) / 1000.0
            else:
                antic_rt_s = np.nan

            early = 1 if 'too soon' in response_check.lower() else 0
            # response_antic: any response (premature or valid), matching
            # check_responses() -- independent of RT, so trialresp
            # counting doesn't depend on response_time.
            responded = (not pd.isna(row.get('probe_resp', np.nan))
                        or not pd.isna(row.get('text_display_resp', np.nan))
                        or not pd.isna(row.get('feedback_resp', np.nan)))
            events.append({
                'onset': 'n/a' if _isnan(rel_onset) else round(rel_onset, 3),
                'duration': _na_or(antic_dur_s),
                'response_time': _na_or(antic_rt_s),
                'trial_type': condname + '_antic',
                'run_num': run_num,
                'response_correct': 'n/a',
                'money': 'n/a',
                'response_early': early,
                'response_antic': int(responded),
                'feedback_window': int(fbflag),
            })

            # pos_feedback/neg_feedback classification is driven by accuracy
            # alone, independent of whether feedback_onset survived -- so the
            # event is always emitted; only its onset/duration fall back to
            # 'n/a' when feedback_onset itself is missing (same
            # truncated-file case as cue_onset above).
            fb_onset_ms = _num(row.get('feedback_onset', np.nan))
            fb_offset_ms = _num(row.get('feedback_offset', np.nan))

            # premature response recodes as incorrect (response_check already extracted above)
            if 'too soon' in response_check.lower():
                fb_acc = 0
            else:
                fb_acc = int(_num(cue_acc)) if not _isnan(cue_acc) else 0

            fb_type = 'pos_feedback' if fb_acc == 1 else 'neg_feedback'

            if not _isnan(fb_onset_ms) and fb_onset_ms != 0:
                fb_onset_val = round((fb_onset_ms - start_ms) / 1000.0, 3)
                fb_dur_s = ((fb_offset_ms - fb_onset_ms) / 1000.0
                           if not _isnan(fb_offset_ms) and fb_offset_ms != 0 else np.nan)
            else:
                fb_onset_val = 'n/a'
                fb_dur_s = np.nan

            events.append({
                'onset': fb_onset_val,
                'duration': _na_or(fb_dur_s),
                'response_time': _na_or(antic_rt_s),
                'trial_type': condname + '_' + fb_type,
                'run_num': run_num,
                'response_correct': fb_acc,
                'money': round(float(money_val), 2) if not _isnan(money_val) else 'n/a',
                'response_early': 'n/a',
                'response_antic': 'n/a',
                'feedback_window': int(fbflag),
            })

    if not events:
        return None, timing_info, {}, 1, 'process_mid: no events produced'

    events_df = pd.DataFrame(events)
    meta = {
        'task': 'MID',
        'experiment_name': exp_name,
        'experiment_version': exp_ver,
        'runtime_version': runtime_ver,
        'behavioral_only': behavioral_only,
        'scanner_type': scanner_type,
        'switch_flag': False,
        'runs_ok': runs_ok,
        'runs_invalid': runs_invalid,
    }
    return events_df, timing_info, meta, 0, ''


def _normalize_block_type(bt):
    """Map a raw WM block-type/procedure string to '2_back'/'0_back'."""
    if pd.isna(bt):
        return ''
    bt = str(bt).strip()
    if not bt:
        return ''
    if '2' in bt:
        return '2_back'
    elif '0' in bt:
        return '0_back'
    return bt.lower()


def _build_wm_stim_lookup(df_wm):
    """Build {stim basename: '0_back'/'2_back'} from a raw WM DataFrame.

    Uses ALL rows unconditionally, independent of any run-validity/
    truncation determination in the NBACK processing itself.

    Blank filler/instruction rows (block_type, stim_type, or target_type
    empty) are skipped -- a block_type present with stim_type or
    target_type empty can reflect a truncated file (premature file
    transfer).
    """
    lookup = {}
    if df_wm is None or 'stim' not in df_wm.columns:
        return lookup
    for _, row in df_wm.iterrows():
        if (pd.isna(row.get('block_type')) or pd.isna(row.get('stim_type')) or
                pd.isna(row.get('target_type'))):
            continue
        stim_raw = row.get('stim', '')
        stim_s = '' if pd.isna(stim_raw) else str(stim_raw).strip()
        if not stim_s:
            continue
        basename = os.path.basename(stim_s)
        if basename in lookup:
            continue
        bt_label = _normalize_block_type(row.get('block_type', ''))
        if bt_label:
            lookup[basename] = bt_label
    return lookup


# ── task processor: NBACK ─────────────────────────────────────────────────

def process_nback(df, outstem, scanner_type, behavioral_only,
                  TR=0.8, skipTRs=0, switch_thresh=0.3, run_offset=1,
                  error_flag=False, extra_flag=False, verbose=True):
    """Extract NBACK events from parsed E-Prime DataFrame.

    Returns (events_df, timing_info, meta, errcode, errmsg).
    events_df extra columns: response_accuracy, block_type, stim_type, target_type
    """
    expected_numtrials = 80
    start_delay_Siemens = 6.4
    expected_start_delay_GE = 12.0

    df = df.copy()

    # parse string responses to numeric: derive the name<->code map from
    # correct_response's own declared values (site button-box wiring
    # varies), matching MATLAB nback.m's own regex + fallback logic
    name_map, ec_map, em_map = _build_resp_name_map(
        df['correct_response'].tolist(), comma_required=True, letter_fallback=True)
    if ec_map:
        return None, [], {}, ec_map, 'process_nback: ' + em_map

    raw_stim_resp = df['stim_resp']
    has_string_stim_resp = raw_stim_resp.apply(
        lambda v: str(v).strip() != '' and str(v).strip().lower() not in ('nan', 'none')
        and not _sst_resp_is_plain_numeric(str(v).strip())).any()
    if has_string_stim_resp and not name_map:
        errmsg = 'string stim_resp values without string correct_response values'
        if verbose:
            print('process_nback: {}'.format(errmsg))
        return None, [], {}, 1, errmsg

    df['stim_resp'] = raw_stim_resp.apply(lambda v: _apply_resp_map(v, name_map))
    df['correct_response'] = df['correct_response'].apply(_parse_cresp)

    for col in ['stim_onset', 'stim_offset', 'stim_dur', 'stim_rt', 'stim_acc',
                'cuefix_onset', 'cuefix_dur', 'cue2back_offset', 'cuetarget_offset',
                'trig_time', 'getready_rttime', 'getready2_rttime']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    exp_name, exp_ver, runtime_ver = _get_experiment_info(df)
    block_col = df['block'].fillna('').astype(str)

    # identify TRSync trigger rows (run 1 and run 2)
    is_trsync1 = (block_col.str.contains('TRSyncPROC', case=False) &
                  ~block_col.str.contains('R2', case=False))
    is_trsync2 = (block_col.str.contains('TRSyncPROC', case=False) &
                  block_col.str.contains('R2', case=False))

    trig_times = []
    start_times = []
    init_times = []

    def _run_timing_nback(trsync_mask, grrt_col, run_label):
        trsync_df = df[trsync_mask]
        if len(trsync_df) == 0:
            return None, None, None
        grrt = pd.to_numeric(trsync_df[grrt_col], errors='coerce').dropna()
        grrt = grrt[grrt > 0]
        if len(grrt) == 0:
            return None, None, None
        trig_ms = float(grrt.min())
        if scanner_type == 'GE':
            start_ms = float(grrt.max())
            # check for missed GE triggers (same logic as _find_run_times_mid_sst):
            # matches MATLAB's start_time(idx_miss) = start_time(idx_miss) -
            # start_delay_diff(idx_miss) -- subtracts the excess delay (not
            # adds it) to correct for 1-3 missed scanner pulses.
            delta = start_ms - trig_ms
            expected_delta = expected_start_delay_GE * 1000.0
            diff = delta - expected_delta
            n_missed = int(round(diff / (TR * 1000.0)))
            if 1 <= n_missed <= 3 and abs(diff - n_missed * TR * 1000.0) < 0.01 * TR * 1000.0:
                start_ms -= diff
                if verbose:
                    print('process_nback: run {}: adjusted start_time for '
                          '{} missed GE trigger(s)'.format(run_label, n_missed))
        elif behavioral_only:
            start_ms = trig_ms
        else:
            # single-trigger files: use skipTRs*TR instead of the fixed
            # start_delay_Siemens constant (numerically the same for
            # genuine Siemens/Philips hardware)
            start_ms = trig_ms + skipTRs * TR * 1000.0
        return trig_ms, start_ms, trsync_df.index[-1]

    trig1, start1, last_trsync1_idx = _run_timing_nback(is_trsync1, 'getready_rttime', 1)
    trig2, start2, last_trsync2_idx = _run_timing_nback(is_trsync2, 'getready2_rttime', 2)

    # fallback: cluster approach if TRSync rows not found
    if trig1 is None:
        trig_list, start_list = _find_run_times_mid_sst(
            df, scanner_type, behavioral_only,
            trig_col='trig_time', getready_col='getready_rttime',
            start_delay_Siemens=start_delay_Siemens,
            expected_start_delay_GE=expected_start_delay_GE, TR=TR,
            skipTRs=skipTRs)
        if not trig_list:
            return None, [], {}, 1, 'process_nback: no trigger times found'
        trig_times = trig_list
        start_times = start_list
        init_times = _get_init_times(df, 'cuefix_onset', start_times)
    else:
        trig_times = [trig1]
        start_times = [start1]
        # init_time for run 1: cuefix_onset at the row immediately after
        # the last TRSyncPROC-matching row -- a positional lookup
        # (matching MATLAB's ind_init = ind_start + 1), not a
        # value-threshold search. cuefix_onset can carry a stale/repeated
        # value forward into the TRSync rows themselves, which a
        # "first value > start_time anywhere in df" search would
        # wrongly match instead of the real next-cue value.
        cuefix = pd.to_numeric(df['cuefix_onset'], errors='coerce')
        init_times.append(_next_row_value(cuefix, last_trsync1_idx, start1))

        if trig2 is not None:
            trig_times.append(trig2)
            start_times.append(start2)
            init_times.append(_next_row_value(cuefix, last_trsync2_idx, start2))

    n_runs = len(start_times)

    # expected init delays for GE: run1=0.5s, run2=0.0s
    if scanner_type == 'GE':
        exp_init = [0.5 if r == 0 else 0.0 for r in range(n_runs)]
    else:
        exp_init = [0.0] * n_runs

    timing_info = check_eprime_timing(
        trig_times, start_times, init_times, scanner_type,
        behavioral_only=behavioral_only,
        expected_start_delay_GE=expected_start_delay_GE,
        expected_start_delay_Siemens=start_delay_Siemens,
        expected_init_delays=exp_init,
        verbose=verbose)

    # button switch on stimulus trials
    stim_mask = (df['stim_onset'].notna() & (df['stim_onset'] > 0) &
                 df['block_type'].notna() & (df['block_type'].astype(str).str.strip() != '') &
                 df['stim_type'].notna() & (df['stim_type'].astype(str).str.strip() != '') &
                 df['target_type'].notna() & (df['target_type'].astype(str).str.strip() != ''))

    df_stim = df[stim_mask].copy().reset_index(drop=True)
    switched = False
    if len(df_stim) > 0:
        df_stim, switched = button_switch(df_stim, 'stim_resp', 'correct_response',
                                          switch_thresh=switch_thresh, verbose=verbose)
        if switched:
            # button_switch only swaps stim_resp; stim_acc must be
            # recomputed against the swapped response for every trial.
            # correct_response already encodes whichever comparison rule
            # applies (0-back target vs. 2-back-ago stimulus), so one
            # formula covers both block types. Previously stim_acc kept the
            # raw, pre-swap E-Prime Stim.ACC value for the whole session.
            resp = df_stim['stim_resp']
            cresp = df_stim['correct_response']
            valid = resp.notna() & cresp.notna()
            df_stim.loc[valid, 'stim_acc'] = (resp[valid] == cresp[valid]).astype(int)
            try:
                df_stim.to_csv(outstem + '_events_switched.csv', index=False)
            except Exception:
                pass

    run_idx = set_ref(df_stim['stim_onset'].values, np.array(start_times))
    df_stim['_run_idx'] = run_idx
    with np.errstate(invalid='ignore'):
        cue_run_idx = set_ref(df['cuefix_onset'].values, np.array(start_times))

    events = []
    runs_ok = []
    runs_invalid = []

    for r in range(n_runs):
        run_df = df_stim[df_stim['_run_idx'] == r]
        start_ms = start_times[r]

        if len(run_df) == 0:
            if verbose:
                print('process_nback: run {} has no stimulus trials'.format(r + run_offset))
            continue

        # the E-Prime run, not a count of the runs kept, so a skipped run
        #   leaves the numbering of the rest, by which run_start_times,
        #   timing_info and run_numbers are indexed, in line
        run_num = r + 1

        if len(run_df) != expected_numtrials:
            if verbose:
                print('process_nback: run {}: {} trials (expected {}), marking invalid (truncated)'.format(
                    run_num, len(run_df), expected_numtrials))
            runs_invalid.append(run_num)
        else:
            runs_ok.append(run_num)

        # cue events: find cuefix_onset rows for this run that are TRSync-adjacent
        # Use cuefix_onset values not belonging to a stim row, restricted to the
        # two real cue procedures -- Procedure[Block] ('block') must be
        # Cue0BackPROC/Cue2BackPROC, matching MATLAB's parms.cues exactly.
        # Without this, the scanner-trigger-sync procedures (TRSyncPROC/
        # TRSyncPROCR2, 16 rows/run) carry a stuck CueFix.OnsetTime/.Duration
        # value from E-Prime's nested logging that otherwise satisfies this
        # mask too, producing 16 spurious duplicate 'cue' rows per run.
        # Run assignment uses set_ref() (nearest preceding start_time,
        # defaulting to run 0 if none precedes) rather than an explicit
        # start_ms/start_times[r+1] window, matching MATLAB's set_ref()
        # exactly -- a cue landing at/before start_ms still belongs to its run.
        cue_mask = (df['cuefix_onset'].notna() & (cue_run_idx == r) &
                    df['stim_onset'].isna() &
                    df['block'].isin(['Cue0BackPROC', 'Cue2BackPROC']))
        cue_rows = df[cue_mask].copy()

        for _, crow in cue_rows.iterrows():
            cfms = crow['cuefix_onset']
            proc = crow.get('block', '')
            # cue duration: CueFix.OnsetTime to the block-type instructional
            # screen's OWN offset (Cue2Back.OffsetTime for Cue2BackPROC,
            # CueTarget.OffsetTime for Cue0BackPROC/CueTargetPROC) -- not
            # CueFix.Duration alone, which covers only the brief fixation
            # preceding that screen. Matches MATLAB's write_files() caller
            # (tmp_offset = event_info_proc.([cuenames{i} '_offset'])).
            if proc == 'Cue2BackPROC':
                cf_offset = crow.get('cue2back_offset', np.nan)
            else:
                cf_offset = crow.get('cuetarget_offset', np.nan)
            cfdur_s = ((cf_offset - cfms) / 1000.0
                       if not _isnan(cf_offset) and not _isnan(cfms) else np.nan)
            events.append({
                'onset': round((cfms - start_ms) / 1000.0, 3),
                'duration': _na_or(cfdur_s),
                'response_time': 'n/a',
                'trial_type': 'cue',
                'run_num': run_num,
                'response_accuracy': 'n/a',
                'block_type': 'n/a',
                'stim_type': 'n/a',
                'target_type': 'n/a',
                'stim': 'n/a',
            })

        for _, row in run_df.iterrows():
            st_ms = row['stim_onset']
            st_offset = row.get('stim_offset', np.nan)
            # duration = Stim.OffsetTime-Stim.OnsetTime, matching MATLAB's
            # write_files() (which never reads Stim.Duration at all) -- the
            # actual measured on-screen elapsed time, not the nominal
            # Stim.Duration design value (confirmed a systematic ~100ms
            # shorter than nominal in real data, not jitter)
            st_dur = ((st_offset - st_ms)
                      if not _isnan(st_offset) and not _isnan(st_ms) else np.nan)
            st_rt = row.get('stim_rt', np.nan)
            st_resp = row.get('stim_resp', np.nan)
            st_acc = row.get('stim_acc', np.nan)
            bt = str(row.get('block_type', '')).strip()
            stype = str(row.get('stim_type', '')).strip()
            ttype = str(row.get('target_type', '')).strip()

            bt_label = _normalize_block_type(bt)

            stype_l = stype.lower()
            ttype_l = ttype.lower()

            # build trial_type: {block}_{stim}  (e.g. 2_back_posface)
            base_type = bt_label + '_' + stype_l if stype_l else bt_label

            # error variant
            if error_flag and not _isnan(st_acc) and int(st_acc) == 0:
                trial_type = base_type + '_err'
            else:
                trial_type = base_type

            _stim_raw = row.get('stim', '')
            _stim_s = '' if pd.isna(_stim_raw) else str(_stim_raw).strip()
            stim_basename = os.path.basename(_stim_s) if _stim_s else 'n/a'
            # 'n/a' when stim_resp is empty, not just when Stim.RT is blank --
            # Stim.RT=0 for a true omission (E-Prime's convention), not blank
            st_rt_s = (st_rt / 1000.0 if not _isnan(st_resp) and not _isnan(st_rt)
                       else np.nan)
            events.append({
                'onset': round((st_ms - start_ms) / 1000.0, 3),
                'duration': _na_or(st_dur / 1000.0 if not _isnan(st_dur) else np.nan),
                'response_time': _na_or(st_rt_s),
                'trial_type': trial_type,
                'run_num': run_num,
                'response_accuracy': int(st_acc) if not _isnan(st_acc) else 'n/a',
                'block_type': bt_label,
                'stim_type': stype_l if stype_l else 'n/a',
                'target_type': ttype_l if ttype_l else 'n/a',
                'stim': stim_basename,
            })

            # extra target/lure/nonlure events
            if extra_flag and ttype_l:
                events.append({
                    'onset': round((st_ms - start_ms) / 1000.0, 3),
                    'duration': _na_or(st_dur / 1000.0 if not _isnan(st_dur) else np.nan),
                    'response_time': _na_or(st_rt_s),
                    'trial_type': bt_label + '_' + ttype_l,
                    'run_num': run_num,
                    'response_accuracy': int(st_acc) if not _isnan(st_acc) else 'n/a',
                    'block_type': bt_label,
                    'stim_type': stype_l if stype_l else 'n/a',
                    'target_type': ttype_l,
                    'stim': stim_basename,
                })

    if not events:
        return None, timing_info, {}, 1, 'process_nback: no events produced'

    events_df = pd.DataFrame(events)
    meta = {
        'task': 'NBACK',
        'experiment_name': exp_name,
        'experiment_version': exp_ver,
        'runtime_version': runtime_ver,
        'behavioral_only': behavioral_only,
        'scanner_type': scanner_type,
        'switch_flag': switched,
        'runs_ok': runs_ok,
        'runs_invalid': runs_invalid,
    }
    return events_df, timing_info, meta, 0, ''


# ── task processor: SST ────────────────────────────────────────────────────

def _revise_sst_types(df, ste_recode_flag=False):
    """Recode SST trial types from raw responses.

    Sets trial_type to one of: CorrectGo, CorrectLateGo, NoRespGo,
    IncorrectGo, IncorrectLateGo, CorrectStop, IncorrectStop.
    STE recoding: stop trials with ssd_rt>0 are either recoded as Go trials
    (ste_recode_flag=True) or labeled 'STE' (ste_recode_flag=False, default).
    A zero-SSD stop trial with an SSD-window response (ssd_resp) is labeled
    'ZeroSSDResp', neither go nor stop.
    Returns df with added 'trial_type_sst' column.
    """
    df = df.copy()
    trial_types = []

    for _, row in df.iterrows():
        is_stop = not _isnan(row.get('stop_onset', np.nan))
        go_resp = row.get('go_resp', np.nan)
        go_cresp = row.get('go_cresp', np.nan)
        go_rt = row.get('go_rt', np.nan)
        go_dur = row.get('go_dur', np.nan)
        ssd_rt = row.get('ssd_rt', np.nan)
        ssd_dur = row.get('ssd_dur', np.nan)
        fix_resp = row.get('fix_resp', np.nan)
        fix_rt = row.get('fix_rt', np.nan)

        go_acc = (not _isnan(go_resp) and not _isnan(go_cresp) and
                  float(go_resp) == float(go_cresp))

        if is_stop:
            # STE detection: responded during the SSD window (ssd_dur > 0
            # AND ssd_rt > 0; zero-SSD trials are not STE).
            if not _isnan(ssd_rt) and float(ssd_rt) > 0 and not _isnan(ssd_dur) and float(ssd_dur) > 0:
                if ste_recode_flag:
                    # legacy: recode as a Go trial using SSD values. Uses
                    # ssd_resp for response; re-derives go_cresp from the
                    # stimulus filename (LeftArrow->1, RightArrow->2) since
                    # Go.CRESP is 0/empty for stop trials. Case-sensitive
                    # exact substring; if neither matches, ste_cresp=NaN so
                    # the trial falls through to IncorrectGo.
                    ssd_resp_val = row.get('ssd_resp', np.nan)
                    stim = str(row.get('stim_file', '') or '')
                    if 'LeftArrow' in stim:
                        ste_cresp = 1.0
                    elif 'RightArrow' in stim:
                        ste_cresp = 2.0
                    else:
                        ste_cresp = np.nan
                    ste_resp = _num(ssd_resp_val)
                    go_acc_ste = (not _isnan(ste_resp) and not _isnan(ste_cresp) and
                                  ste_resp == ste_cresp)
                    trial_types.append('CorrectGo' if go_acc_ste else 'IncorrectGo')
                else:
                    # default: label as 'STE' nuisance condition; not recoded as Go
                    trial_types.append('STE')
            else:
                stop_resp = row.get('stop_resp', np.nan)
                ssd_resp = row.get('ssd_resp', np.nan)
                if (not _isnan(ssd_dur) and float(ssd_dur) == 0 and
                        not _isnan(ssd_resp) and float(ssd_resp) != 0):
                    # zero-SSD stop trial with an SSD-window response:
                    # neither go nor stop
                    trial_types.append('ZeroSSDResp')
                elif not _isnan(ssd_dur) and float(ssd_dur) == 0:
                    # Zero-SSD stop trial: same response-field logic as
                    # non-zero-SSD trials.
                    if (_isnan(stop_resp) or float(stop_resp) == 0) and (_isnan(fix_resp) or float(fix_resp) == 0):
                        trial_types.append('CorrectStop')
                    else:
                        trial_types.append('IncorrectStop')
                elif _isnan(stop_resp) and (_isnan(fix_resp) or float(fix_resp) == 0):
                    trial_types.append('CorrectStop')
                else:
                    trial_types.append('IncorrectStop')
        else:
            # Go trial
            no_resp = (_isnan(go_resp) or float(go_resp) == 0)
            if no_resp:
                if not _isnan(fix_resp) and float(fix_resp) != 0:
                    # late response in fix window
                    fix_acc = (not _isnan(fix_resp) and not _isnan(go_cresp) and
                               float(fix_resp) == float(go_cresp))
                    trial_types.append('CorrectLateGo' if fix_acc else 'IncorrectLateGo')
                else:
                    trial_types.append('NoRespGo')
            else:
                if go_acc:
                    trial_types.append('CorrectGo')
                else:
                    trial_types.append('IncorrectGo')

    df['trial_type_sst'] = trial_types
    return df


def process_sst(df, outstem, scanner_type, behavioral_only,
                TR=0.8, skipTRs=0, switch_thresh=0.3, run_offset=1,
                ste_recode_flag=False, verbose=True):
    """Extract SST events from parsed E-Prime DataFrame.

    Returns (events_df, timing_info, meta, errcode, errmsg).
    events_df extra columns: response_correct_go, correct_stop, duration_ssd,
                       response_time_go_adjusted
    """
    expected_numtrials = 182
    start_delay_Siemens = 6.4
    expected_start_delay_GE = 12.0

    df = df.copy()

    for col in ['go_onset', 'go_offset', 'go_dur', 'go_rt',
                'ssd_dur', 'ssd_rt',
                'fix_onset', 'fix_dur', 'fix_rt',
                'stop_onset', 'stop_dur', 'stop_rt',
                'ssd_onset', 'stop_start',
                'getready_rttime', 'wait1_rttime', 'wait2_rttime',
                'beginfix_onset', 'trig_time', 'block', 'trial_run']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # go_resp/go_cresp/fix_resp/stop_resp/ssd_resp may be non-numeric strings
    # (e.g. '{LEFTARROW}', '1,{LEFTARROW}') on some devices/sites. The
    # resp<->number correspondence is site-specific (some pair LEFTARROW=1,
    # others LEFTARROW=2), so it must be read from go_cresp's own declared
    # pairing and applied uniformly to every response column, rather than
    # assumed from a fixed table.
    name_map, ec_map, em_map = _build_resp_name_map(
        df['go_cresp'].tolist(), comma_required=False, letter_fallback=True)
    if ec_map:
        return None, [], {}, ec_map, 'process_sst: ' + em_map
    for col in ['go_resp', 'go_cresp', 'fix_resp', 'stop_resp', 'ssd_resp']:
        if col in df.columns:
            df[col] = df[col].apply(lambda v: _apply_resp_map(v, name_map))

    exp_name, exp_ver, runtime_ver = _get_experiment_info(df)

    # run timing: GE uses getready_rttime; Siemens uses wait1/wait2
    if scanner_type == 'GE':
        trig_times, start_times = _find_run_times_mid_sst(
            df, scanner_type, behavioral_only,
            trig_col='trig_time', getready_col='getready_rttime',
            start_delay_Siemens=start_delay_Siemens,
            expected_start_delay_GE=expected_start_delay_GE, TR=TR,
            skipTRs=skipTRs)
    else:
        # Siemens: wait1_rttime (run1), wait2_rttime (run2)
        trig_times = []
        start_times = []
        for wcol in ['wait1_rttime', 'wait2_rttime']:
            if wcol not in df.columns:
                continue
            vals = pd.to_numeric(df[wcol], errors='coerce')
            vals = vals[(vals > 0) & vals.notna()].values
            if len(vals) > 0:
                tm = float(np.min(vals))
                # single-trigger files: use skipTRs*TR instead of the fixed
                # start_delay_Siemens constant (numerically the same for
                # genuine Siemens/Philips hardware)
                sm = tm if behavioral_only else tm + skipTRs * TR * 1000.0
                trig_times.append(tm)
                start_times.append(sm)

        if not trig_times:
            trig_times, start_times = _find_run_times_mid_sst(
                df, scanner_type, behavioral_only,
                trig_col='trig_time', getready_col='getready_rttime',
                start_delay_Siemens=start_delay_Siemens,
                expected_start_delay_GE=expected_start_delay_GE, TR=TR,
                skipTRs=skipTRs)

    if not start_times:
        return None, [], {}, 1, 'process_sst: no trigger times found'

    n_runs = len(start_times)
    init_times = _get_init_times(df, 'beginfix_onset', start_times)

    timing_info = check_eprime_timing(
        trig_times, start_times, init_times, scanner_type,
        behavioral_only=behavioral_only,
        expected_start_delay_GE=expected_start_delay_GE,
        expected_start_delay_Siemens=start_delay_Siemens,
        verbose=verbose)

    # trial rows: go trials have go_onset; stop trials have ssd_onset instead
    trial_mask = (df['go_onset'].notna() & (df['go_onset'] > 0)) | \
                 (df['ssd_onset'].notna() & (df['ssd_onset'] > 0))
    df_trials = df[trial_mask].copy().reset_index(drop=True)

    if len(df_trials) == 0:
        return None, timing_info, {}, 1, 'process_sst: no trial rows found'

    # button switch on go_resp / go_cresp (+ fix_resp, stop_resp, ssd_resp)
    df_trials, switched = button_switch(df_trials, 'go_resp', 'go_cresp',
                                        switch_thresh=switch_thresh, verbose=verbose,
                                        extra_resp_cols=['fix_resp', 'stop_resp', 'ssd_resp'])
    if switched:
        try:
            df_trials.to_csv(outstem + '_events_switched.csv', index=False)
        except Exception:
            pass

    df_trials = _revise_sst_types(df_trials, ste_recode_flag=ste_recode_flag)

    # Run assignment: pick ONE raw column as the run indicator for the whole
    # file -- raw Trial when it has <= 2 unique values (some SST versions
    # encode the run/block number directly in Trial, e.g. 1/2), otherwise
    # raw Block for versions where Trial is a per-trial counter instead.
    # Applied uniformly to go and stop trials alike. Neither source is
    # timestamp-based, so there's no nominal-start_time drift exposure
    # (set_ref/start_times above are only used for relative onset time).
    n_trial_unique = df_trials['trial_run'].dropna().nunique()
    if 0 < n_trial_unique <= 2:
        run_col = df_trials['trial_run'].values
    else:
        run_col = df_trials['block'].values
    run_idx = run_col - 1
    # a trial can have a genuine trial type but a missing run-indicator
    # value when the raw file was archived before the export finished
    # writing; assign such trials to the last run.
    missing = np.isnan(run_col)
    if missing.any():
        run_idx = np.where(missing, n_runs - 1, run_idx)
    df_trials['_run_idx'] = run_idx.astype(int)

    # Structural per-run trial count, matching MATLAB's ind_events filter
    # (~isempty(event_info.type) on raw TrialCode, sst.m ~368-369) --
    # computed from the FULL df, not df_trials, since 'BeginFix'/'EndFix'
    # rows carry a non-empty TrialCode but no go_onset/ssd_onset and are
    # excluded from df_trials's onset-based trial_mask. df_trials/run_df
    # itself is unchanged and still drives all actual trial extraction;
    # only the run-validity count below uses this separate value.
    if 0 < n_trial_unique <= 2:
        full_run_col = df['trial_run'].values
    else:
        full_run_col = df['block'].values
    full_run_idx = full_run_col - 1
    full_missing = np.isnan(full_run_col)
    if full_missing.any():
        full_run_idx = np.where(full_missing, n_runs - 1, full_run_idx)
    has_trial_code = df['trial_code'].notna() & (df['trial_code'].astype(str).str.strip() != '')
    structural_counts = pd.Series(full_run_idx[has_trial_code.values]).value_counts()

    # map trial_type_sst → output trial_type
    _SST_TYPE_MAP = {
        'CorrectGo': 'correct_go',
        'CorrectLateGo': 'correctlate_go',
        'NoRespGo': 'noresp_go',
        'IncorrectGo': 'incorrect_go',
        'IncorrectLateGo': 'incorrectlate_go',
        'CorrectStop': 'correct_stop',
        'IncorrectStop': 'incorrect_stop',
        'STE': 'early_stop',
        'ZeroSSDResp': 'zero_ssd_response',
    }

    events = []
    runs_ok = []
    runs_invalid = []

    for r in range(n_runs):
        run_df = df_trials[df_trials['_run_idx'] == r]
        start_ms = start_times[r]

        if len(run_df) == 0:
            if verbose:
                orig = r + run_offset
                print('process_sst: run {} has no trials'.format(orig))
            continue

        # the E-Prime run, not a count of the runs kept, so a skipped run
        #   leaves the numbering of the rest, by which run_start_times,
        #   timing_info and run_numbers are indexed, in line
        run_num = r + 1

        structural_count = int(structural_counts.get(r, 0))
        if structural_count != expected_numtrials:
            if verbose:
                print('process_sst: run {}: {} trials (expected {}), marking invalid (truncated)'.format(
                    run_num, structural_count, expected_numtrials))
            runs_invalid.append(run_num)
        else:
            runs_ok.append(run_num)

        for _, row in run_df.iterrows():
            # a missing value is read as NaN (a float), not an empty string,
            # and the .lower() default is evaluated before the lookup
            sst_type = row.get('trial_type_sst', '')
            if not isinstance(sst_type, str):
                sst_type = ''
            trial_type = _SST_TYPE_MAP.get(sst_type, sst_type.lower())

            is_stop_trial = sst_type in ('CorrectStop', 'IncorrectStop')
            # zero-SSD stop trial with an SSD-window response: neither go nor stop
            is_zr = sst_type == 'ZeroSSDResp'
            is_go_trial = not is_stop_trial and not is_zr

            go_ms = row['go_onset']
            go_offset_ms = row.get('go_offset', np.nan)
            go_dur = row.get('go_dur', np.nan)
            go_rt = row.get('go_rt', np.nan)
            fix_rt = row.get('fix_rt', np.nan)
            fix_resp = row.get('fix_resp', np.nan)
            stop_rt = row.get('stop_rt', np.nan)
            stop_ms = row.get('stop_onset', np.nan)
            stop_start_ms = row.get('stop_start', np.nan)
            ssd_d = row.get('ssd_dur', np.nan)
            ssd_rt_val = row.get('ssd_rt', np.nan)

            # STE-recoded go trial: was a stop trial (stop_onset > 0, ssd_dur > 0)
            # but responded during the SSD window (ssd_rt > 0). Requires
            # ssd_dur > 0 AND ssd_rt > 0 -- zero-SSD trials are not STE even
            # if ssd_rt > 0.
            is_ste = (is_go_trial and
                      not _isnan(stop_ms) and float(stop_ms) > 0 and
                      not _isnan(ssd_d) and float(ssd_d) > 0 and
                      not _isnan(ssd_rt_val) and float(ssd_rt_val) > 0)

            if is_stop_trial:
                # stop event onset: StopSignal.StartTime; offset = onset +
                # 20 ms (effectively a point event)
                onset_ms = stop_start_ms if not _isnan(stop_start_ms) else stop_ms
                dur_ms = 20.0
            elif is_ste or is_zr:
                # STE or ZeroSSDResp: use SSD onset; duration = stop_dur + ssd_dur
                onset_ms = stop_ms  # SSD.OnsetTime (ssd_onset)
                stop_dur_val = row.get('stop_dur', np.nan)
                dur_ms = (stop_dur_val + ssd_d) if (not _isnan(stop_dur_val) and not _isnan(ssd_d)) else np.nan
            else:
                # duration = Go.OffsetTime - Go.OnsetTime, matching MATLAB's
                # write_files() (which reads go_onset_time/go_offset_time,
                # not Go.Duration). This is normally response-terminated
                # (a real, per-trial value), not the nominal fixed Go.Duration
                # -- Go.Duration is a separate field, used only for the
                # late-fix-window response_time_go_adjusted calc below, never for the
                # event's own timing/duration.
                onset_ms = go_ms
                dur_ms = ((go_offset_ms - go_ms)
                         if (not _isnan(go_offset_ms) and not _isnan(go_ms))
                         else np.nan)

            # adjusted go RT. Gated on go_resp being empty, not on
            # go_rt==0: E-Prime's PreRelease can log a response against the
            # wrong (following) object, leaving that trial's own RT at 0
            # despite a real go_resp -- a documented E-Prime possibility,
            # not something confirmed in this project's own data.
            go_resp = row.get('go_resp', np.nan)
            go_resp_empty = _isnan(go_resp) or float(go_resp) == 0
            if is_ste:
                go_rt_adj = ssd_rt_val / 1000.0 if not _isnan(ssd_rt_val) else np.nan
            elif is_go_trial:
                if go_resp_empty and not _isnan(fix_resp) and float(fix_resp) != 0:
                    # late response in fix window: go_dur + fix_rt
                    go_rt_adj = (go_dur + fix_rt) / 1000.0 if (not _isnan(go_dur) and not _isnan(fix_rt)) else np.nan
                elif sst_type == 'NoRespGo':
                    # omission (no response anywhere): F3 will replace with max_go_rt
                    go_rt_adj = np.nan
                else:
                    go_rt_adj = go_rt / 1000.0
            else:
                go_rt_adj = np.nan

            # response_correct_go
            if is_stop_trial or is_zr or sst_type == 'STE':
                cgr = 'n/a'
            elif sst_type == 'NoRespGo':
                cgr = 'omission'
            else:
                cgr = 1 if sst_type == 'CorrectGo' else 0

            # correct_stop=1 if stop_resp, fix_resp, and ssd_resp are all empty;
            # n/a for go and ZeroSSDResp trials
            if is_go_trial or is_zr:
                cs = 'n/a'
            else:
                stop_resp_v = row.get('stop_resp', np.nan)
                ssd_resp_v  = row.get('ssd_resp', np.nan)
                _nr = lambda v: not _isnan(v) and float(v) != 0
                cs = 0 if (_nr(stop_resp_v) or _nr(fix_resp) or _nr(ssd_resp_v)) else 1

            # response_time: ssd_rt (real SSD-window RT) for STE/early_stop
            # -- checked first, since is_ste implies is_go_trial too, and
            # go_rt is empty for these trials (the response happened during
            # the SSD window, not the nominal go window); go_rt for other
            # go trials; stop_rt (StopSignal.RT) for incorrect_stop when
            # stop_resp was recorded -- gated on stop_resp presence, not
            # stop_rt==0, since stop_rt can be genuinely 0 for a real
            # stop-window response (timing glitch), matching response_time_go_adjusted's
            # go_resp-presence gate below.
            if is_ste or is_zr:
                rt_s = ssd_rt_val / 1000.0 if not _isnan(ssd_rt_val) else np.nan
            elif is_go_trial and not _isnan(go_rt):
                rt_s = go_rt / 1000.0
            elif (sst_type == 'IncorrectStop' and not _isnan(stop_rt)
                    and _nr(row.get('stop_resp', np.nan))):
                rt_s = stop_rt / 1000.0
            else:
                rt_s = np.nan

            # response_time_stop_adjusted: for any stop-type trial, used by F3 for
            #   violatorflag (MATLAB computes this for CorrectStop and
            #   IncorrectStop alike, gated on correct_stop, not sst_type)
            #   stop-window response (stop_resp present): stop_rt + ssd_dur
            #     -- gated on stop_resp presence, not stop_rt!=0, since
            #     stop_rt can be genuinely 0 (timing glitch) for a real
            #     stop-window response
            #   fix-window response (stop_resp absent, fix_resp present): stop_dur + fix_rt + ssd_dur
            if is_stop_trial:
                stop_dur_val = row.get('stop_dur', np.nan)
                has_stop_resp = _nr(stop_resp_v)
                has_fix_resp = not _isnan(fix_resp) and float(fix_resp) != 0
                if has_stop_resp and not _isnan(ssd_d):
                    stop_rt_base = stop_rt if not _isnan(stop_rt) else 0.0
                    stop_rt_adj_ms = stop_rt_base + ssd_d
                elif (not has_stop_resp) and has_fix_resp and not _isnan(fix_rt) and not _isnan(stop_dur_val) and not _isnan(ssd_d):
                    stop_rt_adj_ms = stop_dur_val + fix_rt + ssd_d
                else:
                    stop_rt_adj_ms = np.nan
                stop_rt_adj_s = _na_or(stop_rt_adj_ms / 1000.0 if not _isnan(stop_rt_adj_ms) else np.nan)
            else:
                stop_rt_adj_s = 'n/a'

            events.append({
                'onset': round((onset_ms - start_ms) / 1000.0, 3),
                'duration': _na_or(dur_ms / 1000.0 if not _isnan(dur_ms) else np.nan),
                'response_time': _na_or(rt_s),
                'trial_type': trial_type,
                'run_num': run_num,
                'response_correct_go': cgr,
                'correct_stop': cs,
                # duration_ssd: set for stop trials AND STE-recoded go trials so
                # compute_behav_sst can include STE SSDs in SSD_mean.
                'duration_ssd': _na_or(ssd_d / 1000.0 if not _isnan(ssd_d) else np.nan) if (is_stop_trial or is_zr or (is_go_trial and not _isnan(ssd_d) and ssd_d > 0)) else 'n/a',
                'response_time_go_adjusted': _na_or(go_rt_adj) if is_go_trial else 'n/a',
                'response_time_stop_adjusted': stop_rt_adj_s,
                # response_time_ssd: response time during SSD window; written for stop
                # trials and STE-recoded go trials (glitchflag needs it from
                # both).
                'response_time_ssd': _na_or(ssd_rt_val / 1000.0 if (is_stop_trial or is_ste or is_zr) and not _isnan(ssd_rt_val) else np.nan),
            })

    if not events:
        return None, timing_info, {}, 1, 'process_sst: no events produced'

    events_df = pd.DataFrame(events)
    meta = {
        'task': 'SST',
        'experiment_name': exp_name,
        'experiment_version': exp_ver,
        'runtime_version': runtime_ver,
        'behavioral_only': behavioral_only,
        'scanner_type': scanner_type,
        'switch_flag': switched,
        'runs_ok': runs_ok,
        'runs_invalid': runs_invalid,
    }
    return events_df, timing_info, meta, 0, ''


# ── task processor: REC ────────────────────────────────────────────────────

def _rec_switch_resp_format(df):
    """Map stim_cresp_block/stim_resp_block to numeric via the shared
    resp-name-map machinery (letter_fallback=False -- REC never uses
    r/g/b/y; comma_required is moot since REC never uses digit-brace
    format either). Returns (df, errcode, errmsg).
    """
    df = df.copy()
    name_map, errcode, errmsg = _build_resp_name_map(
        df['stim_cresp_block'].tolist(), comma_required=False, letter_fallback=False)
    if errcode:
        return df, errcode, errmsg
    for col in ['stim_cresp_block', 'stim_resp_block']:
        if col not in df.columns:
            continue
        df[col] = df[col].apply(lambda v: _apply_resp_map(v, name_map))
    return df, 0, ''


def process_rec(df_rec, outstem, switch_thresh=0.3, wm_stim_lookup=None, verbose=True):
    """Extract REC (recognition memory) events from REC DataFrame.

    Behavioral only; no scanner timing. wm_stim_lookup (stim basename ->
    '0_back'/'2_back', built from the raw WM E-Prime file independently of
    NBACK's own run-validity filtering) is embedded as the `block_type`
    column so Function 3 can compute the WM-condition breakdown without
    needing the NBACK events TSV.

    Returns (events_df, timing_info=None, meta, errcode, errmsg).
    events_df extra columns: response_accuracy, stim_type, stim, block_type
    """
    wm_stim_lookup = wm_stim_lookup or {}
    exp_name_rec, exp_ver_rec, runtime_ver_rec = _get_experiment_info(df_rec)

    df_rec, ec_switch, em_switch = _rec_switch_resp_format(df_rec)
    if ec_switch:
        return None, None, {}, ec_switch, 'process_rec: ' + em_switch

    for col in ['stim_acc_block', 'stim_rt_block']:
        if col in df_rec.columns:
            df_rec[col] = pd.to_numeric(df_rec[col], errors='coerce')

    # button switch on REC responses
    switched = False
    if 'stim_resp_block' in df_rec.columns and 'stim_cresp_block' in df_rec.columns:
        df_rec, switched = button_switch(df_rec, 'stim_resp_block', 'stim_cresp_block',
                                         switch_thresh=switch_thresh, verbose=verbose)
        if switched:
            # button_switch only swaps stim_resp_block; stim_acc_block must be
            # recomputed against the swapped response.
            resp = pd.to_numeric(df_rec['stim_resp_block'], errors='coerce')
            cresp = pd.to_numeric(df_rec['stim_cresp_block'], errors='coerce')
            valid = resp.notna() & cresp.notna()
            df_rec.loc[valid, 'stim_acc_block'] = (resp[valid] == cresp[valid]).astype(int)
            try:
                df_rec.to_csv(outstem + '_events_switched.csv', index=False)
            except Exception:
                pass

    # build events (onset n/a — behavioral only, no scanner timing)
    events = []
    for _, row in df_rec.iterrows():
        rt_ms = row.get('stim_rt_block', np.nan)
        acc = row.get('stim_acc_block', np.nan)
        stim_type_val = row.get('stim_type_block', '')
        stype_raw = '' if pd.isna(stim_type_val) else str(stim_type_val).strip().lower()
        stim_path = row.get('stim_block', '')
        stim_path_s = '' if pd.isna(stim_path) else str(stim_path).strip()
        stim_basename = os.path.basename(stim_path_s) if stim_path_s else 'n/a'
        events.append({
            'onset': 'n/a',
            'duration': 'n/a',
            'response_time': _na_or(rt_ms / 1000.0 if not _isnan(rt_ms) else np.nan),
            'trial_type': stype_raw if stype_raw else 'n/a',
            'run_num': 1,
            'response_accuracy': int(acc) if not _isnan(acc) else 'n/a',
            'stim_type': stype_raw if stype_raw else 'n/a',
            'stim': stim_basename,
            'block_type': wm_stim_lookup.get(stim_basename, 'n/a'),
        })

    if not events:
        return None, None, {}, 1, 'process_rec: no events produced'

    events_df = pd.DataFrame(events)
    meta = {
        'task': 'REC',
        'experiment_name': exp_name_rec,
        'experiment_version': exp_ver_rec,
        'runtime_version': runtime_ver_rec,
        'behavioral_only': True,
        'scanner_type': 'behavioral',
        'switch_flag': switched,
        'runs_ok': [1],
    }
    return events_df, None, meta, 0, ''


# ── output helpers ─────────────────────────────────────────────────────────

def _write_events_tsv(events_df, fname, verbose=True):
    """Write events DataFrame to TSV, using 'n/a' for missing values.

    Enforces BIDS column order (onset first, duration second) and sorts
    rows by onset ascending (n/a values sorted to end).
    """
    df = events_df.copy()

    # BIDS: onset must be column 1, duration must be column 2
    cols = list(df.columns)
    lead = [c for c in ('onset', 'duration') if c in cols]
    rest = [c for c in cols if c not in lead]
    df = df[lead + rest]

    # BIDS: rows must be sorted by onset; treat n/a as infinity (sort last)
    if 'onset' in df.columns:
        def _onset_key(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return float('inf')
        sort_keys = df['onset'].apply(_onset_key)
        df = df.iloc[sort_keys.argsort(kind='stable').values].reset_index(drop=True)

    # format numeric columns to 3 decimal places
    def _fmt(v):
        if v == 'n/a' or (not isinstance(v, str) and _isnan(v)):
            return 'n/a'
        try:
            return '{:.3f}'.format(float(v))
        except (TypeError, ValueError):
            return 'n/a'

    for col in ('onset', 'duration', 'response_time'):
        if col in df.columns:
            df[col] = df[col].apply(_fmt)

    # replace remaining NaN with 'n/a'
    df = df.fillna('n/a')

    df.to_csv(fname, sep='\t', index=False)
    if verbose:
        print('wrote {}'.format(fname))


def _write_events_json(json_data, fname, verbose=True):
    """Write events JSON sidecar."""
    with open(fname, 'w') as fh:
        json.dump(json_data, fh, indent=2)
    if verbose:
        print('wrote {}'.format(fname))


def _build_json_sidecar(task, error_flag=False, extra_flag=False):
    """Assemble the _events.json column-description sidecar dict.

    trial_type.Levels lists every trial type the task can write, with
    eprime_to_events' error_flag and extra_flag (NBACK only).
    """
    _special = {
        'dummy': 'Pre-task period covering discarded TRs',
        'scanner_trigger': 'Scanner trigger signal onset (negative onset = before acquisition start)',
        'task_init': 'Task initialization / first fixation onset',
    }
    _stim_desc = {
        'posface': 'positive face', 'neutface': 'neutral face',
        'negface': 'negative face', 'place': 'place/scene',
    }

    def _nback_levels():
        lev = dict(_special)
        lev['cue'] = 'Block cue'
        for n in ('0', '2'):
            for stim in _stim_desc:
                t = '{}_back_{}'.format(n, stim)
                lev[t] = '{}-back {}'.format(n, _stim_desc[stim])
                if error_flag:
                    lev[t + '_err'] = '{}-back {}, incorrect response'.format(
                        n, _stim_desc[stim])
            if extra_flag:
                for ttype in ('target', 'lure', 'nonlure'):
                    lev['{}_back_{}'.format(n, ttype)] = '{}-back {}'.format(n, ttype)
        return lev

    def _rec_levels():
        lev = {'dummy': 'Pre-task period covering discarded TRs'}
        lev.update(_rec_stim_levels())
        return lev

    def _rec_stim_levels():
        lev = {}
        for age in ('old', 'new'):
            for stim in _stim_desc:
                lev[age + stim] = '{} {}'.format(age.capitalize(), _stim_desc[stim])
        return lev

    task_levels = {
        'MID': dict(_special, **{
            'small_reward_antic': 'Small reward anticipation',
            'large_reward_antic': 'Large reward anticipation',
            'small_loss_antic': 'Small loss anticipation',
            'large_loss_antic': 'Large loss anticipation',
            'neutral_antic': 'Neutral anticipation',
            'small_reward_pos_feedback': 'Small reward positive feedback',
            'small_reward_neg_feedback': 'Small reward negative feedback',
            'large_reward_pos_feedback': 'Large reward positive feedback',
            'large_reward_neg_feedback': 'Large reward negative feedback',
            'small_loss_pos_feedback': 'Small loss positive feedback',
            'small_loss_neg_feedback': 'Small loss negative feedback',
            'large_loss_pos_feedback': 'Large loss positive feedback',
            'large_loss_neg_feedback': 'Large loss negative feedback',
            'neutral_pos_feedback': 'Neutral positive feedback',
            'neutral_neg_feedback': 'Neutral negative feedback',
        }),
        'NBACK': _nback_levels(),
        'SST': dict(_special, **{
            'correct_go': 'Correct go trial',
            'correctlate_go': 'Late correct go trial',
            'noresp_go': 'Go trial with no response',
            'incorrect_go': 'Incorrect go trial',
            'incorrectlate_go': 'Late incorrect go trial',
            'correct_stop': 'Successfully inhibited stop trial',
            'incorrect_stop': 'Failed stop trial',
            'early_stop': 'Stop trial with a response during the stop signal '
                          'delay, before the stop signal appeared',
            'zero_ssd_response': 'Stop trial with a stop signal delay of 0 and a '
                                 'response recorded during the SSD window',
        }),
        'REC': _rec_levels(),
    }

    levels = task_levels.get(task, dict(_special))

    expt_examples = {
        'MID': 'MID_20161218',
        'NBACK': 'Emo_Nback_20160913',
        'SST': 'SST_SP_20161208',
        'REC': 'Emo_Nback_RecMem_CDbox_20160815',
    }

    sidecar = {
        'onset': {'Description': 'Event onset relative to the first volume of '
                                 'the bold file, dummy volumes included (s); for '
                                 'a run with no scan (in beh), relative to the '
                                 'run start less the dummy time',
                  'Units': 's'},
        'duration': {'Description': 'Event duration (s)', 'Units': 's'},
        'response_time': {
            'Description': 'Response time relative to stimulus onset (s)',
            'Units': 's',
        },
        'trial_type': {
            'Description': 'Event condition label',
            'Levels': levels,
        },
        'run_num': {'Description': 'Run number within the E-Prime file'},
        'task_name': {
            'Description': 'E-Prime experiment name, for the task version and '
                           'variant (e.g., {})'.format(
                               expt_examples.get(task, 'MID_20161218')),
        },
        'task_version': {'Description': 'E-Prime experiment version string'},
        'runtime_version': {'Description': 'E-Prime runtime version string'},
        'switch_flag': {
            'Description': 'Whether the button mapping was swapped',
            'Levels': {'1': 'Button mapping swapped', '0': 'Button mapping not swapped'},
        },
        # Per-run timing QC. Written for every task; constant across all rows of
        # a run, since the value describes the run rather than the event. 'n/a'
        # when no timing info is available for the run (e.g. REC, which is
        # behavioral-only). Summarized over a session's runs by
        # events_to_behavioral() in cmig_abcd_events_to_behavioral.py.
        'start_delay': {
            'Description': 'Delay from the scanner trigger to the task start '
                           '(s); n/a if not measured (the behavioral version, '
                           'or a non-GE version with no scanner information)',
            'Units': 's',
        },
        'start_delay_diff': {
            'Description': 'start_delay minus the expected start delay for the '
                           'scanner type (s); n/a when start_delay is',
            'Units': 's',
        },
        'init_delay_diff': {
            'Description': 'Delay from the task start to the task '
                           'initialization, minus the expected delay (s); n/a '
                           'if unavailable',
            'Units': 's',
        },
        'start_time': {
            'Description': 'Run start time on the E-Prime computer '
                           '(YYYY-MM-DDThh:mm:ss); n/a if unavailable',
        },
        'timing_ok': {
            'Description': "Whether the run's timing checks gave no warnings; "
                           'n/a if no timing information is available or '
                           'start_delay is n/a',
            'Levels': {'1': 'No timing warnings', '0': 'One or more timing warnings'},
        },
    }

    if task == 'MID':
        sidecar['response_correct'] = {
            'Description': 'Response accuracy',
            'Levels': {'1': 'Correct', '0': 'Incorrect'},
        }
        sidecar['money'] = {
            'Description': 'Cumulative earnings at end of trial (dollars)',
            'Units': 'USD',
        }
        sidecar['response_early'] = {
            'Description': 'Whether the button was pressed before the target '
                           'appeared (too soon): 1 yes, 0 no; n/a except for '
                           'anticipation events',
            'Levels': {'1': 'Pressed too soon', '0': 'No premature press'},
        }
        sidecar['response_antic'] = {
            'Description': 'Whether any button press was recorded for the '
                           'trial, too soon or in a response window: 1 yes, 0 '
                           'no; n/a except for anticipation events',
            'Levels': {'1': 'Response recorded', '0': 'No response'},
        }
        sidecar['feedback_window'] = {
            'Description': 'Whether this MID version has the additional '
                           'response windows (TextDisplay, Feedback)',
            'Levels': {'1': 'Newer MID, with the additional response windows',
                       '0': 'Older MID, without them'},
        }
    elif task == 'NBACK':
        sidecar['response_accuracy'] = {
            'Description': 'Whether the response to the stimulus was correct; '
                           'n/a for cue and dummy events',
            'Levels': {'1': 'Correct', '0': 'Incorrect'},
        }
        sidecar['block_type'] = {
            'Description': 'Working memory block type',
            'Levels': {'0_back': '0-back block', '2_back': '2-back block'},
        }
        sidecar['stim_type'] = {
            'Description': 'Stimulus face/scene category',
            'Levels': {
                'posface': 'Positive face',
                'neutface': 'Neutral face',
                'negface': 'Negative face',
                'place': 'Place/scene',
            },
        }
        sidecar['target_type'] = {
            'Description': 'Target classification',
            'Levels': {
                'target': 'Target stimulus (match)',
                'lure': 'Lure stimulus (foil match)',
                'nonlure': 'Non-lure stimulus',
            },
        }
        sidecar['stim'] = {
            'Description': 'Stimulus image filename (basename, no path)',
        }
    elif task == 'SST':
        sidecar['response_correct_go'] = {
            'Description': 'Go trial response accuracy',
            'Levels': {'1': 'Correct', '0': 'Incorrect', 'omission': 'No response'},
        }
        sidecar['correct_stop'] = {
            'Description': 'Stop trial inhibition; n/a for go, early_stop and '
                           'zero_ssd_response trials',
            'Levels': {'1': 'Inhibited', '0': 'Failed'},
        }
        sidecar['duration_ssd'] = {
            'Description': 'Stop signal delay duration (s)',
            'Units': 's',
        }
        sidecar['response_time_go_adjusted'] = {
            'Description': 'Adjusted go reaction time (s)',
            'Units': 's',
        }
        sidecar['response_time_stop_adjusted'] = {
            'Description': 'Response time on incorrect_stop trials from the go '
                           'stimulus onset, stop signal delay included (s); '
                           'n/a for other trial types',
            'Units': 's',
        }
        sidecar['response_time_ssd'] = {
            'Description': 'Response time during the stop signal delay, before '
                           'the stop signal appeared (s): above 0 for early_stop '
                           'trials, 0 for other stop trials; n/a for go trials',
            'Units': 's',
        }
    elif task == 'REC':
        sidecar['response_accuracy'] = {
            'Description': 'Whether the response to the stimulus was correct',
            'Levels': {'1': 'Correct', '0': 'Incorrect'},
        }
        sidecar['stim_type'] = {
            'Description': 'Recognition stimulus category (e.g., oldposface, newneutface)',
            'Levels': _rec_stim_levels(),
        }
        sidecar['stim'] = {
            'Description': 'Stimulus image filename (basename, no path)',
        }
        sidecar['block_type'] = {
            'Description': 'n-back block type the stimulus was shown in during '
                           'the scanner task; n/a for stimuli not shown in the '
                           'n-back task (new recognition foils)',
            'Levels': {'0_back': '0-back block', '2_back': '2-back block'},
        }

    return sidecar


# ── no-responses warning ────────────────────────────────────────────────────

def _count_responses(df, task, available_fields):
    """Tally non-empty response fields across the whole session, using the
    same per-task field lists as MATLAB's check_no_responses() (MID:
    probe_resp/text_display_resp/feedback_resp, the latter two only when
    feedback_resp's raw column exists; SST: go_resp/fix_resp/stop_resp/
    ssd_resp; NBACK: stim_resp, restricted to rows with a real stimulus
    via the same block_type/stim_type/target_type filter as MATLAB's
    remove_non_events(); REC: stim_resp_block). Unlike MATLAB, a zero
    count here is surfaced as a warning, not a rejection -- Python still
    writes output for these sessions."""
    if task == 'MID':
        if 'feedback_resp' in available_fields:
            nresp = int(df['probe_resp'].notna().sum()
                       + df['text_display_resp'].notna().sum()
                       + df['feedback_resp'].notna().sum())
        else:
            nresp = int(df['probe_resp'].notna().sum())
    elif task == 'SST':
        nresp = int(df['go_resp'].notna().sum() + df['fix_resp'].notna().sum()
                   + df['stop_resp'].notna().sum() + df['ssd_resp'].notna().sum())
    elif task == 'NBACK':
        mask = df['block_type'].notna() & df['stim_type'].notna() & df['target_type'].notna()
        nresp = int(df.loc[mask, 'stim_resp'].notna().sum())
    elif task == 'REC':
        nresp = int(df['stim_resp_block'].notna().sum())
    else:
        return None
    return nresp


# ── public function ────────────────────────────────────────────────────────

def _run_start_times(df, task):
    """Wall-clock start time of each run, as 'YYYY-MM-DDThh:mm:ss'.

    As eprime_sprdsht_get.py computes eprime.csv's eprime_date/eprime_time:
    SessionDate + SessionTime of the first row, plus each run's start offset
    on E-Prime's clock, truncated to whole seconds. The offset is
    PrepTime.OnsetTime (MID) or BeginFix.StartTime (SST), every non-empty
    row in order, or CueFix.StartTime in the row after the last TRSyncPROC
    (run 1) or TRSyncPROCR2 (run 2) of Procedure[Block] (NBACK). Returns a
    list in E-Prime run order, empty when any of these is missing.
    """
    if task == 'MID':
        offsets = list(pd.to_numeric(df['initfix_onset'], errors='coerce').dropna())
    elif task == 'SST':
        offsets = list(pd.to_numeric(df['beginfix_start'], errors='coerce').dropna())
    elif task == 'NBACK':
        offsets = []
        blocks = df['block'].tolist()
        cuefix = pd.to_numeric(df['cuefix_start'], errors='coerce').tolist()
        for signal in ('TRSyncPROC', 'TRSyncPROCR2'):
            ind = [i for i, b in enumerate(blocks) if b == signal]
            if ind and ind[-1] + 1 < len(cuefix) and not math.isnan(cuefix[ind[-1] + 1]):
                offsets.append(cuefix[ind[-1] + 1])
    else:
        return []
    if not offsets:
        return []
    try:
        session_t0 = dateutil.parser.parse('{} {}'.format(
            str(df['session_date'].iloc[0]).strip(),
            str(df['session_time'].iloc[0]).strip()))
    except (ValueError, OverflowError, TypeError):
        return []
    return [(session_t0 + datetime.timedelta(seconds=int(int(o) / 1000)))
            .strftime('%Y-%m-%dT%H:%M:%S') for o in offsets]


def eprime_to_events(fname_eprime, outdir='.', outstem=None,
                     task=None, TR=0.8, skipTRs=None,
                     Manufacturer=None, SoftwareVersion=None,
                     run_offset=1, switch_thresh=0.3, error_flag=False,
                     extra_flag=False, ste_recode_flag=None,
                     fname_wm=None, run_numbers=None, run_times=None,
                     run_times_tolerance=4.0,
                     forceflag=False, verbose=True):
    """Read an ABCD E-Prime file and write BIDS events TSVs.

    Parameters
    ----------
    fname_eprime : str
        Path to E-Prime exported spreadsheet.
    outdir : str
        Output directory (created if absent) -- the BIDS func/ or beh/
        directory this run's files belong under; the caller is
        responsible for the sub-/ses- nesting above it. A non-standard
        'eprime' subdirectory is created under outdir for files that
        aren't official BIDS suffixes (currently just
        _events_switched.csv, when a button-swap was kept).
    outstem : str or None
        Base BIDS stem with no run/desc entity (e.g.
        sub-XXXXXXXX_ses-XXA_task-<task>); Function 1 appends
        _run-<NN> (and _desc-invalid, for a run with the wrong trial count or
        no responses)
        itself. Derived from fname_eprime if None (non-BIDS fallback).
    task : str or None
        Task the file is expected to hold ('MID', 'SST', 'NBACK', 'REC');
        an error if ExperimentName names another task, as in the MATLAB
        extractors. Auto-detected from ExperimentName if None.
    TR : float
        Repetition time in seconds (default 0.8).
    skipTRs : int or None
        Number of TRs to skip at run start. If None (default): 0 for a
        behavioral-only run (detected from ExperimentName); else computed
        from Manufacturer/SoftwareVersion via set_abcd_skipTRs() if given;
        else 0 with a warning (no information available).
    Manufacturer : str or None
        DICOM/BIDS Manufacturer field, used to compute skipTRs when
        skipTRs is None (ignored otherwise).
    SoftwareVersion : str or None
        DICOM SoftwareVersion(s) field, used to compute skipTRs when
        skipTRs is None (ignored otherwise).
    run_offset : int
        Starting run number for output file naming (default 1).
    switch_thresh : float
        Button-switch accuracy threshold (default 0.3).
    error_flag : bool
        Append '_err' to incorrect NBACK trial types (default False).
    extra_flag : bool
        Add target/lure/nonlure extra events for NBACK (default False).
    ste_recode_flag : bool or None
        SST only; controls how Stop-Too-Early (STE) trials are labeled.
        Default None: auto-detect from the SST experiment version --
        True if the first digit is >=2, else False. An explicit True or
        False always overrides the auto-detected value unconditionally
        (e.g. to compare both interpretations for the same session).
    fname_wm : str or None
        REC only: path to the WM (n-back) E-Prime file, read independently
        of any NBACK processing, to resolve each REC stimulus's 0-back/
        2-back condition into the `block_type` column. If None, `block_type`
        is 'n/a' for all rows.
    run_numbers : list of int or None
        Required for MID/SST/NBACK (ignored for REC, which is always a
        single behavioral administration with no BIDS run entity). One
        entry per run Function 1 is expected to detect internally, in
        chronological order -- matched strictly positionally, never
        reassigned. Each entry is either a positive BIDS run number
        (this detected run is written as that real run-NN) or 0 (this
        detected run is already known, by the caller/matching CSV, to
        have no corresponding scan -- dropped with a note, not written).
        A trailing 0 (or run of trailing 0s) may be omitted from the
        list entirely. If Function 1 detects fewer runs than
        len(run_numbers), the missing run(s) are warned about (nothing
        written for them); if it detects more, the extra run(s) are
        dropped with a warning.
    run_times : list of float or None
        Optional, same length as run_numbers if given (raises if not).
        Per-run expected E-Prime run timestamp (HHMMSS-style, matching
        eprime_time), used only as a cross-check against Function 1's
        own detected trigger-time gaps between runs -- catches a silent
        mislabeling if Function 1's own detected run count/order ever
        differs from what the caller assumed. Never used to resolve a
        position mismatch. Because eprime_time is a real wall-clock
        time and Function 1's own detected times are E-Prime-internal
        relative offsets, only the *gap* between consecutive runs can
        be cross-checked (not an absolute alignment) -- so this check
        is a no-op for a single-run match.
    run_times_tolerance : float
        Tolerance in seconds for the run_times gap cross-check (default
        4.0 -- eprime_time is whole-seconds precision, so two
        independently-rounded endpoints can differ from the true gap by
        up to ~1s each, needing a looser tolerance than a single
        timestamp comparison would; a GE run with a missed-trigger
        start_time correction (see check_eprime_timing) can add up to a
        few more seconds of legitimate deviation on top of that, since
        the corrected start_time need not land on the same whole-second
        boundary the caller's eprime_time metadata implies).
    forceflag : bool
        Overwrite existing output files (default False).
    verbose : bool
        Print progress (default True).

    Returns
    -------
    n_runs : int
        Count of trial-count-valid runs written.
    runs_ok : list of int
        BIDS run numbers written as valid (i.e. run_numbers values, not
        Function 1's internal run-detection order); for REC, always [1]
        (no BIDS run entity).
    errcode : int    0=success 1=error
    errmsg : str
    """
    os.makedirs(outdir, exist_ok=True)
    outdir_eprime = os.path.join(outdir, 'eprime')
    os.makedirs(outdir_eprime, exist_ok=True)

    if outstem is None:
        outstem = clean_fstem(fname_eprime)
    outstem_path = os.path.join(outdir, outstem)
    # non-standard files (not real BIDS suffixes) go in outdir/eprime/
    outstem_eprime_path = os.path.join(outdir_eprime, outstem)

    if run_times is not None and run_numbers is not None and len(run_times) != len(run_numbers):
        return 0, [], 1, ('eprime_to_events: run_times length ({}) does not match '
                          'run_numbers length ({})'.format(len(run_times), len(run_numbers)))

    # ── read E-Prime file ──────────────────────────────────────────────────
    # quick partial read first, to detect the task from ExperimentName
    # before choosing the full column set
    _peek_cols = ['ExperimentName', 'ExperimentVersion']
    _peek_fields = ['experiment_name', 'experiment_version']
    df_peek, ec, em, _ = read_eprime(fname_eprime, outstem_path,
                                      _peek_cols, _peek_fields, verbose=False)
    if ec:
        return 0, [], 1, em

    # auto-detect task from ExperimentName
    exp_name_raw = ''
    if df_peek is not None and 'experiment_name' in df_peek.columns:
        vals = df_peek['experiment_name'].dropna()
        if len(vals) > 0:
            exp_name_raw = str(vals.iloc[0])

    # a practice experiment is never a real run (eprime_sprdsht_get.py
    # rejects it too), whatever task its name contains
    if 'practice' in exp_name_raw.lower():
        return 0, [], 1, 'eprime_to_events: practice experiment, ExperimentName={!r}'.format(exp_name_raw)

    detected_task, behavioral_only, scanner_type = detect_task(exp_name_raw)
    if task is None:
        task = detected_task
    if task is None:
        return 0, [], 1, 'eprime_to_events: could not detect task from ExperimentName={!r}'.format(exp_name_raw)

    task = task.upper()
    # a name that names another task is a wrong file; one that names none
    #   (e.g. an nBack file eprime_sprdsht_get.py classified by its file name)
    #   is taken as the task passed
    if detected_task is None:
        print('eprime_to_events: WARNING: no task in ExperimentName={!r}; '
              'taking task {}'.format(exp_name_raw, task))
    elif detected_task != task:
        return 0, [], 1, ('eprime_to_events: wrong experiment for task {}: '
                          'ExperimentName={!r}'.format(task, exp_name_raw))

    # resolve skipTRs=None: behavioral-only runs need no dummy time;
    # otherwise compute from Manufacturer/SoftwareVersion if given, else
    # fall back to 0 with a warning (missing imaging data is not
    # special-cased here -- it isn't behavioral, so it takes this same path).
    skipTRs_unknown = False
    if skipTRs is None:
        if scanner_type == 'behavioral':
            skipTRs = 0
        elif Manufacturer is not None or SoftwareVersion is not None:
            skipTRs = set_abcd_skipTRs(Manufacturer, SoftwareVersion, verbose=verbose)
        else:
            skipTRs = 0
            skipTRs_unknown = True
            if verbose:
                print('eprime_to_events: WARNING: skipTRs not supplied and no '
                     'Manufacturer/SoftwareVersion given -- defaulting skipTRs=0')

    # the start delay is measured only by the GE version; elsewhere the start
    # time is constructed from the trigger, as trig + skipTRs*TR. With no
    # scanner (behavioral version) or no scanner info, that start delay and
    # the timing check on it say nothing, so they are written as n/a
    start_delay_na = behavioral_only or (scanner_type == 'Siemens' and skipTRs_unknown)

    # choose column set for full read
    if task == 'MID':
        cols, fields = _MID_COLS, _MID_FIELDS
    elif task == 'NBACK':
        cols, fields = _NBACK_COLS, _NBACK_FIELDS
    elif task == 'SST':
        cols, fields = _SST_COLS, _SST_FIELDS
    elif task == 'REC':
        cols, fields = _REC_COLS, _REC_FIELDS
    else:
        return 0, [], 1, 'eprime_to_events: unknown task {!r}'.format(task)

    df, ec, em, available_fields = read_eprime(
        fname_eprime, outstem_path, cols, fields, verbose=verbose)
    if ec:
        return 0, [], 1, em

    run_start_times = _run_start_times(df, task)

    nresp = _count_responses(df, task, available_fields)
    if nresp == 0:
        print('eprime_to_events: WARNING: no responses recorded anywhere '
             'in {} (task={})'.format(fname_eprime, task))

    if task == 'SST' and ste_recode_flag is None:
        _, exp_ver, _ = _get_experiment_info(df)
        try:
            ste_recode_flag = int(exp_ver.split('.')[0]) >= 2
        except (ValueError, IndexError):
            ste_recode_flag = False

    # ── process task ───────────────────────────────────────────────────────
    common_kw = dict(TR=TR, skipTRs=skipTRs, run_offset=run_offset, verbose=verbose)

    if task == 'MID':
        events_df, timing_info, meta, ec, em = process_mid(
            df, outstem_eprime_path, scanner_type, behavioral_only,
            available_fields=available_fields, **common_kw)

    elif task == 'NBACK':
        events_df, timing_info, meta, ec, em = process_nback(
            df, outstem_eprime_path, scanner_type, behavioral_only,
            switch_thresh=switch_thresh,
            error_flag=error_flag, extra_flag=extra_flag, **common_kw)

    elif task == 'SST':
        events_df, timing_info, meta, ec, em = process_sst(
            df, outstem_eprime_path, scanner_type, behavioral_only,
            switch_thresh=switch_thresh, ste_recode_flag=ste_recode_flag, **common_kw)

    elif task == 'REC':
        wm_stim_lookup = {}
        if fname_wm and os.path.exists(fname_wm):
            df_wm, ec_wm, em_wm, _ = read_eprime(
                fname_wm, outstem_eprime_path + '_WM', _WM_COLS, _WM_FIELDS, verbose=verbose)
            if ec_wm:
                if verbose:
                    print('eprime_to_events: WARNING: could not read WM file '
                          '{}: {}'.format(fname_wm, em_wm))
            else:
                wm_stim_lookup = _build_wm_stim_lookup(df_wm)
        events_df, timing_info, meta, ec, em = process_rec(
            df, outstem_eprime_path, wm_stim_lookup=wm_stim_lookup,
            switch_thresh=switch_thresh, verbose=verbose)

    if ec:
        return 0, [], ec, em

    runs_ok = meta.get('runs_ok', [])
    runs_invalid = meta.get('runs_invalid', [])
    all_detected = sorted(runs_ok + runs_invalid)

    if not all_detected:
        return 0, [], 1, 'eprime_to_events: no valid runs'

    # ── match internally-detected runs to caller-supplied BIDS run numbers ──
    if task == 'REC':
        # behavioral-only, always a single administration, no BIDS run
        # entity, no run_numbers/run_times matching.
        run_num_to_bids = {r: None for r in all_detected}
    else:
        if run_numbers is None:
            return 0, [], 1, ('eprime_to_events: run_numbers is required for '
                              'task={!r}'.format(task))

        # run_numbers and run_times are by E-Prime run, as the detected runs'
        #   numbers are
        n_expected = len(run_numbers)
        run_num_to_bids = {}
        for run_num in all_detected:
            if run_num > n_expected:
                if verbose:
                    print('eprime_to_events: WARNING: detected run {} has no '
                          'corresponding run_numbers entry -- dropped'.format(run_num))
                continue
            bids_rn = run_numbers[run_num - 1]
            if bids_rn == 0:
                if verbose:
                    print('eprime_to_events: run {} dropped (run_numbers[{}]=0, '
                          'known non-match)'.format(run_num, run_num - 1))
                continue
            run_num_to_bids[run_num] = bids_rn
        for i, bids_rn in enumerate(run_numbers):
            if bids_rn != 0 and (i + 1) not in all_detected:
                print('eprime_to_events: WARNING: expected run_numbers[{}]={} '
                      'but no run detected at that position -- nothing written '
                      'for it'.format(i, bids_rn))

        # run_times gap cross-check: verification only, never resolves a
        # position mismatch. eprime_time is wall-clock but Function 1's own
        # detected times are E-Prime-internal relative offsets, so only the
        # *gap* between consecutive matched runs can be compared -- a no-op
        # for a single matched run. Uses start_time_s (last trigger), not
        # trig_time_s (first trigger): eprime_time tracks true scan start,
        # and an atypical gap within a run's own trigger-wait sequence
        # shifts trig_time_s but not start_time_s.
        if run_times is not None and timing_info:
            _timing_by_internal_run = {ti.get('run'): ti for ti in timing_info}
            matched = sorted(r for r in run_num_to_bids if r <= len(run_times))
            for pa, pb in zip(matched[:-1], matched[1:]):
                t_a = _hhmmss_to_seconds(run_times[pa - 1])
                t_b = _hhmmss_to_seconds(run_times[pb - 1])
                ti_a = _timing_by_internal_run.get(pa)
                ti_b = _timing_by_internal_run.get(pb)
                if ti_a is None or ti_b is None or math.isnan(t_a) or math.isnan(t_b):
                    continue
                gap_expected = t_b - t_a
                gap_detected = (ti_b.get('start_time_s', float('nan')) -
                               ti_a.get('start_time_s', float('nan')))
                if math.isnan(gap_detected):
                    continue
                if abs(gap_expected - gap_detected) > run_times_tolerance:
                    return 0, [], 1, (
                        'eprime_to_events: run_times cross-check failed -- expected '
                        'gap between runs {} and {} is {:.3f}s but detected gap is '
                        '{:.3f}s (tolerance {:.3f}s); Function 1 may have detected '
                        'a different run count/order than the caller assumed'.format(
                            pa, pb, gap_expected,
                            gap_detected, run_times_tolerance))

    if not run_num_to_bids:
        return 0, [], 1, 'eprime_to_events: no runs to write after run_numbers matching'

    # a run with no responses is invalid as a wrong trial count is: its
    #   response-defined trial types (e.g. SST correct_go/noresp_go) would
    #   misrepresent it; the rule is the run selection's (has_responses)
    if task != 'REC':
        for run_num in sorted(run_num_to_bids):
            if run_num in runs_ok and not has_responses(
                    events_df[events_df['run_num'] == run_num].to_dict('records'), task):
                print('eprime_to_events: WARNING: run {} has no responses -- '
                      'writing it as desc-invalid'.format(run_num))
                runs_ok = [r for r in runs_ok if r != run_num]
                runs_invalid = sorted(runs_invalid + [run_num])

    n_runs = sum(1 for r in run_num_to_bids if r in runs_ok)
    task_lower = task.lower()

    # constant-valued metadata columns added to every row
    meta_cols = {
        'task_name': meta.get('experiment_name', ''),
        'task_version': meta.get('experiment_version', ''),
        'runtime_version': meta.get('runtime_version', ''),
        'switch_flag': int(meta.get('switch_flag', False)),
    }

    # build a lookup from run number to timing info entry
    timing_by_run = {}
    if timing_info:
        for ti in timing_info:
            timing_by_run[ti.get('run')] = ti

    # ── write per-run TSV files ────────────────────────────────────────────
    for run_num in sorted(run_num_to_bids):
        run_df = events_df[events_df['run_num'] == run_num].copy()

        # add metadata columns to all event rows
        for col, val in meta_cols.items():
            run_df[col] = val

        # imaging-task-only timing QC columns (n/a for REC, whose
        # timing_info is None) -- let events_to_behavioral() summarize them
        # from _events.tsv alone, without needing the E-Prime file or
        # re-deriving scanner-type/task-specific "expected delay" rules.
        run_ti = timing_by_run.get(run_num)
        if run_ti is not None and start_delay_na:
            run_df['start_delay'] = 'n/a'
            run_df['start_delay_diff'] = 'n/a'
            run_df['init_delay_diff'] = _na_or(run_ti.get('init_delay_diff_s', np.nan))
            run_df['timing_ok'] = 'n/a'
        elif run_ti is not None:
            run_df['start_delay'] = _na_or(run_ti.get('start_delay_s', np.nan))
            run_df['start_delay_diff'] = _na_or(run_ti.get('start_delay_diff_s', np.nan))
            run_df['init_delay_diff'] = _na_or(run_ti.get('init_delay_diff_s', np.nan))
            run_df['timing_ok'] = int(bool(run_ti.get('timing_ok', False)))
        else:
            run_df['start_delay'] = 'n/a'
            run_df['start_delay_diff'] = 'n/a'
            run_df['init_delay_diff'] = 'n/a'
            run_df['timing_ok'] = 'n/a'
        # wall-clock run start, compared with the bold acq_time downstream
        if 1 <= run_num <= len(run_start_times):
            run_df['start_time'] = run_start_times[run_num - 1]
        else:
            run_df['start_time'] = 'n/a'

        # onsets so far are from the end of the dummy volumes; BIDS times them
        # from the first stored volume, and ABCD stores its dummy volumes, so
        # shift by their duration (0 for REC and the behavioral version)
        skip_s = round(skipTRs * TR, 3)
        if skip_s:
            run_df['onset'] = (pd.to_numeric(run_df['onset'], errors='coerce')
                               + skip_s).round(3)

        # dummy row
        dummy = {
            'onset': 0.0,
            'duration': round(skipTRs * TR, 3),
            'response_time': 'n/a',
            'trial_type': 'dummy',
            'run_num': run_num,
        }
        dummy.update(meta_cols)
        for col in run_df.columns:
            if col not in dummy:
                dummy[col] = 'n/a'

        special_rows = [dummy]

        # scanner_trigger and task_init rows from timing info
        run_ti = timing_by_run.get(run_num)
        if run_ti is not None:
            start_delay = run_ti.get('start_delay_s', 0.0)
            init_delay = run_ti.get('init_delay_s', 0.0)

            trig_row = {
                'onset': round(skip_s - start_delay, 3),
                'duration': 0.0,
                'response_time': 'n/a',
                'trial_type': 'scanner_trigger',
                'run_num': run_num,
            }
            trig_row.update(meta_cols)
            for col in run_df.columns:
                if col not in trig_row:
                    trig_row[col] = 'n/a'

            init_row = {
                'onset': round(skip_s + init_delay, 3),
                'duration': 0.0,
                'response_time': 'n/a',
                'trial_type': 'task_init',
                'run_num': run_num,
            }
            init_row.update(meta_cols)
            for col in run_df.columns:
                if col not in init_row:
                    init_row[col] = 'n/a'

            special_rows.extend([trig_row, init_row])

        special_df = pd.DataFrame(special_rows)
        run_df = pd.concat([special_df, run_df], ignore_index=True, sort=False)

        # REC is behavioral-only and always has exactly one run: no run
        # entity at all. Other tasks use the caller-assigned BIDS run
        # number, with a desc-invalid entity for a run with the wrong trial
        # count or no responses.
        if task == 'REC':
            tsv_name = '{}_events.tsv'.format(outstem)
        else:
            bids_rn = run_num_to_bids[run_num]
            if run_num in runs_invalid:
                tsv_name = '{}_run-{:02d}_desc-invalid_events.tsv'.format(outstem, bids_rn)
            else:
                tsv_name = '{}_run-{:02d}_events.tsv'.format(outstem, bids_rn)
        tsv_path = os.path.join(outdir, tsv_name)

        if not forceflag and os.path.exists(tsv_path):
            if verbose:
                print('eprime_to_events: {} exists, skipping (forceflag=False)'.format(tsv_path))
            continue

        _write_events_tsv(run_df, tsv_path, verbose=verbose)

    # ── remove stale run TSVs from previous runs no longer written ────────────
    # (e.g. a run dropped by a new run_numbers mapping, or one that flipped
    # between valid/desc-invalid across reruns). Only ever deletes a run
    # number that's a candidate value in THIS call's own run_numbers
    # (candidate_bids) -- a subject split across multiple raw files gets one
    # eprime_to_events() call per file, each with its own run_numbers slice,
    # so a run number outside that slice could belong to a sibling file's
    # current output rather than being genuinely stale.
    if task != 'REC' and forceflag and run_num_to_bids:
        import glob as _glob
        written_valid_bids = set(run_num_to_bids[r] for r in run_num_to_bids if r in runs_ok)
        written_invalid_bids = set(run_num_to_bids[r] for r in run_num_to_bids if r in runs_invalid)
        candidate_bids = set(rn for rn in (run_numbers or []) if rn)
        pattern = os.path.join(outdir, '{}_run-*_events.tsv'.format(outstem))
        for stale in _glob.glob(pattern):
            base = os.path.basename(stale)
            is_invalid_variant = base.endswith('_desc-invalid_events.tsv')
            if is_invalid_variant:
                num_str = base[len(outstem + '_run-'):-len('_desc-invalid_events.tsv')]
            else:
                num_str = base[len(outstem + '_run-'):-len('_events.tsv')]
            try:
                bids_n = int(num_str)
            except ValueError:
                continue
            if bids_n not in candidate_bids:
                continue
            keep = ((is_invalid_variant and bids_n in written_invalid_bids) or
                    (not is_invalid_variant and bids_n in written_valid_bids))
            if not keep:
                os.remove(stale)
                if verbose:
                    print('eprime_to_events: removed stale {}'.format(base))

    if task == 'REC':
        return n_runs, runs_ok, 0, ''
    runs_ok_bids = sorted(run_num_to_bids[r] for r in run_num_to_bids if r in runs_ok)
    return n_runs, runs_ok_bids, 0, ''


# ── root sidecars ──────────────────────────────────────────────────────────

# task label in the BIDS file names, as cmig_convert_abcd_eprime writes them
_TASK_LABELS = {'MID': 'mid', 'NBACK': 'nback', 'SST': 'sst', 'REC': 'rec'}


def write_root_sidecars(rootdir, tasks=None, error_flag=False,
                        extra_flag=False, forceflag=False, verbose=True):
    """Write task-<task>_events.json for each task into rootdir.

    At the rawdata root, each describes the columns of every _events.tsv of
    its task, func/ and beh/, by BIDS inheritance; any other directory gets
    the same files (e.g. to send them). tasks: any of 'MID', 'NBACK', 'SST',
    'REC' (default all); error_flag and extra_flag as given to
    eprime_to_events for NBACK.

    Returns (errcode, errmsg).
    """
    tasks = tasks or list(_TASK_LABELS)
    try:
        os.makedirs(rootdir, exist_ok=True)
        for task in tasks:
            task = task.upper()
            if task not in _TASK_LABELS:
                return 1, 'write_root_sidecars: unknown task {}'.format(task)
            fname = os.path.join(rootdir,
                                 'task-{}_events.json'.format(_TASK_LABELS[task]))
            if os.path.exists(fname) and not forceflag:
                if verbose:
                    print('write_root_sidecars: {} exists, skipping '
                          '(forceflag=False)'.format(fname))
                continue
            sidecar = _build_json_sidecar(task, error_flag=error_flag,
                                          extra_flag=extra_flag)
            # written whole, then renamed, so a reader never sees a partial file
            fname_tmp = fname + '.tmp{}'.format(os.getpid())
            _write_events_json(sidecar, fname_tmp, verbose=False)
            os.replace(fname_tmp, fname)
            if verbose:
                print('wrote {}'.format(fname))
    except OSError as err:
        return 1, 'write_root_sidecars: {}'.format(err)
    return 0, ''


def parse_args():
    parser = argparse.ArgumentParser(
        description='Write the task-<task>_events.json sidecars, describing '
                    'the columns of the E-Prime events files, into a BIDS '
                    'rawdata root (or any directory).')
    parser.add_argument('--rootdir', required=True,
                        help='output directory, e.g. the rawdata root')
    parser.add_argument('--task', nargs='+', default=list(_TASK_LABELS),
                        choices=list(_TASK_LABELS),
                        help='tasks (default: all)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='overwrite existing files (default: 0)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1])
    return parser.parse_args()


def main():
    args = parse_args()
    errcode, errmsg = write_root_sidecars(
        args.rootdir, tasks=args.task, forceflag=bool(args.forceflag),
        verbose=bool(args.verbose))
    if errcode:
        print('ERROR: ' + errmsg, file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
