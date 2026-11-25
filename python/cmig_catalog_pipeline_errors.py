#!/usr/bin/env python3
"""
Keep a cumulative record of processing errors by session and step, across
batches, reruns and rounds.

For each session given, every status file the script is expected to write
(cmig_mkjobs_bids_pipeline.status_rows) is read. A status row with errcode
!= 0, or errcode 0 with an errmsg (a warning), is an error occurrence. Where
the job has ended, a missing expected status file is one too: the first
missing step gets the job's .err (a crash, with its MATLAB signature), and
any later one "no status file: not run after an earlier step failed". Each
occurrence is classed by the first matching rule of the rules file
(errors_rules.tsv: rule_id, error_class, module glob, pattern, note; a
python regex searched in the errmsg) and merged into the catalog:

  errors_log.tsv      every occurrence, appended; a repeat of the same step
                      and errmsg_norm updates only its n_seen and last_seen
  errors_current.tsv  one row per (participant_id, session_id, step) that
                      ever failed, with its latest error and whether a later
                      run of the step succeeded (resolved, resolved_date,
                      resolved_by, resolution rerun_ok); hand notes from
                      errors_notes.tsv (read, never written) merged in
  errors_log.json, errors_current.json   the column descriptions
  history/            weekly copies of the two tables,
                      errors_log_<YYMMDD>.tsv and errors_current_<YYMMDD>.tsv,
                      made when the newest copy's date stamp is 7 or more
                      days old; never pruned

Classes: not_processed (a refusal, fixed by the inputs), missing_input,
failure, file_read, infrastructure, resource (out of memory, disk full,
quota, I/O error, stale file handle), warning (errcode 0 with an errmsg) and
unknown (no rule matched). Successful steps are not recorded. The rounds
(cmig_run_bids_rounds.py) retry a session whose status errmsg matches a
resource rule, even with errcode 0 (class_regex).

Sessions come from one or more mkjobs batchdirs (--batchdir: scriptlist.txt,
the job scripts for the script name, roots and options, and pbsout/ for the
.out/.err), or from --participant_id/--session_id or --seslist_file with the
script name and roots given as options (a backfill of an existing tree,
where every job is taken as ended).

The catalog directory is locked (fcntl) while it is updated, and each table
is written to .new and renamed, so a round and a hand run do not interleave.

Usage:
  cmig_catalog_pipeline_errors.py --catalog_dir DIR --batchdir DIR [DIR ...]
  cmig_catalog_pipeline_errors.py --catalog_dir DIR --seslist_file TSV \\
      --script_name cmig_pipeline_bids --raw_rootdir ... --proc_rootdir ...

Created:  11/06/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Last mod: 11/20/2026 by Don Hagler
"""

import argparse
import csv
import datetime
import fcntl
import fnmatch
import json
import os
import re
import shutil
import socket
import sys

import cmig_check_batchdir as cbd
import cmig_mkjobs_bids_pipeline as mkjobs
import cmig_sing_utils

TIME_FORMAT = '%Y-%m-%dT%H:%M:%S'
NA = 'n/a'

# days between the history/ copies of the tables
HISTORY_DAYS = 7
HISTORY_TABLES = ['errors_log', 'errors_current']

CLASSES = ['not_processed', 'missing_input', 'failure', 'file_read',
           'infrastructure', 'resource', 'warning', 'unknown']

# the errmsg of a missing status file, by why it is missing; the rules class
#   them (mi_not_run, in_lost)
MSG_NOT_RUN = 'no status file: not run after an earlier step failed'
MSG_LOST = 'no status file'

ROOT_ARGS = ['raw_rootdir', 'proc_rootdir', 'fsurf_rootdir',
             'analysis_rootdir', 'postqc_rootdir', 'rawqc_rootdir',
             'pc_rootdir', 'manualqc_rootdir']
# job script options that change the expected status files
OPT_ARGS = ['func_flag', 'dwi_flag', 'task']

LOG_COLUMNS = [
    'participant_id', 'session_id', 'modality', 'datatype', 'module', 'task',
    'acq', 'step', 'errcode', 'errmsg', 'errmsg_norm', 'error_class',
    'error_rule', 'err_signature', 'source', 'startdate', 'procdate',
    'first_seen', 'last_seen', 'n_seen', 'batchdir', 'jobname', 'jobid',
    'host', 'app_name', 'app_version', 'software_version', 'fs_version']
CURRENT_COLUMNS = LOG_COLUMNS + ['resolved', 'resolved_date', 'resolved_by',
                                 'resolution', 'note']
NOTES_COLUMNS = ['participant_id', 'session_id', 'step', 'errmsg_norm',
                 'note', 'resolution']

DESCRIPTIONS = {
    'participant_id': 'Participant identifier',
    'session_id': 'Session identifier',
    'modality': 'sMRI, fMRI or dMRI; none for a session-level step',
    'datatype': 'The status file\'s directory: anat, func, dwi, fmap, fsurf, '
                'beh, or session',
    'module': 'rawqc, proc, fsurf, analyze, postqc, behav, pc or pcqc',
    'task': 'Task label from the status file name',
    'acq': 'acq label from the status file name',
    'step': 'datatype/status file stem, e.g. anat/desc-proc',
    'errcode': 'errcode in the status file; n/a where there is none',
    'errmsg': 'errmsg in the status file, else the job\'s .err',
    'errmsg_norm': 'errmsg with paths, sub-/ses-/run- labels and numbers '
                   'replaced, or the MATLAB error signature of a crash; '
                   'the grouping key',
    'error_class': 'not_processed, missing_input, failure, file_read, '
                   'infrastructure, resource, warning or unknown',
    'error_rule': 'rule_id of the rule that set error_class',
    'err_signature': 'MATLAB error identifier and innermost frame in the .err',
    'source': 'status (a status file), err (the job\'s .err, no status file) '
              'or missing (no status file, no .err)',
    'startdate': 'startdate in the status file',
    'procdate': 'procdate in the status file',
    'first_seen': 'When the catalog first recorded this error',
    'last_seen': 'When the catalog last saw this error',
    'n_seen': 'Runs of the step that gave this error (distinct procdate)',
    'batchdir': 'Batch directory of the job',
    'jobname': 'Job name',
    'jobid': 'Scheduler job id (JOB_ID in the .out)',
    'host': 'Host the job ran on (HOSTNAME in the .out)',
    'app_name': 'app_name in the status file',
    'app_version': 'app_version in the status file',
    'software_version': 'software_version in the status file',
    'fs_version': 'fs_version in the status file',
    'resolved': '1 where a later run of the step succeeded',
    'resolved_date': 'procdate of that successful run',
    'resolved_by': 'software_version of that successful run',
    'resolution': 'rerun_ok, or a resolution given by hand in errors_notes.tsv',
    'note': 'note given by hand in errors_notes.tsv',
}


def _now():
    return datetime.datetime.now().strftime(TIME_FORMAT)


def _na(val):
    if val is None:
        return NA
    val = str(val).replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
    return val if val.strip() else NA


def _val(row, key):
    """A table value, '' for n/a or absent."""
    val = row.get(key)
    return '' if val is None or val == NA else val


#############################################################################
# rules and classes
#############################################################################

def default_rules_file():
    return os.path.join(cmig_sing_utils.cmig_getenv('CMIG_PIPELINE_DIR'),
                        'parms', 'errors_rules.tsv')


def load_rules(fname=None):
    """The rules, in file order, each with its pattern compiled."""
    fname = fname or default_rules_file()
    rules = []
    with open(fname, newline='') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if not row.get('rule_id') or row['rule_id'].startswith('#'):
                continue
            if row['error_class'] not in CLASSES:
                raise ValueError('{0}: rule {1}: unknown error_class {2}'.format(
                    fname, row['rule_id'], row['error_class']))
            row['regex'] = re.compile(row['pattern'], re.MULTILINE)
            rules.append(row)
    return rules


def class_regex(rules, error_class, rule_ids=None):
    """One compiled regex joining the patterns of the rules of a class (or
    of those of rule_ids), or None where there are none."""
    pats = [r['pattern'] for r in rules if r['error_class'] == error_class and
            (rule_ids is None or r['rule_id'] in rule_ids)]
    if not pats:
        return None
    return re.compile('|'.join('(?:{0})'.format(p) for p in pats), re.MULTILINE)


def classify(errmsg, module, rules, errcode=None):
    """(error_class, rule_id) of one error: warning for errcode 0, else the
    first rule whose module glob matches and whose pattern is found in the
    errmsg, else unknown."""
    if errcode is not None and str(errcode) == '0':
        return 'warning', 'errcode0'
    for rule in rules:
        if not fnmatch.fnmatch(module or '', rule.get('module') or '*'):
            continue
        if rule['regex'].search(errmsg or ''):
            return rule['error_class'], rule['rule_id']
    return 'unknown', NA


def normalize(errmsg):
    """errmsg with absolute paths reduced to <path>/<basename>, sub-/ses-/run-
    labels and numbers replaced and spaces collapsed, so two sessions with the
    same fault share it."""
    text = errmsg or ''
    text = re.sub(r'(?<![\w.])/[^\s:;,\'"]+',
                  lambda m: '<path>/' + os.path.basename(m.group(0)), text)
    text = re.sub(r'sub-[A-Za-z0-9]+', 'sub-<id>', text)
    text = re.sub(r'ses-[A-Za-z0-9]+', 'ses-<id>', text)
    text = re.sub(r'run-\d+', 'run-<n>', text)
    text = re.sub(r'(?<![\w<])\d+(\.\d+)?', '<n>', text)
    return re.sub(r'\s+', ' ', text).strip()


#############################################################################
# steps
#############################################################################

def step_fields(step, task=''):
    """modality, datatype, module, task and acq of a step name."""
    datatype, stem = step.split('/', 1)
    ents = dict(e.split('-', 1) for e in stem.split('_') if '-' in e)
    module = ents.get('desc', '')
    acq = ents.get('acq', '')
    if datatype in ('anat', 'fsurf'):
        modality = 'sMRI'
    elif datatype in ('func', 'beh'):
        modality = 'fMRI'
    elif datatype == 'dwi':
        modality = 'dMRI'
    elif datatype == 'fmap':
        modality = {'func': 'fMRI', 'dwi': 'dMRI'}.get(acq, NA)
    else:
        modality = 'none'
    return {'modality': modality, 'datatype': datatype, 'module': module,
            'task': task or ents.get('task', ''), 'acq': acq}


#############################################################################
# sessions to read
#############################################################################

def parse_job_script(fpath):
    """(script_name, participant_id, session_id, dirs, opts) from an mkjobs
    job script, or None where it names no cmig script."""
    try:
        with open(fpath) as f:
            text = f.read()
    except (IOError, OSError):
        return None
    m = re.search(r'\$pydir/(\w+)\.py(?:\s+(\w+))?', text)
    if not m:
        return None
    script_name = m.group(1)
    if script_name == 'cmig_run_bids_module' and m.group(2):
        script_name = m.group(2)
    args = {}
    for line in text.splitlines():
        line = line.strip().rstrip('\\').strip()
        if not line.startswith('--'):
            continue
        parts = line[2:].split()
        if parts:
            args[parts[0]] = parts[1:]
    pid = (args.get('participant_id') or [''])[0]
    sid = (args.get('session_id') or [''])[0]
    dirs = dict((k, v[0]) for k, v in args.items() if k.endswith('_rootdir') and v)
    # a module script names its output root output_rootdir
    if 'output_rootdir' in args and script_name in mkjobs.SCRIPT_CONFIG:
        key = mkjobs.SCRIPT_CONFIG[script_name].get('output_rootdir_from')
        if key:
            dirs[key] = args['output_rootdir'][0]
    opts = {}
    for key in OPT_ARGS:
        if key in args and args[key]:
            vals = args[key]
            if key == 'task':
                opts[key] = vals
            else:
                opts[key] = int(vals[0]) if vals[0].lstrip('-').isdigit() else vals[0]
    return script_name, pid, sid, dirs, opts


def batch_sessions(batchdir):
    """One dict per job of an mkjobs batchdir: the session, its script, roots
    and options, and the job's .out/.err and whether it ended."""
    pbsout = os.path.join(batchdir, 'pbsout')
    entries = os.listdir(pbsout) if os.path.isdir(pbsout) else []
    has_markers = cbd.markers_present(pbsout)
    sessions = []
    for job in cbd.read_scriptlist(batchdir):
        parsed = parse_job_script(os.path.join(batchdir, job + '.csh'))
        if parsed is None:
            print('WARNING: no cmig script in {0}'.format(job), file=sys.stderr)
            continue
        script_name, pid, sid, dirs, opts = parsed
        out_path, err_path = cbd.resolve_job_files(pbsout, job, entries)
        state = cbd.job_state(pbsout, job, out_path, err_path, has_markers)
        err_text = cbd.read_file(err_path) if err_path else ''
        # a job killed by the scheduler reads as started, its marker never
        #   renamed; its .err says it ended
        ended = (state == cbd.STATE_FINISHED or
                 (state == cbd.STATE_STARTED and 'CANCELLED' in err_text))
        header = cbd.out_header(out_path)
        sessions.append({
            'participant_id': pid, 'session_id': sid,
            'script_name': script_name, 'dirs': dirs, 'opts': opts,
            'ended': ended, 'err_text': err_text,
            'batchdir': batchdir, 'jobname': job,
            'jobid': header.get('JOB_ID', ''), 'host': header.get('HOSTNAME', '')})
    return sessions


def listed_sessions(args):
    """One dict per session given by --participant_id/--session_id or
    --seslist_file, with the roots and options given; every job is taken as
    ended."""
    pairs = []
    if args.seslist_file:
        with open(args.seslist_file, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                pairs.append((row['participant_id'], row['session_id']))
    if args.participant_id and args.session_id:
        pairs.append((args.participant_id, args.session_id))
    dirs = dict((k, getattr(args, k)) for k in ROOT_ARGS if getattr(args, k))
    opts = {'func_flag': args.func_flag, 'dwi_flag': args.dwi_flag}
    if args.task:
        opts['task'] = args.task
    return [{'participant_id': pid, 'session_id': sid,
             'script_name': args.script_name, 'dirs': dirs, 'opts': opts,
             'ended': bool(args.ended_flag), 'err_text': '',
             'batchdir': '', 'jobname': '', 'jobid': '', 'host': ''}
            for pid, sid in pairs]


#############################################################################
# occurrences
#############################################################################

def err_crash(err_text):
    """(errmsg, err_signature) of a job's .err, with benign lines and the
    wrapper's "status file not found" lines dropped; ('', '') where nothing
    is left."""
    lines = [l for l in cbd.strip_benign(err_text).splitlines()
             if 'status file not found' not in l]
    text = '\n'.join(lines).strip()
    if not text:
        return '', ''
    sig = cbd.err_signature(text)
    if sig == 'no MATLAB error in .err':
        sig = ''
    return cbd.tsv_errmsg(text), sig


def session_occurrences(ses, rules):
    """(occurrences, successes) of one session: a row dict per error, and
    {step: status row} for each step whose status is errcode 0 with no
    errmsg."""
    rows = mkjobs.status_rows(ses['script_name'], ses['participant_id'],
                              ses['session_id'], ses['dirs'], ses['opts'])
    occs = []
    successes = {}
    any_failed = any(r['status'] is not None and
                     str(r['status'].get('errcode', '')) != '0' for r in rows)
    crash_used = False
    for r in rows:
        st = r['status']
        fields = step_fields(r['step'], r['task'])
        if st is not None:
            errcode = str(st.get('errcode', ''))
            errmsg = st.get('errmsg', '')
            if errcode == '0' and not errmsg.strip():
                successes[r['step']] = st
                continue
            source, sig = 'status', ''
        elif ses['ended']:
            errcode = NA
            errmsg, sig = ('', '') if crash_used else err_crash(ses['err_text'])
            if errmsg:
                source = 'err'
                crash_used = True
            elif any_failed or crash_used:
                source, errmsg = 'missing', MSG_NOT_RUN
            else:
                source, errmsg = 'missing', MSG_LOST
            st = {}
        else:
            continue
        error_class, rule = classify(errmsg, fields['module'], rules,
                                     errcode if errcode != NA else None)
        occ = {'participant_id': ses['participant_id'],
               'session_id': ses['session_id'], 'step': r['step'],
               'errcode': errcode, 'errmsg': errmsg,
               'errmsg_norm': sig or normalize(errmsg),
               'error_class': error_class, 'error_rule': rule,
               'err_signature': sig, 'source': source,
               'batchdir': ses['batchdir'], 'jobname': ses['jobname'],
               'jobid': ses['jobid'], 'host': ses['host']}
        occ.update(fields)
        for key in ('startdate', 'procdate', 'app_name', 'app_version',
                    'software_version', 'fs_version'):
            occ[key] = st.get(key, '')
        occs.append(occ)
    return occs, successes


#############################################################################
# the catalog tables
#############################################################################

def read_table(fname):
    if not os.path.isfile(fname):
        return []
    with open(fname, newline='') as f:
        return [dict(r) for r in csv.DictReader(f, delimiter='\t')]


def write_table(fname, rows, columns):
    rows = sorted(rows, key=lambda r: (_val(r, 'participant_id'),
                                       _val(r, 'session_id'), _val(r, 'step'),
                                       _val(r, 'first_seen')))
    tmp = fname + '.new'
    with open(tmp, 'w', newline='') as f:
        f.write('\t'.join(columns) + '\n')
        for row in rows:
            f.write('\t'.join(_na(row.get(c)) for c in columns) + '\n')
    os.rename(tmp, fname)


def write_json(fname, columns):
    tmp = fname + '.new'
    with open(tmp, 'w') as f:
        json.dump(dict((c, {'Description': DESCRIPTIONS.get(c, '')})
                       for c in columns), f, indent=2)
        f.write('\n')
    os.rename(tmp, fname)


def take_lock(catalog_dir):
    """The open lock file, locked, waiting for another writer to finish."""
    fid = open(os.path.join(catalog_dir, 'errors.lock'), 'a+')
    fcntl.flock(fid, fcntl.LOCK_EX)
    fid.seek(0)
    fid.truncate()
    fid.write('pid {0} host {1} start {2}\n'.format(
        os.getpid(), socket.gethostname(), _now()))
    fid.flush()
    return fid


def copy_history(catalog_dir, today=None):
    """Copy the tables to history/<table>_<YYMMDD>.tsv where the newest copy's
    date stamp is HISTORY_DAYS or more days old; nothing is deleted. Returns
    the stamp of the copies made, or ''."""
    today = today or datetime.date.today()
    histdir = os.path.join(catalog_dir, 'history')
    pat = re.compile(r'^(?:{0})_(\d{{6}})\.tsv$'.format('|'.join(HISTORY_TABLES)))
    stamps = []
    if os.path.isdir(histdir):
        for fname in os.listdir(histdir):
            m = pat.match(fname)
            if m:
                try:
                    stamps.append(datetime.datetime.strptime(m.group(1),
                                                             '%y%m%d').date())
                except ValueError:
                    continue
    if stamps and (today - max(stamps)).days < HISTORY_DAYS:
        return ''
    os.makedirs(histdir, exist_ok=True)
    stamp = today.strftime('%y%m%d')
    for table in HISTORY_TABLES:
        src = os.path.join(catalog_dir, table + '.tsv')
        dst = os.path.join(histdir, '{0}_{1}.tsv'.format(table, stamp))
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)
    return stamp


def merge(log_rows, current, occs, successes, now):
    """Fold one session's occurrences and successes into the log rows and the
    current rows ({(pid, sid, step): row}); returns the counts (new,
    repeated, resolved)."""
    n_new = n_repeat = n_resolved = 0
    for occ in occs:
        key = (occ['participant_id'], occ['session_id'], occ['step'])
        cur = current.get(key)
        if (cur is not None and _val(cur, 'resolved') != '1' and
                _val(cur, 'errmsg_norm') == (occ['errmsg_norm'] or '')):
            # the same open error: a new run of the step counts once more
            if _val(cur, 'procdate') != (occ['procdate'] or ''):
                cur['n_seen'] = str(int(_val(cur, 'n_seen') or 0) + 1)
                cur['procdate'] = occ['procdate']
            cur['last_seen'] = now
            for row in reversed(log_rows):
                if ((row['participant_id'], row['session_id'], row['step']) == key
                        and _val(row, 'errmsg_norm') == _val(cur, 'errmsg_norm')):
                    row['n_seen'] = cur['n_seen']
                    row['last_seen'] = now
                    break
            n_repeat += 1
            continue
        row = dict(occ, first_seen=now, last_seen=now, n_seen='1')
        log_rows.append(row)
        current[key] = dict(row, resolved='0', resolved_date='',
                            resolved_by='', resolution='', note='')
        n_new += 1
    for step, st in successes.items():
        key = (st.get('participant_id', ''), st.get('session_id', ''), step)
        cur = current.get(key)
        if cur is None or _val(cur, 'resolved') == '1':
            continue
        if (st.get('procdate') or '') > _val(cur, 'procdate'):
            cur.update({'resolved': '1', 'resolved_date': st.get('procdate', ''),
                        'resolved_by': st.get('software_version', ''),
                        'resolution': 'rerun_ok'})
            n_resolved += 1
    return n_new, n_repeat, n_resolved


def merge_notes(current, notes):
    """Hand notes: a row applies to the current rows that match each of its
    keys given, the session's step (participant_id, session_id, step) and the
    errmsg_norm; a row naming only an errmsg_norm applies to every session
    with it."""
    for note in notes:
        step = _val(note, 'step')
        norm = _val(note, 'errmsg_norm')
        if not step and not norm:
            continue
        for key, cur in current.items():
            if step and key != (note.get('participant_id'),
                                note.get('session_id'), step):
                continue
            if norm and _val(cur, 'errmsg_norm') != norm:
                continue
            if _val(note, 'note'):
                cur['note'] = note['note']
            if _val(note, 'resolution'):
                cur['resolution'] = note['resolution']


def catalog_errors(catalog_dir, sessions, rules_file=None, verbose=True):
    """Read the sessions' status files and jobs and update the catalog in
    catalog_dir. Returns (errcode, errmsg)."""
    try:
        rules = load_rules(rules_file)
    except (IOError, OSError, ValueError, re.error) as e:
        return 1, 'rules not loaded: {0}'.format(e)
    os.makedirs(catalog_dir, exist_ok=True)
    now = _now()
    # read everything before taking the lock, so the lock is held briefly
    found = []
    for ses in sessions:
        # a tree that has been removed would read as every step missing
        gone = sorted(set(d for d in ses['dirs'].values()
                          if not os.path.isdir(d)))
        if gone:
            print('WARNING: {0} {1}: skipped, root not found: {2}'.format(
                ses['participant_id'], ses['session_id'], ', '.join(gone)),
                file=sys.stderr)
            continue
        try:
            found.append(session_occurrences(ses, rules))
        except (IOError, OSError, KeyError, ValueError) as e:
            print('WARNING: {0} {1}: not read: {2}'.format(
                ses['participant_id'], ses['session_id'], e), file=sys.stderr)
    fid = take_lock(catalog_dir)
    try:
        log_name = os.path.join(catalog_dir, 'errors_log.tsv')
        cur_name = os.path.join(catalog_dir, 'errors_current.tsv')
        log_rows = read_table(log_name)
        current = dict(((r['participant_id'], r['session_id'], r['step']), r)
                       for r in read_table(cur_name))
        counts = [0, 0, 0]
        for occs, successes in found:
            for i, n in enumerate(merge(log_rows, current, occs, successes, now)):
                counts[i] += n
        merge_notes(current, read_table(os.path.join(catalog_dir,
                                                      'errors_notes.tsv')))
        write_table(log_name, log_rows, LOG_COLUMNS)
        write_table(cur_name, list(current.values()), CURRENT_COLUMNS)
        write_json(os.path.join(catalog_dir, 'errors_log.json'), LOG_COLUMNS)
        write_json(os.path.join(catalog_dir, 'errors_current.json'),
                   CURRENT_COLUMNS)
        stamp = copy_history(catalog_dir)
    finally:
        fid.close()
    if verbose:
        n_open = sum(1 for r in current.values() if _val(r, 'resolved') != '1')
        print('{0} sessions read: {1} new errors, {2} repeated, {3} resolved; '
              '{4} open in {5}'.format(len(found), counts[0], counts[1],
                                       counts[2], n_open, cur_name))
        if stamp:
            print('history copies made: {0}'.format(
                os.path.join(catalog_dir, 'history', '*_{0}.tsv'.format(stamp))))
    return 0, ''


#############################################################################
# command line
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Record processing errors by session and step in a '
                    'cumulative catalog.')
    parser.add_argument('--catalog_dir', required=True,
                        help='directory of errors_log.tsv, errors_current.tsv '
                             'and errors_notes.tsv')
    parser.add_argument('--rules', default=None,
                        help='rules file (default: '
                             '$CMIG_PIPELINE_DIR/parms/errors_rules.tsv)')
    parser.add_argument('--batchdir', nargs='+', default=None,
                        help='mkjobs batch directories to read')
    parser.add_argument('--participant_id', default=None,
                        help='one session, with --session_id')
    parser.add_argument('--session_id', default=None)
    parser.add_argument('--seslist_file', default=None,
                        help='tsv with participant_id and session_id columns')
    parser.add_argument('--script_name', default='cmig_pipeline_bids',
                        choices=sorted(mkjobs.SCRIPT_CONFIG),
                        help='script whose status files the listed sessions '
                             'have (default: cmig_pipeline_bids)')
    for key in ROOT_ARGS:
        parser.add_argument('--' + key, default=None)
    parser.add_argument('--func_flag', type=int, default=1, choices=[0, 1])
    parser.add_argument('--dwi_flag', type=int, default=1, choices=[0, 1])
    parser.add_argument('--task', nargs='+', default=None)
    parser.add_argument('--ended_flag', type=int, default=1, choices=[0, 1],
                        help='take the listed sessions\' jobs as ended, so a '
                             'missing status file is an error (default: 1)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1])
    return parser.parse_args()


def run(args):
    sessions = []
    for batchdir in args.batchdir or []:
        if not os.path.isdir(batchdir):
            print('ERROR: batchdir not found: {0}'.format(batchdir),
                  file=sys.stderr)
            return 1, 'batchdir not found: {0}'.format(batchdir)
        sessions += batch_sessions(batchdir)
    if args.seslist_file or args.participant_id:
        if not args.raw_rootdir:
            return 1, '--raw_rootdir is required with listed sessions'
        sessions += listed_sessions(args)
    if not sessions:
        return 1, 'no sessions given (--batchdir, --seslist_file or ' \
                  '--participant_id/--session_id)'
    return catalog_errors(args.catalog_dir, sessions, args.rules, args.verbose)


def main():
    errcode, errmsg = run(parse_args())
    if errcode:
        print('ERROR: {0}'.format(errmsg), file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
