#!/usr/bin/env python3
"""
Run a cmig pipeline (or module) over many BIDS sessions across several
storage devices and compute clusters, one "round" per call, as a cron job
would: settle the jobs submitted earlier, top each cluster up to its job
limit, and submit new jobs for the sessions not yet done.

A round:
   1. takes the run directory's lock (a round still running: exit 0)
   2. reads config.json and the tables (devices, clusters, sessions,
      assignments, jobs)
   3. asks each enabled cluster for the user's jobs; a cluster that does not
      answer is left alone this round
   4. settles the jobs that have left a queue, from their status files: done;
      failed (retried with errflag 1, then 2, then held; held at once where
      every failing step is a refusal, class not_processed in the error
      rules, since a rerun gives the same result); lost (status
      incomplete: resubmitted up to max_lost, then held); stale (a step the
      job left stale, logged as a WARNING: a candidate up to max_stale, then
      held); a job of this run found in a queue but not in jobs.tsv is
      adopted; the memory, CPU and disk I/O of the jobs that ended done or
      failed are read from the scheduler accounting (one query per cluster);
      each settled job's .out and .err are copied into every output root of
      the script, as logs/<pid>/<sid>/<pid>_<sid>_desc-<step><YYMMDDHHMM>.out
      and .err (step: the script name without cmig_ and _bids; the time: the
      .out's mtime), every run's kept (not with --dry_run 1)
   5. returns a done session to the candidates (state changed) where its
      rawdata scans.tsv or an events file changed after its last job was
      submitted or it was last found done; mkjobs then makes a job only for
      the steps that change made stale
   6. measures df of each device and its projected use
   7. cancels own pending jobs older than move_pending_hours on a cluster
      while another has room, returning their sessions to the pool
   8. works out a quota per cluster: job_limit + pending_buffer minus the
      user's jobs there
   9. takes candidate sessions (ready in sessions.tsv and, with ready_file,
      in that list; no job running, not done or held), assigns a device to each new participant, and writes a
      job script for each that has work left into the round's batchdir; the
      first time a round uses a device's roots, it creates them and writes
      their dataset_description.json files and the atlas where missing
      (mkjobs write_root_files, with the opts dataset_name and FS_version);
      with an image, the round first computes its sha256 (passed to the
      jobs for the status files) and reads its versions.json (for
      GeneratedBy)
  10. hands the jobs out one at a time to the cluster with the most quota
      left, and submits them (one ssh per cluster)
  11. writes the tables (each to .new, then renamed; changed ones copied to
      history/), a rounds.tsv row and the round's log

With --dry_run 1 nothing is submitted or cancelled and the tables are
written to <rundir>/dry_run/ instead.

Run directory:
  config.json     settings (see CONFIG_DEFAULTS; --init writes a template)
  devices.tsv     cmig_bids_devices.py (--add_device, after a device is cleared)
  clusters.tsv    cluster, scheduler, enabled, job_limit, pending_buffer,
                  clustmem, clustcpu, cstr
  sessions.tsv    cmig_list_bids_sessions.py
  assignments.tsv participant -> device, written as participants are assigned
  jobs.tsv        one row per session considered: device, state, attempts,
                  failures, lost, stale, cluster, jobid, jobname, batchdir, dates,
                  the last job's accounting (mem_gb, mem_type: rss on slurm,
                  vmem on SGE; cores, cpu_hours, read_gb: slurm only;
                  write_gb), errmsg, checked_date (found done without a
                  job); state queued, running, done, failed, lost, held,
                  missing, changed (done, then its rawdata changed), stale
                  (a job ended with steps stale)
  rounds.tsv      one row per round; per cluster (cluster:value), the jobs
                  that ended done or failed this round (finished), their
                  median mem_gb and cores, and their total write_gb and
                  read_gb (slurm only)

Created:  10/29/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/09/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/16/2026 by Don Hagler
"""

import argparse
import datetime
import fcntl
import getpass
import json
import math
import os
import re
import shutil
import socket
import statistics
import sys
import time

import cmig_bids_devices as devmod
import cmig_catalog_pipeline_errors as errcat
import cmig_cluster_utils as clu
import cmig_mkjobs_bids_pipeline as mkjobs
import cmig_select_bids_func_runs

TIME_FORMAT = '%Y-%m-%dT%H:%M:%S'
STAMP_FORMAT = '%Y%m%d_%H%M%S'

CONFIG_DEFAULTS = {
    'raw_rootdir': None,
    'seslist': ['ses-00A', 'ses-02A', 'ses-04A', 'ses-06A', 'ses-08A'],
    'script_name': 'cmig_pipeline_bids',
    # output roots under each device's rootdir, by the script's dir argument
    'device_roots': {'proc_rootdir': 'cmig_proc', 'fsurf_rootdir': 'freesurfer',
                     'analysis_rootdir': 'cmig_analysis',
                     'postqc_rootdir': 'cmig_postqc', 'rawqc_rootdir': 'cmig_rawqc'},
    # make_opts options for the job scripts (e.g. func_flag, cleanup_flag,
    #   appdir, appver, pydir; image, image_binds and fs_license for jobs
    #   that run the image through run.py); dataset_name (e.g. "ABCD Study")
    #   names the roots' dataset_description.json, and is passed to run.py
    'opts': {},
    'tracking_csv': None,
    'ready_only': 1,
    # a participant_id,session_id csv or tsv (e.g. the ready.csv of a rawdata
    #   readiness check): only its sessions are candidates; None: no list
    'ready_file': None,
    'avg_session_gb': 12.0,
    'transient_gb': 10.0,
    'max_frac': 0.70,
    'hard_frac': 0.90,
    'max_lost': 2,
    'max_stale': 2,
    'max_failures': 2,
    'max_submit_per_round': 500,
    # the pending buffer of a cluster whose pending_buffer is blank: the jobs
    #   it finishes between two rounds (job_limit x round_interval_hours /
    #   job_hours) times buffer_factor, so it still has jobs queued when the
    #   next round comes
    'round_interval_hours': 24.0,
    'job_hours': 21.0,
    'buffer_factor': 1.2,
    # a cluster with job_limit blank has no limit of ours (its scheduler
    #   shares it among users): keep queued what its running jobs finish
    #   between rounds, at least min_pending
    'min_pending': 20,
    'move_pending_hours': 0,
    'order': 'participant',
    'prefix': 'abcd_full',
    'batchroot': None,
    'logdir': None,
    'user': None,
    'history_keep': 30,
    # the error rules file (default: $CMIG_PIPELINE_DIR/parms/errors_rules.tsv)
    'error_rules': None,
}

JOB_COLUMNS = ['participant_id', 'session_id', 'device', 'state', 'attempts',
               'failures', 'lost', 'stale', 'cluster', 'jobid', 'jobname', 'batchdir',
               'submit_date', 'state_date', 'checked_date', 'mem_gb', 'mem_type',
               'cores', 'cpu_hours', 'read_gb', 'write_gb', 'errmsg']
# jobs.tsv columns from the accounting: (column, job_usage field, format)
USAGE_COLUMNS = [('mem_gb', 'mem_gb', '{0:.3f}'), ('mem_type', 'mem_type', '{0}'),
                 ('cores', 'cores', '{0:.3f}'), ('cpu_hours', 'cpu_h', '{0:.4f}'),
                 ('read_gb', 'read_gb', '{0:.3f}'),
                 ('write_gb', 'write_gb', '{0:.3f}')]
CLUSTER_COLUMNS = ['cluster', 'scheduler', 'enabled', 'job_limit',
                   'pending_buffer', 'clustmem', 'clustcpu', 'cstr']
ROUND_COLUMNS = ['date', 'duration_s', 'dry_run', 'queued', 'running', 'done',
                 'failed', 'lost', 'held', 'missing', 'changed', 'stale',
                 'settled', 'adopted', 'reopened', 'cancelled', 'submitted', 'assigned', 'unassigned',
                 'clusters', 'devices', 'finished', 'mem_gb_median',
                 'cores_median', 'write_gb_total', 'read_gb_total', 'errors']
IN_QUEUE = ('queued', 'running')
NOT_CANDIDATE = ('queued', 'running', 'done', 'held', 'missing')
JOBNAME_RE = re.compile(r'_job_\d+_(sub-[^_]+)_(ses-[^_]+)$')


def _now():
    return datetime.datetime.now().strftime(TIME_FORMAT)


def _parse_time(val):
    try:
        return datetime.datetime.strptime(val, TIME_FORMAT)
    except (TypeError, ValueError):
        return None


def _int(val, default=0):
    try:
        return int(val)
    except (TypeError, ValueError):
        return default


class Log(object):
    """Print to stdout and append to the round's log file."""

    def __init__(self, fname=None):
        self.fid = open(fname, 'a') if fname else None

    def __call__(self, msg):
        line = '{0} {1}'.format(_now(), msg)
        print(line)
        sys.stdout.flush()
        if self.fid:
            self.fid.write(line + '\n')
            self.fid.flush()

    def close(self):
        if self.fid:
            self.fid.close()


#############################################################################
# run directory
#############################################################################

def take_lock(rundir):
    """The open lock file, locked, or None if another round holds it."""
    fname = os.path.join(rundir, 'round.lock')
    fid = open(fname, 'a+')
    try:
        fcntl.flock(fid, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        fid.seek(0)
        holder = fid.read().strip()
        fid.close()
        print('round still running ({0}); exiting'.format(holder or 'unknown'))
        return None
    fid.seek(0)
    fid.truncate()
    fid.write('pid {0} host {1} start {2}\n'.format(
        os.getpid(), socket.gethostname(), _now()))
    fid.flush()
    return fid


def read_config(rundir):
    fname = os.path.join(rundir, 'config.json')
    with open(fname) as fid:
        config = dict(CONFIG_DEFAULTS, **json.load(fid))
    if not config['raw_rootdir']:
        raise ValueError('config.json: raw_rootdir not set')
    if config['script_name'] not in mkjobs.SCRIPT_CONFIG:
        raise ValueError('config.json: unknown script_name {0}'.format(
            config['script_name']))
    if mkjobs.SCRIPT_CONFIG[config['script_name']].get('image_only') and \
            not config['opts'].get('image'):
        raise ValueError('config.json: ' + mkjobs.IMAGE_ONLY_MSG.format(
            config['script_name']))
    if config['ready_file'] and not os.path.isfile(config['ready_file']):
        raise ValueError('config.json: ready_file not found: {0}'.format(
            config['ready_file']))
    config['user'] = config['user'] or getpass.getuser()
    config['batchroot'] = config['batchroot'] or os.path.join(rundir, 'batchdirs')
    config['logdir'] = config['logdir'] or os.path.join(rundir, 'logs')
    # read once a round; settle() classes the failures with them
    config['error_rules_loaded'] = errcat.load_rules(config['error_rules'])
    return config


def write_template(rundir):
    """config.json and clusters.tsv templates, where absent."""
    os.makedirs(rundir, exist_ok=True)
    fname = os.path.join(rundir, 'config.json')
    if not os.path.exists(fname):
        with open(fname, 'w') as fid:
            json.dump(CONFIG_DEFAULTS, fid, indent=2)
            fid.write('\n')
        print('wrote {0}'.format(fname))
    fname = os.path.join(rundir, 'clusters.tsv')
    if not os.path.exists(fname):
        # pending_buffer blank: computed from the round interval and job time
        # brightcluster1 has no per-user limit (its scheduler shares it among
        #   users): 250, the running jobs to expect there, to start
        rows = [{'cluster': 'brightcluster1', 'scheduler': 'slurm', 'enabled': 0,
                 'job_limit': 250, 'pending_buffer': '', 'clustmem': '60G',
                 'clustcpu': 8, 'cstr': 'qcshjob'}]
        for c in ['mmilcluster8', 'mmilcluster9', 'mmilcluster10']:
            rows.append({'cluster': c, 'scheduler': 'sge', 'enabled': 0,
                         'job_limit': 120, 'pending_buffer': '', 'clustmem': '60G',
                         'clustcpu': 8, 'cstr': 'qcshjob'})
        devmod.write_tsv(fname, rows, CLUSTER_COLUMNS)
        print('wrote {0}'.format(fname))


class Tables(object):
    """The run's tables, read at the start of a round and written at its end."""

    def __init__(self, rundir):
        self.rundir = rundir
        path = lambda name: os.path.join(rundir, name)
        self.devices = devmod.read_devices(path('devices.tsv'))
        self.assignments = devmod.read_assignments(path('assignments.tsv'))
        self.sessions = devmod.read_tsv(path('sessions.tsv'))
        self.clusters = [c for c in devmod.read_tsv(path('clusters.tsv'))
                         if _int(c.get('enabled'))]
        self.jobs = {(r['participant_id'], r['session_id']): r
                     for r in devmod.read_tsv(path('jobs.tsv'))}

    def write(self, outdir, history_keep, stamp):
        os.makedirs(outdir, exist_ok=True)
        path = lambda name: os.path.join(outdir, name)
        devmod.write_devices(path('devices.tsv'), self.devices)
        devmod.write_assignments(path('assignments.tsv'), self.assignments)
        devmod.write_tsv(path('jobs.tsv'),
                         [self.jobs[k] for k in sorted(self.jobs)], JOB_COLUMNS)
        histdir = os.path.join(outdir, 'history')
        os.makedirs(histdir, exist_ok=True)
        for name in ['devices', 'assignments', 'jobs']:
            shutil.copy2(path(name + '.tsv'),
                         os.path.join(histdir, '{0}_{1}.tsv'.format(name, stamp)))
            old = sorted(f for f in os.listdir(histdir) if f.startswith(name + '_'))
            for f in old[:-history_keep] if history_keep > 0 else []:
                os.remove(os.path.join(histdir, f))

    def job(self, pid, sid):
        """The jobs.tsv row of a session, created if new."""
        key = (pid, sid)
        if key not in self.jobs:
            self.jobs[key] = {c: '' for c in JOB_COLUMNS}
            self.jobs[key].update({'participant_id': pid, 'session_id': sid,
                                   'attempts': 0, 'failures': 0, 'lost': 0,
                                   'stale': 0})
        return self.jobs[key]


def set_state(row, state, errmsg=None):
    row['state'] = state
    row['state_date'] = _now()
    if errmsg is not None:
        row['errmsg'] = (errmsg or '').replace('\t', ' ').replace('\n', ' ')


def device_dirs(config, device):
    """The root dirs of the script for one device."""
    dirs = {'raw_rootdir': config['raw_rootdir']}
    for key, sub in config['device_roots'].items():
        dirs[key] = os.path.join(device['rootdir'], sub)
    return dirs


#############################################################################
# round steps
#############################################################################

def query_clusters(tables, config, log):
    """{cluster: {jobid: (jobname, state)}} for the clusters that answered."""
    answers = {}
    for c in tables.clusters:
        jobs = clu.query_jobs(c['cluster'], c['scheduler'], config['user'])
        if jobs is None:
            log('{0}: no answer; left alone this round'.format(c['cluster']))
            continue
        answers[c['cluster']] = jobs
        nrun, nqueue = clu.count_jobs(jobs)
        log('{0}: {1} running, {2} queued (all of {3}\'s jobs)'.format(
            c['cluster'], nrun, nqueue, config['user']))
    return answers


def not_processed(config, pid, sid, dirs, opts, rules):
    """True where the session has a failing status file and every one of
    them is a refusal (class not_processed): a rerun gives the same result
    until its inputs change."""
    classes = []
    for r in mkjobs.status_rows(config['script_name'], pid, sid, dirs, opts):
        st = r['status']
        if st is None or str(st.get('errcode', '')) == '0':
            continue
        module = errcat.step_fields(r['step'])['module']
        classes.append(errcat.classify(st.get('errmsg', ''), module, rules)[0])
    return bool(classes) and all(c == 'not_processed' for c in classes)


def log_step(script_name):
    """The desc- label of a script's job logs: cmig_pipeline_bids_rawqc ->
    pipelinerawqc."""
    name = re.sub(r'^cmig_', '', script_name).replace('_bids', '')
    return re.sub(r'[^A-Za-z0-9]', '', name)


def copy_job_logs(row, script_name, dirs, log):
    """Copy a settled job's .out and .err into each output root's logs/.
    Done here rather than in the job script, since a job's .out is still open
    at its own end and a job the scheduler kills never reaches a last line."""
    batchdir = row.get('batchdir') or ''
    prefix = os.path.basename(os.path.normpath(batchdir)) + '_'
    jobname = row.get('jobname') or ''
    if not batchdir or not jobname.startswith(prefix):
        return
    # the paths clu.submit_lines gives the scheduler
    src = os.path.join(batchdir, 'pbsout', jobname[len(prefix):])
    if not os.path.isfile(src + '.out'):
        return
    pid, sid = row['participant_id'], row['session_id']
    stamp = datetime.datetime.fromtimestamp(
        os.path.getmtime(src + '.out')).strftime('%y%m%d%H%M')
    stem = '{0}_{1}_desc-{2}{3}'.format(pid, sid, log_step(script_name), stamp)
    for key in mkjobs._output_roots(mkjobs.SCRIPT_CONFIG[script_name]):
        root = dirs.get(key)
        if not root or not os.path.isdir(root):
            continue
        logdir = os.path.join(root, 'logs', pid, sid)
        try:
            os.makedirs(logdir, exist_ok=True)
            for ext in ('.out', '.err'):
                dst = os.path.join(logdir, stem + ext)
                if os.path.isfile(src + ext) and not os.path.exists(dst):
                    shutil.copy2(src + ext, dst)
        except OSError as exc:
            log('WARNING: {0} {1}: job logs not copied to {2}: {3}'.format(
                pid, sid, logdir, exc))


def settle(tables, config, answers, opts, delay_flags, devices_by_name, log,
           dry_run=False):
    """Update the jobs in a queue; settle those that left one. Returns
    (n_settled, n_adopted, the rows that ended done or failed)."""
    n_settled = 0
    finished = []
    in_queue = {}
    for cluster, jobs in answers.items():
        for jobid, (name, state) in jobs.items():
            in_queue[(cluster, jobid)] = (name, state)
    for key, row in sorted(tables.jobs.items()):
        if row.get('state') not in IN_QUEUE or row.get('cluster') not in answers:
            continue
        found = in_queue.get((row['cluster'], row['jobid']))
        if found:
            state = found[1]
            if state in IN_QUEUE and state != row['state']:
                set_state(row, state)
            continue
        # gone from the queue: the status files say how it ended
        device = devices_by_name.get(row['device'])
        if device is None:
            continue
        dirs = device_dirs(config, device)
        if not dry_run:
            copy_job_logs(row, config['script_name'], dirs, log)
        status, errmsg = mkjobs.check_session(
            config['script_name'], row['participant_id'], row['session_id'],
            dirs, opts, delay_flags=delay_flags)
        n_settled += 1
        if status in ('done', 'failed'):
            finished.append(row)
        if status == 'done':
            set_state(row, 'done', '')
        elif status == 'failed':
            row['failures'] = _int(row.get('failures')) + 1
            if row['failures'] > config['max_failures']:
                set_state(row, 'held', errmsg)
            elif not_processed(config, row['participant_id'],
                               row['session_id'], dirs, opts,
                               config['error_rules_loaded']):
                # no errflag 1/2 retry, which would delete and recompute
                #   outputs whose outcome the inputs fix
                set_state(row, 'held', errmsg)
            else:
                set_state(row, 'failed', errmsg)
        elif status == 'stale':
            # the job left steps it did not run stale (or the rawdata changed
            #   while it ran): a candidate, whose next job reruns them
            row['stale'] = _int(row.get('stale')) + 1
            if row['stale'] > config['max_stale']:
                set_state(row, 'held', 'ended stale {0} times: {1}'.format(
                    row['stale'], errmsg))
            else:
                set_state(row, 'stale', errmsg)
            log('WARNING: {0} {1}: job ended stale: {2}'.format(
                row['participant_id'], row['session_id'], errmsg))
        else:
            row['lost'] = _int(row.get('lost')) + 1
            if row['lost'] > config['max_lost']:
                set_state(row, 'held', 'lost {0} times'.format(row['lost']))
            else:
                set_state(row, 'lost', 'job ended with status {0}'.format(status))
        log('{0} {1}: {2}{3}'.format(row['participant_id'], row['session_id'],
                                     row['state'],
                                     ': ' + row['errmsg'] if row.get('errmsg') else ''))
    # a job of this run in a queue but not in jobs.tsv (a round that crashed
    #   after submitting, or a job submitted by hand) is adopted
    n_adopted = 0
    known = set((r.get('cluster'), r.get('jobid')) for r in tables.jobs.values())
    prefix = config['prefix'] + '_'
    for (cluster, jobid), (name, state) in sorted(in_queue.items()):
        if (cluster, jobid) in known or not name.startswith(prefix):
            continue
        m = JOBNAME_RE.search(name)
        if not m or state not in IN_QUEUE:
            continue
        row = tables.job(m.group(1), m.group(2))
        if row.get('state') in IN_QUEUE:
            continue
        row.update({'cluster': cluster, 'jobid': jobid, 'jobname': name})
        set_state(row, state)
        if not row.get('device'):
            a = tables.assignments.get(m.group(1))
            row['device'] = a['device'] if a else ''
        n_adopted += 1
        log('adopted {0} on {1} ({2})'.format(name, cluster, jobid))
    return n_settled, n_adopted, finished


def record_usage(tables, finished, log):
    """Write each finished job's accounting into its jobs.tsv row, one
    query per cluster; return the per-cluster summary columns of
    rounds.tsv."""
    sched = {c['cluster']: c['scheduler'] for c in tables.clusters}
    by_cluster = {}
    for row in finished:
        by_cluster.setdefault(row['cluster'], []).append(row)
    summary = {k: [] for k in ['finished', 'mem_gb_median', 'cores_median',
                               'write_gb_total', 'read_gb_total']}
    for cluster in sorted(by_cluster):
        rows = by_cluster[cluster]
        usage = None
        if cluster in sched:
            try:
                usage = clu.job_usage(cluster, sched[cluster],
                                      [r['jobid'] for r in rows])
            except ValueError as e:
                log('{0}: accounting not queried: {1}'.format(cluster, e))
        if usage is None:
            log('{0}: no accounting for {1} finished jobs'.format(cluster, len(rows)))
            usage = {}
        found = [usage[r['jobid']] for r in rows if r['jobid'] in usage]
        for row in rows:
            u = usage.get(row['jobid'], {})
            for col, key, fmt in USAGE_COLUMNS:
                row[col] = '' if u.get(key) is None else fmt.format(u[key])

        def values(key):
            return [u[key] for u in found if u.get(key) is not None]

        def median(key):
            v = values(key)
            return '{0:.2f}'.format(statistics.median(v)) if v else 'na'
        summary['finished'].append('{0}:{1}'.format(cluster, len(rows)))
        summary['mem_gb_median'].append('{0}:{1}'.format(cluster, median('mem_gb')))
        summary['cores_median'].append('{0}:{1}'.format(cluster, median('cores')))
        summary['write_gb_total'].append('{0}:{1:.1f}'.format(
            cluster, sum(values('write_gb'))))
        # read is counted on slurm only
        if sched.get(cluster) == 'slurm':
            summary['read_gb_total'].append('{0}:{1:.1f}'.format(
                cluster, sum(values('read_gb'))))
        log('{0}: accounting for {1} of {2} finished jobs'.format(
            cluster, len(found), len(rows)))
    return {k: ' '.join(v) for k, v in summary.items()}


def move_pending(tables, config, answers, quotas, dry_run, log):
    """Cancel own pending jobs older than move_pending_hours on clusters
    while another cluster has room and nothing pending; their sessions
    return to the pool. Returns the number cancelled."""
    hours = config['move_pending_hours']
    if not hours:
        return 0
    room = sum(q for c, q in quotas.items()
               if not any(r.get('state') == 'queued' and r.get('cluster') == c
                          for r in tables.jobs.values()))
    if room <= 0:
        return 0
    now = datetime.datetime.now()
    old = []
    for row in tables.jobs.values():
        if row.get('state') != 'queued' or row.get('cluster') not in answers:
            continue
        t = _parse_time(row.get('submit_date'))
        if t and (now - t).total_seconds() > hours * 3600:
            old.append(row)
    # the newest first: the oldest keep their place in their queue
    old.sort(key=lambda r: r.get('submit_date'), reverse=True)
    old = old[:room]
    n = 0
    sched = {c['cluster']: c['scheduler'] for c in tables.clusters}
    for cluster in sorted(set(r['cluster'] for r in old)):
        rows = [r for r in old if r['cluster'] == cluster]
        ok = clu.cancel_jobs(cluster, sched[cluster], [r['jobid'] for r in rows],
                             dry_run=dry_run)
        if ok or dry_run:
            for r in rows:
                set_state(r, 'lost', 'pending over {0} h on {1}: moved'.format(
                    hours, cluster))
                n += 1
    if n:
        log('cancelled {0} pending jobs to move them'.format(n))
    return n


def pending_buffer(cluster, config):
    """The jobs to keep queued beyond job_limit: the cluster's
    pending_buffer, or, where blank, the jobs it finishes between two rounds
    times buffer_factor."""
    if str(cluster.get('pending_buffer', '')).strip() != '':
        return _int(cluster['pending_buffer'])
    per_round = (_int(cluster.get('job_limit')) * config['round_interval_hours'] /
                 max(config['job_hours'], 0.1))
    return int(math.ceil(per_round * config['buffer_factor']))


def cluster_quotas(tables, answers, config):
    """{cluster: jobs it can take this round}, so that as running jobs
    finish, queued ones take their place until the next round:
      job_limit set: job_limit + its pending buffer minus all of the user's
        jobs there (running and queued)
      job_limit blank (no limit of ours; the scheduler shares the cluster):
        the jobs the running ones finish between rounds, times
        buffer_factor and at least min_pending (and, with none left queued,
        at least as many as run), minus those queued"""
    quotas = {}
    for c in tables.clusters:
        if c['cluster'] not in answers:
            continue
        nrun, nqueue = clu.count_jobs(answers[c['cluster']])
        if str(c.get('job_limit', '')).strip() == '':
            per_round = nrun * config['round_interval_hours'] / max(config['job_hours'], 0.1)
            target = max(config['min_pending'],
                         int(math.ceil(per_round * config['buffer_factor'])))
            # nothing left queued: the cluster took everything; queue at
            #   least as many again as run, so the count grows until a queue
            #   builds
            if nqueue == 0:
                target = max(target, nrun)
            quotas[c['cluster']] = max(0, target - nqueue)
        else:
            quotas[c['cluster']] = max(0, _int(c.get('job_limit')) +
                                       pending_buffer(c, config) - nrun - nqueue)
    return quotas


def _rawdata_mtime(raw_rootdir, pid, sid):
    """The newest mtime of a session's rawdata scans.tsv and events files, or
    None."""
    times = []
    fscans = os.path.join(raw_rootdir, pid, sid, '{0}_{1}_scans.tsv'.format(pid, sid))
    if os.path.isfile(fscans):
        times.append(datetime.datetime.fromtimestamp(os.path.getmtime(fscans)))
    events = cmig_select_bids_func_runs.events_mtime(raw_rootdir, pid, sid)
    if events is not None:
        times.append(events)
    return max(times) if times else None


def reopen_changed(tables, config, log):
    """Return to the candidates (state changed) each done session whose
    rawdata scans.tsv or events changed after its last job's submit_date, or
    after it was last found done without a job (checked_date), so mkjobs
    compares it again; returns their number."""
    n = 0
    for (pid, sid), row in sorted(tables.jobs.items()):
        if row.get('state') != 'done':
            continue
        refs = [t for t in (_parse_time(row.get('submit_date')),
                            _parse_time(row.get('checked_date'))) if t is not None]
        if not refs:
            continue
        newest = _rawdata_mtime(config['raw_rootdir'], pid, sid)
        # the dates are kept to the second
        if newest is not None and newest.replace(microsecond=0) > max(refs):
            set_state(row, 'changed', 'rawdata scans.tsv or events changed {0}'.format(
                newest.strftime(TIME_FORMAT)))
            n += 1
    if n:
        log('{0} done sessions returned to the candidates: rawdata changed'.format(n))
    return n


def candidates(tables, config, eligible):
    """The (participant, session) pairs to consider, in the config order."""
    pairs = [(pid, sid) for pid in eligible for sid in eligible[pid]]
    if config['order'] == 'session':
        pairs.sort(key=lambda p: (p[1], p[0]))
    else:
        pairs.sort()
    return [p for p in pairs
            if tables.jobs.get(p, {}).get('state') not in NOT_CANDIDATE]


def make_jobs(tables, config, eligible, quota_total, opts, delay_flags,
              batchdir, df_failed, running_by_device, log):
    """Write job scripts for up to quota_total candidates; return their
    [(job_name, row)] and the numbers assigned and left unassigned."""
    script_name = config['script_name']
    devices_by_name = {d['device']: d for d in tables.devices}
    checked = set()
    made = []
    n_assigned = n_unassigned = 0
    no_device = set()
    pairs = candidates(tables, config, eligible)
    # the image's sha256 and versions, once per round with jobs to make
    if pairs and quota_total > 0:
        mkjobs.add_image_info(opts, log)
    for pid, sid in pairs:
        if len(made) >= quota_total:
            break
        if pid in no_device:
            continue
        if pid not in tables.assignments:
            name = devmod.assign(
                tables.devices, tables.assignments, pid, len(eligible[pid]),
                config['avg_session_gb'], config['max_frac'], config['hard_frac'],
                running=running_by_device, transient_gb=config['transient_gb'])
            if name is None:
                no_device.add(pid)
                n_unassigned += 1
                continue
            n_assigned += 1
        device = devices_by_name[tables.assignments[pid]['device']]
        if device['device'] in df_failed:
            continue
        dirs = device_dirs(config, device)
        if device['device'] not in checked:
            mkjobs.check_dirs(script_name, dirs)
            # kept where they exist, so only the first round writes them
            for msg in mkjobs.write_root_files(
                    script_name, dirs, dataset_name=opts.get('dataset_name'),
                    FS_version=opts['FS_version'],
                    versions=opts.get('image_versions')):
                log('{0}: {1}'.format(device['device'], msg))
            checked.add(device['device'])
        row = tables.job(pid, sid)
        row['device'] = device['device']
        # a retry after a failure deletes outputs: errflag 1, then 2
        errflag = min(_int(row.get('failures')), 2) if row.get('state') == 'failed' else 0
        job_name = 'job_{0:04d}_{1}_{2}'.format(len(made) + 1, pid, sid)
        status, errmsg, fpath = mkjobs.make_session_job(
            script_name, pid, sid, dirs, opts, batchdir, job_name,
            errflag=errflag, newflag=True, delay_flags=delay_flags, verbose=False)
        if fpath is None:
            if status == 'done':
                set_state(row, 'done', '')
                row['checked_date'] = _now()
            elif status == 'missing':
                set_state(row, 'missing', errmsg)
            elif status == 'failed':
                # failed before the rounds saw it: retried with errflag 1 next
                row['failures'] = _int(row.get('failures')) + 1
                set_state(row, 'failed', errmsg)
            continue
        made.append((job_name, row))
    if no_device:
        log('{0} participants found no device (max_frac {1}); their sessions '
            'wait'.format(len(no_device), config['max_frac']))
    return made, n_assigned, n_unassigned


def distribute(made, quotas):
    """{cluster: [(job_name, row)]}: one job at a time to the cluster with
    the most quota left."""
    left = dict(quotas)
    out = {c: [] for c in quotas}
    for item in made:
        c = max(sorted(left), key=lambda k: left[k]) if left else None
        if c is None or left[c] <= 0:
            break
        out[c].append(item)
        left[c] -= 1
    return out


def submit(tables, batchdir, by_cluster, dry_run, log):
    """Submit each cluster's jobs; record those with a job id as queued."""
    n = 0
    clusters = {c['cluster']: c for c in tables.clusters}
    batchname = os.path.basename(batchdir)
    for cluster, items in sorted(by_cluster.items()):
        if not items:
            continue
        c = clusters[cluster]
        names = [j for j, _ in items]
        with open(os.path.join(batchdir, 'scriptlist_{0}.txt'.format(cluster)), 'w') as fid:
            fid.write(''.join(j + '\n' for j in names))
        ids = clu.submit_jobs(cluster, c['scheduler'], batchdir, names,
                              c.get('clustmem') or '60G', _int(c.get('clustcpu'), 1),
                              cstr=c.get('cstr') or 'qcshjob', batchname=batchname,
                              dry_run=dry_run)
        for job_name, row in items:
            jobid = ids.get(job_name)
            if jobid is None and not dry_run:
                log('{0}: no job id from {1}; back to the pool'.format(job_name, cluster))
                continue
            row.update({'cluster': cluster, 'jobid': jobid or 'dry_run',
                        'jobname': '{0}_{1}'.format(batchname, job_name),
                        'batchdir': batchdir, 'submit_date': _now()})
            # the previous job's accounting
            row.update({col: '' for col, _, _ in USAGE_COLUMNS})
            row['attempts'] = _int(row.get('attempts')) + 1
            set_state(row, 'queued')
            n += 1
        log('{0}: {1} {2}'.format(cluster, len(items),
                                   'jobs (dry run)' if dry_run else 'jobs submitted'))
    return n


def round_row(tables, stamp, t0, dry_run, counts, answers, quotas, errors):
    states = {}
    for row in tables.jobs.values():
        states[row.get('state')] = states.get(row.get('state'), 0) + 1
    row = {'date': stamp, 'duration_s': int(time.time() - t0), 'dry_run': int(dry_run)}
    for s in ['queued', 'running', 'done', 'failed', 'lost', 'held', 'missing',
              'changed', 'stale']:
        row[s] = states.get(s, 0)
    row.update(counts)
    row['clusters'] = ' '.join(
        '{0}:{1}'.format(c['cluster'], 'q{0}'.format(quotas[c['cluster']])
                         if c['cluster'] in answers else 'unreachable')
        for c in tables.clusters)
    row['devices'] = ' '.join(
        '{0}:{1:.3f}/{2}'.format(d['device'], d.get('projected_frac') or 0,
                                 '{0:.3f}'.format(d['df_frac'])
                                 if d.get('df_frac') is not None else 'na')
        for d in tables.devices)
    row['errors'] = '; '.join(errors)
    return row


def run_round(rundir, dry_run=False, max_submit=None):
    config = read_config(rundir)
    stamp = datetime.datetime.now().strftime(STAMP_FORMAT)
    os.makedirs(config['logdir'], exist_ok=True)
    log = Log(os.path.join(config['logdir'], '{0}_round_{1}{2}.log'.format(
        config['prefix'], stamp, '_dry_run' if dry_run else '')))
    t0 = time.time()
    errors = []
    log('round {0} in {1}{2}'.format(stamp, rundir, ' (dry run)' if dry_run else ''))

    tables = Tables(rundir)
    script_name = config['script_name']
    opts = mkjobs.make_opts(**config['opts'])
    delay_flags = mkjobs.load_delay_flags(script_name, config['tracking_csv'])
    # the script's required input dirs that sessions.tsv records
    required = [m for m in mkjobs.SCRIPT_CONFIG[script_name].get('required_subdirs', [])
                if m in ('anat', 'func', 'dwi', 'fmap', 'beh')]
    ready_pairs = None
    if config['ready_file']:
        ready_pairs = devmod.read_ready_file(config['ready_file'])
        log('{0} sessions in the ready file {1}'.format(
            len(ready_pairs), config['ready_file']))
    eligible = devmod.eligible_sessions(tables.sessions, config['seslist'],
                                        required, bool(config['ready_only']),
                                        ready_pairs)
    log('{0} eligible sessions of {1} participants'.format(
        sum(len(v) for v in eligible.values()), len(eligible)))

    answers = query_clusters(tables, config, log)
    devices_by_name = {d['device']: d for d in tables.devices}
    n_settled, n_adopted, finished = settle(tables, config, answers, opts,
                                            delay_flags, devices_by_name, log,
                                            dry_run=dry_run)
    usage_summary = record_usage(tables, finished, log)
    n_reopened = reopen_changed(tables, config, log)

    df_failed = devmod.update_df(tables.devices)
    for name in df_failed:
        errors.append('df failed: {0}'.format(name))
        log('{0}: df failed; no jobs for its participants this round'.format(name))
    devmod.update_projected(tables.devices, config['avg_session_gb'])
    running_by_device = {}
    for row in tables.jobs.values():
        if row.get('state') == 'running':
            running_by_device[row['device']] = running_by_device.get(row['device'], 0) + 1

    quotas = cluster_quotas(tables, answers, config)
    n_cancelled = move_pending(tables, config, answers, quotas, dry_run, log)
    quota_total = sum(quotas.values())
    limit = config['max_submit_per_round'] if max_submit is None else max_submit
    quota_total = min(quota_total, limit)
    log('quota: {0} ({1})'.format(quota_total, ' '.join(
        '{0}:{1}'.format(c, q) for c, q in sorted(quotas.items()))))

    batchdir = os.path.join(config['batchroot'], '{0}_{1}{2}'.format(
        config['prefix'], stamp, '_dry_run' if dry_run else ''))
    os.makedirs(batchdir, exist_ok=True)
    made, n_assigned, n_unassigned = make_jobs(
        tables, config, eligible, quota_total, opts, delay_flags, batchdir,
        df_failed, running_by_device, log)
    n_submitted = submit(tables, batchdir, distribute(made, quotas), dry_run, log)

    counts = {'settled': n_settled, 'adopted': n_adopted, 'reopened': n_reopened,
              'cancelled': n_cancelled,
              'submitted': n_submitted, 'assigned': n_assigned,
              'unassigned': n_unassigned}
    counts.update(usage_summary)
    outdir = os.path.join(rundir, 'dry_run') if dry_run else rundir
    tables.write(outdir, config['history_keep'], stamp)
    row = round_row(tables, stamp, t0, dry_run, counts, answers, quotas, errors)
    fname = os.path.join(outdir, 'rounds.tsv')
    rows = devmod.read_tsv(fname) + [row]
    devmod.write_tsv(fname, rows, ROUND_COLUMNS)
    log('round done: {0}'.format(' '.join('{0} {1}'.format(k, row[k]) for k in
                                          ROUND_COLUMNS[3:19])))
    log.close()
    return 0, None


#############################################################################
# command line
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='One round of the cmig pipeline over many sessions: settle '
                    'the jobs, top each cluster up, submit new jobs.')
    parser.add_argument('--rundir', required=True,
                        help='run directory (config.json and the tables)')
    parser.add_argument('--dry_run', type=int, default=0, choices=[0, 1],
                        help='write the job scripts but submit and cancel nothing; '
                             'tables to <rundir>/dry_run (default: 0)')
    parser.add_argument('--max_submit', type=int, default=None,
                        help='at most this many jobs this round '
                             '(default: config max_submit_per_round)')
    parser.add_argument('--init', type=int, default=0, choices=[0, 1],
                        help='write config.json and clusters.tsv templates '
                             'where absent, and exit (default: 0)')
    return parser.parse_args()


def run(args):
    if args.init:
        write_template(args.rundir)
        return 0, None
    if not os.path.isdir(args.rundir):
        errmsg = 'rundir not found: {0}'.format(args.rundir)
        print('ERROR: {0}'.format(errmsg), file=sys.stderr)
        return 1, errmsg
    lock = take_lock(args.rundir)
    if lock is None:
        return 0, None
    try:
        return run_round(args.rundir, bool(args.dry_run), args.max_submit)
    except ValueError as e:
        print('ERROR: {0}'.format(e), file=sys.stderr)
        return 1, str(e)
    finally:
        lock.close()


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
