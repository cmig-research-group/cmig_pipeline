#!/usr/bin/env python3
"""
Run a cmig pipeline (or module) over many BIDS sessions across several
storage devices and compute clusters, one "round" per call, as a cron job
would: settle the jobs submitted earlier, top each cluster up to its job
limit, and submit new jobs for the sessions not yet done.

A run directory runs one job type (the config's script_name, prefix, opts
and device_roots), or several (config job_types, e.g. raw QC and the full
pipeline): each type has its own jobs (jobs.tsv's job_type column), job
names (its prefix) and batchdir per round; the types share the devices and
assignments (a participant is assigned once, by the first type that takes
it), the ready list, the clusters' quotas (counted over all of the user's
jobs) and max_submit_per_round. The types are made and submitted in the
order listed, so list the short ones first to give them the free slots.

A round:
   1. takes the run directory's lock (a round still running: exit 0)
   2. reads config.json and the tables (devices, clusters, sessions,
      assignments, jobs)
   3. asks each enabled cluster for the user's jobs; a cluster that does not
      answer is left alone this round
   4. settles the jobs that have left a queue, from their status files: done;
      failed (retried with errflag 1, then 2, then held; held at once where
      every failing step is a refusal, class not_processed in the error
      rules, since a rerun gives the same result; also a job whose status
      is incomplete but whose .out shows it ran to its end with an error,
      log_crash); lost (status incomplete otherwise: resubmitted up to
      max_lost, then held); stale (a step the
      job left stale, logged as a WARNING: a candidate up to max_stale, then
      held); a job of this run found in a queue but not in jobs.tsv is
      adopted; the memory, CPU and disk I/O of the jobs that ended done,
      failed or lost are read from the scheduler accounting (one query per cluster);
      each settled job's .out and .err are copied into every output root of
      the script, as logs/<pid>/<sid>/<pid>_<sid>_desc-<step><YYMMDDHHMM>.out
      and .err (step: the script name without cmig_ and _bids; the time: the
      .out's mtime), every run's kept (not with --dry_run 1); a failed or
      lost job that ran out of memory (slurm state OUT_OF_MEMORY, SGE
      exit_status 137, "Out of memory" in its .out or .err, or a failed
      status errmsg matching the error rule rs_memory) counts in
      jobs.tsv's oom, and its retries get clustmem x retry_mem_factor^oom
      (oom up to max_failures)
   5. returns a done session to the candidates (state changed) where its
      rawdata scans.tsv or an events file changed after its last job was
      submitted or it was last found done; mkjobs then makes a job only for
      the steps that change made stale
   6. measures df of each device and its projected use
   7. cancels own pending jobs older than move_pending_hours on a cluster
      while another has room, returning their sessions to the pool
   8. works out a quota per cluster: job_limit + pending_buffer minus the
      user's jobs there
   9. for each job type in turn: takes candidate sessions (ready in
      sessions.tsv and, with ready_file, in that list; no job of the type
      running, not done or held), assigns a device to each new participant,
      and writes a job script for each that has work left into the type's
      batchdir of the round; the first time a round uses a device's roots,
      it creates them and writes their dataset_description.json files and
      the atlas where missing
      (mkjobs write_root_files, with the opts dataset_name and FS_version);
      with an image, the round first computes its sha256 (passed to the
      jobs for the status files) and reads its versions.json (for
      GeneratedBy), once for the types sharing an image
  10. hands the type's jobs out one at a time to the cluster (of the type's)
      with the most quota left, and submits them (one ssh per cluster and
      type)
  11. writes the tables (each to .new, then renamed; changed ones copied to
      history/), a rounds.tsv row and the round's log

With --dry_run 1 nothing is submitted or cancelled and the tables are
written to <rundir>/dry_run/ instead.

Run directory:
  config.json     settings (see CONFIG_DEFAULTS; --init writes a template)
  devices.tsv     cmig_bids_devices.py (--add_device, after a device is cleared)
  clusters.tsv    cluster, scheduler, enabled, job_limit, pending_buffer,
                  clustmem, clustcpu, cstr (a job type's clustmem and
                  clustcpu, where given, take precedence)
  sessions.tsv    cmig_list_bids_sessions.py
  assignments.tsv participant -> device, written as participants are assigned
  jobs.tsv        one row per session and job type considered: job_type
                  ('' for the one type of a config without job_types; a
                  row without one is read as the type whose prefix starts
                  its jobname, else the type of the top-level script_name,
                  else the first; a row of a type not in the config is kept
                  as it is), device, state, attempts, failures, lost,
                  stale, oom (the job ended out of memory), cluster, jobid,
                  jobname, batchdir, clustmem (the last job's), dates,
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

Job types (config job_types; without it, one type from the top level):
  "job_types": [
    {"name": "rawqc", "script_name": "cmig_pipeline_bids_rawqc",
     "prefix": "abcd_rawqc", "opts": {"pc_flag": 1, ...},
     "device_roots": {"rawqc_rootdir": "cmig_rawqc", ...},
     "clusters": ["mmilcluster8"], "clustmem": "20G", "clustcpu": 1},
    {"name": "full", "script_name": "cmig_pipeline_bids",
     "prefix": "abcd_full", "opts": {...},
     "clusters": ["brightcluster1", "mmilcluster9", "mmilcluster10"],
     "clustmem": "28G", "clustcpu": {"brightcluster1": 8, "mmilcluster9": 2}}]
  name is required; script_name, prefix, opts and device_roots left out
  are the top level's (each replaces it whole, opts are not merged); the
  prefixes (job names) must differ, neither prefix + '_' starting the
  other's; clusters (default: all enabled); clustmem and clustcpu, one
  value or {cluster: value}, else the cluster's in clusters.tsv. The
  top-level prefix names the round's log.

Created:  10/29/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/09/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Prev mod: 11/16/2026 by Don Hagler
Prev mod: 11/18/2026 by Don Hagler
Prev mod: 11/20/2026 by Don Hagler
Prev mod: 11/21/2026 by Don Hagler
Last mod: 11/23/2026 by Don Hagler
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
    # a retry after a job ran out of memory (slurm OUT_OF_MEMORY, SGE
    #   exit_status 137, "Out of memory" in its log or an out of memory status
    #   errmsg, error rule rs_memory) gets clustmem times this
    #   factor per such ending, up to max_failures of them (1.5: 1.5x, 2.25x)
    'retry_mem_factor': 1.5,
    # several job types from this run directory: a list of {"name",
    #   "script_name", "prefix", "opts", "device_roots", "clusters",
    #   "clustmem", "clustcpu"} (see the docstring), made and submitted in
    #   that order; None: one type, from script_name, prefix, opts and
    #   device_roots above
    'job_types': None,
}
# the keys of a job type; those also at the top level default to its value
TYPE_KEYS = ['name', 'script_name', 'prefix', 'opts', 'device_roots',
             'clusters', 'clustmem', 'clustcpu']
TYPE_TOP_KEYS = ['script_name', 'prefix', 'opts', 'device_roots']

JOB_COLUMNS = ['participant_id', 'session_id', 'job_type', 'device', 'state',
               'attempts', 'failures', 'lost', 'stale', 'oom', 'cluster', 'jobid',
               'jobname', 'batchdir', 'clustmem', 'submit_date', 'state_date', 'checked_date', 'mem_gb', 'mem_type',
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
# MATLAB's out of memory error, in a job's .out or .err
OOM_LOG_RE = re.compile(r'\bOut of memory\b')
# the closing line of a job that ran to its end through run.py
RUNPY_END_RE = re.compile(r'^run\.py: (\d+) of \d+ \S+( \S+)? failed\s*$', re.MULTILINE)
ERROR_LINE_RE = re.compile(r': ERROR: ')
# the error rules' resource class, and the rule of it whose match in a failed
#   status errmsg gets the memory raise
RESOURCE_CLASS = 'resource'
OOM_RULES = ('rs_memory',)


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
    config['types'] = job_types(config)
    if config['ready_file'] and not os.path.isfile(config['ready_file']):
        raise ValueError('config.json: ready_file not found: {0}'.format(
            config['ready_file']))
    config['user'] = config['user'] or getpass.getuser()
    config['batchroot'] = config['batchroot'] or os.path.join(rundir, 'batchdirs')
    config['logdir'] = config['logdir'] or os.path.join(rundir, 'logs')
    # read once a round; settle() classes the failures with them
    config['error_rules_loaded'] = errcat.load_rules(config['error_rules'])
    # a failed status errmsg matching an out of memory rule gets the memory raise
    config['oom_re'] = errcat.class_regex(
        config['error_rules_loaded'], RESOURCE_CLASS, OOM_RULES)
    return config


def job_types(config):
    """The run's job types, each a dict of TYPE_KEYS: config job_types, the
    keys a type leaves out among TYPE_TOP_KEYS taken from the top level; or,
    without job_types, one type named '' from the top level. Raises
    ValueError."""
    given = config.get('job_types')
    if given is not None and not isinstance(given, list):
        raise ValueError('config.json: job_types must be a list')
    multi = bool(given)
    if not multi:
        given = [{'name': ''}]
    types = []
    for t in given:
        if not isinstance(t, dict):
            raise ValueError('config.json: job_types: each must be an object')
        unknown = sorted(set(t) - set(TYPE_KEYS))
        if unknown:
            raise ValueError('config.json: job_types: unknown keys {0}'.format(
                ' '.join(unknown)))
        jt = {k: t.get(k) for k in TYPE_KEYS}
        for k in TYPE_TOP_KEYS:
            if k not in t:
                jt[k] = config.get(k)
        jt['name'] = str(jt['name'] or '')
        jt['opts'] = jt['opts'] or {}
        where = 'job type {0}: '.format(jt['name']) if multi else ''
        if multi and not jt['name']:
            raise ValueError('config.json: job_types: each needs a name')
        if not jt['prefix']:
            raise ValueError('config.json: {0}prefix not set'.format(where))
        if jt['script_name'] not in mkjobs.SCRIPT_CONFIG:
            raise ValueError('config.json: {0}unknown script_name {1}'.format(
                where, jt['script_name']))
        if mkjobs.SCRIPT_CONFIG[jt['script_name']].get('image_only') and \
                not jt['opts'].get('image'):
            raise ValueError('config.json: ' + where + mkjobs.IMAGE_ONLY_MSG.format(
                jt['script_name']))
        if jt['clusters'] is not None and not isinstance(jt['clusters'], list):
            raise ValueError('config.json: {0}clusters must be a list'.format(where))
        types.append(jt)
    names = [jt['name'] for jt in types]
    if len(set(names)) != len(names):
        raise ValueError('config.json: job_types: names not unique')
    # adoption tells the types' jobs apart by their prefix
    for a in types:
        for b in types:
            if a is not b and (a['prefix'] + '_').startswith(b['prefix'] + '_'):
                raise ValueError('config.json: job_types: prefix {0} starts with '
                                 'prefix {1}'.format(a['prefix'], b['prefix']))
    return types


def _tag(jt):
    """The log prefix of a job type: '<name>: ', or '' for the one type of a
    config without job_types."""
    return '{0}: '.format(jt['name']) if jt['name'] else ''


def type_clusters(jt, tables):
    """The enabled clusters a job type uses."""
    return [c for c in tables.clusters
            if not jt.get('clusters') or c['cluster'] in jt['clusters']]


def type_cluster_value(jt, cluster, key, default):
    """A job type's clustmem or clustcpu on a cluster: the type's (one
    value, or {cluster: value}), else the cluster's (clusters.tsv), else
    default."""
    val = jt.get(key)
    if isinstance(val, dict):
        val = val.get(cluster['cluster'])
    if val is None or str(val).strip() == '':
        val = cluster.get(key)
    if val is None or str(val).strip() == '':
        val = default
    return val


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

    def __init__(self, rundir, types, script_name):
        self.rundir = rundir
        path = lambda name: os.path.join(rundir, name)
        self.devices = devmod.read_devices(path('devices.tsv'))
        self.assignments = devmod.read_assignments(path('assignments.tsv'))
        self.sessions = devmod.read_tsv(path('sessions.tsv'))
        self.clusters = [c for c in devmod.read_tsv(path('clusters.tsv'))
                         if _int(c.get('enabled'))]
        # jobs by (participant, session, job type); a row without job_type
        #   (jobs.tsv from before the types): the type whose prefix starts
        #   its jobname, else the top-level script_name's type, else the
        #   first; a row of a type not in the config is kept as it is
        names = [jt['name'] for jt in types]
        default = next((jt['name'] for jt in types
                        if jt['script_name'] == script_name), names[0])
        self.jobs = {}
        self.other_jobs = []
        for r in devmod.read_tsv(path('jobs.tsv')):
            name = r.get('job_type') or ''
            if name not in names:
                if name:
                    self.other_jobs.append(r)
                    continue
                name = next((jt['name'] for jt in types
                             if (r.get('jobname') or '').startswith(
                                 jt['prefix'] + '_')), default)
            r['job_type'] = name
            key = (r['participant_id'], r['session_id'], name)
            if key in self.jobs:
                raise ValueError('jobs.tsv: {0} {1} {2} listed twice'.format(
                    key[0], key[1], name or '(no job_type)'))
            self.jobs[key] = r

    def write(self, outdir, history_keep, stamp):
        os.makedirs(outdir, exist_ok=True)
        path = lambda name: os.path.join(outdir, name)
        devmod.write_devices(path('devices.tsv'), self.devices)
        devmod.write_assignments(path('assignments.tsv'), self.assignments)
        devmod.write_tsv(path('jobs.tsv'),
                         [self.jobs[k] for k in sorted(self.jobs)] +
                         self.other_jobs, JOB_COLUMNS)
        histdir = os.path.join(outdir, 'history')
        os.makedirs(histdir, exist_ok=True)
        for name in ['devices', 'assignments', 'jobs']:
            shutil.copy2(path(name + '.tsv'),
                         os.path.join(histdir, '{0}_{1}.tsv'.format(name, stamp)))
            old = sorted(f for f in os.listdir(histdir) if f.startswith(name + '_'))
            for f in old[:-history_keep] if history_keep > 0 else []:
                os.remove(os.path.join(histdir, f))

    def job(self, pid, sid, job_type):
        """The jobs.tsv row of a session's job type, created if new."""
        key = (pid, sid, job_type)
        if key not in self.jobs:
            self.jobs[key] = {c: '' for c in JOB_COLUMNS}
            self.jobs[key].update({'participant_id': pid, 'session_id': sid,
                                   'job_type': job_type, 'attempts': 0,
                                   'failures': 0, 'lost': 0, 'stale': 0,
                                   'oom': 0})
        return self.jobs[key]


def set_state(row, state, errmsg=None):
    row['state'] = state
    row['state_date'] = _now()
    if errmsg is not None:
        row['errmsg'] = (errmsg or '').replace('\t', ' ').replace('\n', ' ')


def device_dirs(config, jt, device):
    """The root dirs of a job type's script for one device."""
    dirs = {'raw_rootdir': config['raw_rootdir']}
    for key, sub in jt['device_roots'].items():
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


def not_processed(script_name, pid, sid, dirs, opts, rules):
    """True where the session has a failing status file and every one of
    them is a refusal (class not_processed): a rerun gives the same result
    until its inputs change."""
    classes = []
    for r in mkjobs.status_rows(script_name, pid, sid, dirs, opts):
        st = r['status']
        if st is None or str(st.get('errcode', '')) == '0':
            continue
        module = errcat.step_fields(r['step'])['module']
        classes.append(errcat.classify(st.get('errmsg', ''), module, rules)[0])
    return bool(classes) and all(c == 'not_processed' for c in classes)


def status_oom(script_name, pid, sid, dirs, opts, oom_re):
    """Where the errmsg of a failed status (errcode != 0) of the session
    matches oom_re (the out of memory resource rules), a reason naming its
    step, else None."""
    if oom_re is None:
        return None
    for r in mkjobs.status_rows(script_name, pid, sid, dirs, opts):
        st = r['status']
        if st is None or str(st.get('errcode', '')).strip() == '0':
            continue
        if oom_re.search(st.get('errmsg') or ''):
            return 'out of memory error in the {0} status'.format(r['step'])
    return None


def log_step(script_name):
    """The desc- label of a script's job logs: cmig_pipeline_bids_rawqc ->
    pipelinerawqc."""
    name = re.sub(r'^cmig_', '', script_name).replace('_bids', '')
    return re.sub(r'[^A-Za-z0-9]', '', name)


def job_log_base(row):
    """<batchdir>/pbsout/<job> of a row's job (the paths clu.submit_lines
    gives the scheduler, + .out and .err), or None."""
    batchdir = row.get('batchdir') or ''
    prefix = os.path.basename(os.path.normpath(batchdir)) + '_'
    jobname = row.get('jobname') or ''
    if not batchdir or not jobname.startswith(prefix):
        return None
    return os.path.join(batchdir, 'pbsout', jobname[len(prefix):])


def log_oom(row):
    """Where a job's .out or .err holds MATLAB's "Out of memory", a reason
    naming the file, else None."""
    src = job_log_base(row)
    if src is None:
        return None
    for ext in ('.out', '.err'):
        try:
            with open(src + ext, errors='replace') as fid:
                for line in fid:
                    if OOM_LOG_RE.search(line):
                        return '"Out of memory" in {0}{1}'.format(
                            os.path.basename(src), ext)
        except OSError:
            continue
    return None


def log_crash(row, script_name):
    """Where a job's .out shows it ran to its end and reported a failure
    (run.py's "N of M sessions failed" with N > 0, or the script's closing
    lines "<script_name>: errcode = N", N not 0, then "<script_name>:
    finished"), a reason, with the last ERROR line of its .err, else None.
    A job that never started, or was killed before its end, has no such
    lines."""
    src = job_log_base(row)
    if src is None:
        return None
    try:
        with open(src + '.out', errors='replace') as fid:
            out_text = fid.read()
    except OSError:
        return None
    ended = [m for m in RUNPY_END_RE.finditer(out_text) if int(m.group(1)) > 0]
    if not ended:
        name = re.escape(script_name)
        for m in re.finditer(r'^{0}: errcode = (\S+)\s*$'.format(name),
                             out_text, re.MULTILINE):
            if m.group(1) != '0' and re.search(
                    r'^{0}: finished\s*$'.format(name),
                    out_text[m.end():], re.MULTILINE):
                ended.append(m)
    if not ended:
        return None
    why = 'job ended with an error and no failed status'
    try:
        with open(src + '.err', errors='replace') as fid:
            errors = [ln.strip() for ln in fid if ERROR_LINE_RE.search(ln)]
        if errors:
            why += ': ' + errors[-1][:500]
    except OSError:
        pass
    return why


def copy_job_logs(row, script_name, dirs, log):
    """Copy a settled job's .out and .err into each output root's logs/.
    Done here rather than in the job script, since a job's .out is still open
    at its own end and a job the scheduler kills never reaches a last line."""
    src = job_log_base(row)
    if src is None or not os.path.isfile(src + '.out'):
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


def settle(tables, config, answers, devices_by_name, log, dry_run=False,
           oom_why=None):
    """Update the jobs in a queue; settle those that left one, each by its
    job type's script, roots, opts and delay flags. Returns (n_settled,
    n_adopted, the rows that ended done or failed, those that ended failed
    or lost). oom_why, a dict, gets {id(row): reason} for
    each failed row with an out of memory status errmsg (status_oom)."""
    types = {jt['name']: jt for jt in config['types']}
    n_settled = 0
    finished = []
    ended = []
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
        jt = types[row['job_type']]
        opts = jt['job_opts']
        dirs = device_dirs(config, jt, device)
        if not dry_run:
            copy_job_logs(row, jt['script_name'], dirs, log)
        status, errmsg = mkjobs.check_session(
            jt['script_name'], row['participant_id'], row['session_id'],
            dirs, opts, delay_flags=jt['delay_flags'])
        n_settled += 1
        if status == 'pending':
            # a job that ran and crashed before writing a status counts as
            #   failed; one that never ran, or was killed, stays lost
            crash = log_crash(row, jt['script_name'])
            if crash:
                status, errmsg = 'failed', crash
        if status in ('done', 'failed'):
            finished.append(row)
        if status != 'done' and status != 'stale':
            ended.append(row)
        if status == 'done':
            set_state(row, 'done', '')
        elif status == 'failed':
            row['failures'] = _int(row.get('failures')) + 1
            if oom_why is not None:
                why = status_oom(jt['script_name'], row['participant_id'],
                                 row['session_id'], dirs, opts,
                                 config.get('oom_re'))
                if why:
                    oom_why[id(row)] = why
            if row['failures'] > config['max_failures']:
                set_state(row, 'held', errmsg)
            elif not_processed(jt['script_name'], row['participant_id'],
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
            log('WARNING: {0}{1} {2}: job ended stale: {3}'.format(
                _tag(jt), row['participant_id'], row['session_id'], errmsg))
        else:
            row['lost'] = _int(row.get('lost')) + 1
            if row['lost'] > config['max_lost']:
                set_state(row, 'held', 'lost {0} times'.format(row['lost']))
            else:
                set_state(row, 'lost', 'job ended with status {0}'.format(status))
        log('{0}{1} {2}: {3}{4}'.format(_tag(jt), row['participant_id'],
                                        row['session_id'], row['state'],
                                        ': ' + row['errmsg'] if row.get('errmsg') else ''))
    # a job of this run in a queue but not in jobs.tsv (a round that crashed
    #   after submitting, or a job submitted by hand) is adopted, as the job
    #   type whose prefix starts its name
    n_adopted = 0
    known = set((r.get('cluster'), r.get('jobid'))
                for r in list(tables.jobs.values()) + tables.other_jobs)
    for (cluster, jobid), (name, state) in sorted(in_queue.items()):
        if (cluster, jobid) in known:
            continue
        jt = next((t for t in config['types']
                   if name.startswith(t['prefix'] + '_')), None)
        if jt is None:
            continue
        m = JOBNAME_RE.search(name)
        if not m or state not in IN_QUEUE:
            continue
        row = tables.job(m.group(1), m.group(2), jt['name'])
        if row.get('state') in IN_QUEUE:
            continue
        row.update({'cluster': cluster, 'jobid': jobid, 'jobname': name})
        # the adoption time stands in for the unknown submit time, so
        #   reopen_changed has a date to compare rawdata with
        if not row.get('submit_date'):
            row['submit_date'] = _now()
        set_state(row, state)
        if not row.get('device'):
            a = tables.assignments.get(m.group(1))
            row['device'] = a['device'] if a else ''
        n_adopted += 1
        log('{0}adopted {1} on {2} ({3})'.format(_tag(jt), name, cluster, jobid))
    return n_settled, n_adopted, finished, ended


def record_usage(tables, finished, log, ended=()):
    """Write the accounting of each job that finished (ended done or
    failed) or ended lost into its jobs.tsv row, one query per cluster;
    return the per-cluster summary columns of rounds.tsv (over the finished
    jobs) and {(cluster, jobid): usage}."""
    sched = {c['cluster']: c['scheduler'] for c in tables.clusters}
    by_cluster = {}
    seen = set()
    for row in list(finished) + list(ended):
        if id(row) in seen:
            continue
        seen.add(id(row))
        by_cluster.setdefault(row['cluster'], []).append(row)
    fin = set(id(r) for r in finished)
    summary = {k: [] for k in ['finished', 'mem_gb_median', 'cores_median',
                               'write_gb_total', 'read_gb_total']}
    usage_all = {}
    for cluster in sorted(by_cluster):
        rows_all = by_cluster[cluster]
        rows = [r for r in rows_all if id(r) in fin]
        usage = None
        if cluster in sched:
            try:
                # a lost job's id may be blank or not a scheduler's
                usage = clu.job_usage(cluster, sched[cluster],
                                      [r['jobid'] for r in rows_all
                                       if id(r) in fin or
                                       re.match(r'^\d+([._]\d+)?$', r['jobid'] or '')])
            except ValueError as e:
                log('{0}: accounting not queried: {1}'.format(cluster, e))
        if usage is None:
            if rows:
                log('{0}: no accounting for {1} finished jobs'.format(
                    cluster, len(rows)))
            usage = {}
        found = [usage[r['jobid']] for r in rows if r['jobid'] in usage]
        for row in rows_all:
            u = usage.get(row['jobid'], {})
            if u:
                usage_all[(cluster, row['jobid'])] = u
            for col, key, fmt in USAGE_COLUMNS:
                row[col] = '' if u.get(key) is None else fmt.format(u[key])
        if not rows:
            continue

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
    return {k: ' '.join(v) for k, v in summary.items()}, usage_all


def oom_exponent(row, config):
    """The power of retry_mem_factor a row's next job gets: its out of
    memory endings, up to max_failures."""
    return min(_int(row.get('oom')), max(_int(config['max_failures']), 0))


def mark_oom(ended, usage, config, log, oom_why=None):
    """Count in oom each job that ended failed or lost out of memory (the
    accounting's oom: slurm OUT_OF_MEMORY or SGE exit_status 137; else
    "Out of memory" in its .out or .err; else an out of memory status errmsg,
    oom_why from settle), noting it in errmsg with the memory factor of its
    retry; returns their number."""
    n = 0
    for row in ended:
        u = usage.get((row.get('cluster'), row.get('jobid'))) or {}
        why = u.get('oom') or log_oom(row) or (oom_why or {}).get(id(row))
        if not why:
            continue
        row['oom'] = _int(row.get('oom')) + 1
        note = 'out of memory ({0})'.format(why)
        if row.get('state') != 'held':
            note += ': next job at clustmem x {0:g}'.format(
                config['retry_mem_factor'] ** oom_exponent(row, config))
        row['errmsg'] = '; '.join(m for m in (row.get('errmsg'), note) if m)
        log('{0}{1} {2}: {3}'.format(
            row['job_type'] + ': ' if row.get('job_type') else '',
            row['participant_id'], row['session_id'], note))
        n += 1
    return n


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
    after it was last found done without a job (checked_date), or, with
    neither (a job adopted before adoption set submit_date), after it became
    done (state_date), so mkjobs compares it again; returns their number."""
    n = 0
    for (pid, sid, _), row in sorted(tables.jobs.items()):
        if row.get('state') != 'done':
            continue
        refs = [t for t in (_parse_time(row.get('submit_date')),
                            _parse_time(row.get('checked_date'))) if t is not None]
        if not refs:
            t = _parse_time(row.get('state_date'))
            if t is None:
                continue
            refs = [t]
        newest = _rawdata_mtime(config['raw_rootdir'], pid, sid)
        # the dates are kept to the second
        if newest is not None and newest.replace(microsecond=0) > max(refs):
            set_state(row, 'changed', 'rawdata scans.tsv or events changed {0}'.format(
                newest.strftime(TIME_FORMAT)))
            n += 1
    if n:
        log('{0} done sessions returned to the candidates: rawdata changed'.format(n))
    return n


def candidates(tables, config, jt, eligible):
    """The (participant, session) pairs to consider for a job type, in the
    config order."""
    pairs = [(pid, sid) for pid in eligible for sid in eligible[pid]]
    if config['order'] == 'session':
        pairs.sort(key=lambda p: (p[1], p[0]))
    else:
        pairs.sort()
    return [p for p in pairs
            if tables.jobs.get((p[0], p[1], jt['name']), {}).get('state')
            not in NOT_CANDIDATE]


def make_jobs(tables, config, jt, quota_total, batchdir, df_failed,
              running_by_device, no_device, image_info, log):
    """Write a job type's job scripts for up to quota_total candidates;
    return their [(job_name, row)] and the number of participants assigned.
    no_device (shared by the types) collects the participants no device
    took; image_info ({image: (sha256, versions)}) shares an image's
    sha256 between the types."""
    script_name = jt['script_name']
    opts = jt['job_opts']
    eligible = jt['eligible']
    delay_flags = jt['delay_flags']
    devices_by_name = {d['device']: d for d in tables.devices}
    checked = set()
    made = []
    n_assigned = 0
    pairs = candidates(tables, config, jt, eligible)
    # the image's sha256 and versions, once per round with jobs to make
    if pairs and quota_total > 0:
        image = opts.get('image')
        if image in image_info and 'image_sha256' not in opts:
            opts['image_sha256'], opts['image_versions'] = image_info[image]
        mkjobs.add_image_info(opts, log)
        if image:
            image_info[image] = (opts.get('image_sha256'),
                                 opts.get('image_versions'))
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
                continue
            n_assigned += 1
        device = devices_by_name[tables.assignments[pid]['device']]
        if device['device'] in df_failed:
            continue
        dirs = device_dirs(config, jt, device)
        if device['device'] not in checked:
            mkjobs.check_dirs(script_name, dirs)
            # kept where they exist, so only the first round writes them
            for msg in mkjobs.write_root_files(
                    script_name, dirs, dataset_name=opts.get('dataset_name'),
                    FS_version=opts['FS_version'],
                    versions=opts.get('image_versions')):
                log('{0}{1}: {2}'.format(_tag(jt), device['device'], msg))
            checked.add(device['device'])
        row = tables.job(pid, sid, jt['name'])
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
    return made, n_assigned


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


def submit(tables, config, jt, batchdir, by_cluster, dry_run, log):
    """Submit each cluster's jobs of a job type, at the type's clustmem and
    clustcpu there (else the cluster's), a job that ran out of memory before
    at clustmem x retry_mem_factor^oom; record those with a job id as
    queued."""
    n = 0
    tag = _tag(jt)
    clusters = {c['cluster']: c for c in tables.clusters}
    batchname = os.path.basename(batchdir)
    for cluster, items in sorted(by_cluster.items()):
        if not items:
            continue
        c = clusters[cluster]
        names = [j for j, _ in items]
        with open(os.path.join(batchdir, 'scriptlist_{0}.txt'.format(cluster)), 'w') as fid:
            fid.write(''.join(j + '\n' for j in names))
        base = str(type_cluster_value(jt, c, 'clustmem', '60G'))
        mems = {}
        for job_name, row in items:
            k = oom_exponent(row, config)
            mems[job_name] = clu.scale_mem(base, config['retry_mem_factor'] ** k) \
                if k else base
        ids = clu.submit_jobs(cluster, c['scheduler'], batchdir, names, mems,
                              _int(type_cluster_value(jt, c, 'clustcpu', 1), 1),
                              cstr=c.get('cstr') or 'qcshjob', batchname=batchname,
                              dry_run=dry_run)
        n_id = 0
        for job_name, row in items:
            jobid = ids.get(job_name)
            if jobid is None and not dry_run:
                log('{0}{1}: no job id from {2}; back to the pool'.format(
                    tag, job_name, cluster))
                continue
            n_id += 1
            row.update({'cluster': cluster, 'jobid': jobid or 'dry_run',
                        'jobname': '{0}_{1}'.format(batchname, job_name),
                        'batchdir': batchdir, 'clustmem': mems[job_name],
                        'submit_date': _now()})
            # the previous job's accounting
            row.update({col: '' for col, _, _ in USAGE_COLUMNS})
            row['attempts'] = _int(row.get('attempts')) + 1
            set_state(row, 'queued')
            n += 1
        if dry_run:
            log('{0}{1}: {2} jobs (dry run)'.format(tag, cluster, len(items)))
        else:
            log('{0}{1}: {2} jobs submitted{3}'.format(
                tag, cluster, n_id, ', {0} without a job id'.format(len(items) - n_id)
                if n_id < len(items) else ''))
        raised = sorted(set(m for m in mems.values() if m != base))
        if raised:
            log('{0}{1}: {2} at a raised memory limit ({3}; clustmem {4}) after '
                'running out of memory'.format(
                    tag, cluster, sum(1 for m in mems.values() if m != base),
                    ' '.join(raised), base))
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

    types = config['types']
    tables = Tables(rundir, types, config['script_name'])
    if tables.other_jobs:
        log('WARNING: {0} jobs.tsv rows of job types not in config.json (kept '
            'as they are): {1}'.format(len(tables.other_jobs), ' '.join(sorted(
                set(r['job_type'] for r in tables.other_jobs)))))
    ready_pairs = None
    if config['ready_file']:
        ready_pairs = devmod.read_ready_file(config['ready_file'])
        log('{0} sessions in the ready file {1}'.format(
            len(ready_pairs), config['ready_file']))
    for jt in types:
        script_name = jt['script_name']
        jt['job_opts'] = mkjobs.make_opts(**jt['opts'])
        jt['delay_flags'] = mkjobs.load_delay_flags(script_name,
                                                    config['tracking_csv'])
        # the script's required input dirs that sessions.tsv records
        required = [m for m in
                    mkjobs.SCRIPT_CONFIG[script_name].get('required_subdirs', [])
                    if m in ('anat', 'func', 'dwi', 'fmap', 'beh')]
        jt['eligible'] = devmod.eligible_sessions(
            tables.sessions, config['seslist'], required,
            bool(config['ready_only']), ready_pairs)
        log('{0}{1} eligible sessions of {2} participants'.format(
            _tag(jt), sum(len(v) for v in jt['eligible'].values()),
            len(jt['eligible'])))

    answers = query_clusters(tables, config, log)
    devices_by_name = {d['device']: d for d in tables.devices}
    oom_why = {}
    n_settled, n_adopted, finished, ended = settle(
        tables, config, answers, devices_by_name, log, dry_run=dry_run,
        oom_why=oom_why)
    usage_summary, usage = record_usage(tables, finished, log, ended)
    mark_oom(ended, usage, config, log, oom_why)
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

    # the types in the order listed, sharing the quotas and the round's limit
    left = dict(quotas)
    no_device = set()
    image_info = {}
    n_assigned = n_submitted = 0
    for jt in types:
        tq = {c['cluster']: left[c['cluster']] for c in type_clusters(jt, tables)
              if c['cluster'] in left}
        type_total = min(sum(tq.values()), quota_total)
        if len(types) > 1:
            log('{0}quota: {1} ({2})'.format(_tag(jt), type_total, ' '.join(
                '{0}:{1}'.format(c, q) for c, q in sorted(tq.items()))))
        batchdir = os.path.join(config['batchroot'], '{0}_{1}{2}'.format(
            jt['prefix'], stamp, '_dry_run' if dry_run else ''))
        os.makedirs(batchdir, exist_ok=True)
        made, n = make_jobs(tables, config, jt, type_total, batchdir, df_failed,
                            running_by_device, no_device, image_info, log)
        n_assigned += n
        by_cluster = distribute(made, tq)
        n_submitted += submit(tables, config, jt, batchdir, by_cluster, dry_run, log)
        for c, items in by_cluster.items():
            left[c] -= len(items)
            quota_total -= len(items)
    n_unassigned = len(no_device)
    if no_device:
        log('{0} participants found no device (max_frac {1}); their sessions '
            'wait'.format(len(no_device), config['max_frac']))

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
