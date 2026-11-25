#!/usr/bin/env python3
"""
Query, submit and cancel cluster jobs over ssh to a head node, for slurm
and SGE clusters.

Jobs are submitted with the same per-job sbatch/qsub lines as MMPS's qjobs
(job name <batchname>_<job>, output in <batchdir>/pbsout, the job wrapper,
e.g. qcshjob), all jobs for a cluster in one ssh call, and each job's id is
read from its submit output. On slurm, --cpus-per-task is requested. On
SGE, a job with clustcpu > 1 asks for that many slots (-pe smp) and is bound
to as many cores (-binding linear), with h_vmem (a per-slot consumable) at
clustmem / clustcpu and mem_free (a host free-memory check) at clustmem; a
1-CPU job is submitted as in qjobs. clustmem may differ from job to job (a
retry at a higher memory limit). ssh never reads the caller's stdin.

Job states from query_jobs are normalized to 'queued', 'running' or
'other:<raw state>'.

job_usage reads each job's memory, CPU and disk I/O from the scheduler's
accounting (sacct on slurm, qacct on SGE), one ssh per cluster, and whether
the scheduler ended it for its memory (oom: slurm state OUT_OF_MEMORY, SGE
exit_status 137).

Created:  10/29/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Last mod: 11/18/2026 by Don Hagler
"""

import argparse
import getpass
import math
import os
import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET

SCHEDULERS = ['slurm', 'sge']
SSH_OPTIONS = ['-o', 'BatchMode=yes', '-o', 'ConnectTimeout=30']
MARKER = '@@JOB'

SLURM_QUEUED = {'PENDING', 'CONFIGURING', 'REQUEUED', 'REQUEUE_FED',
                'REQUEUE_HOLD', 'RESIZING'}
SLURM_RUNNING = {'RUNNING', 'COMPLETING'}
SGE_RUNNING = {'r', 't', 'Rr', 'Rt'}

# compute node host name -> cluster (head node)
HOST_CLUSTERS = [(r'^101[34]s-', 'brightcluster1'),
                 (r'^mmil-compute-8-', 'mmilcluster8'),
                 (r'^mmil-compute-9-', 'mmilcluster9'),
                 (r'^mmil-compute-10-', 'mmilcluster10')]
CLUSTER_SCHEDULERS = {'brightcluster1': 'slurm', 'mmilcluster8': 'sge',
                      'mmilcluster9': 'sge', 'mmilcluster10': 'sge'}

# job_usage fields, per job
USAGE_FIELDS = ['elapsed_h', 'cpu_h', 'cores', 'mem_gb', 'mem_type',
                'read_gb', 'write_gb', 'node', 'state', 'exit_status', 'failed',
                'oom']
SACCT_FIELDS = ['JobID', 'JobIDRaw', 'State', 'Elapsed', 'NodeList',
                'TotalCPU', 'MaxRSS', 'MaxDiskRead', 'MaxDiskWrite']
# sacct and qacct sizes are in multiples of 1024; GB here is 1024^3 bytes
SIZE_UNITS = {'': 1.0, 'B': 1.0, 'K': 1024.0, 'M': 1024.0 ** 2,
              'G': 1024.0 ** 3, 'T': 1024.0 ** 4, 'P': 1024.0 ** 5}
GB = 1024.0 ** 3
# SGE memory units: upper case multiples of 1024, lower case of 1000
SGE_MEM_UNITS = ['', 'K', 'M', 'G', 'T']


def run_ssh(host, command, timeout=120):
    """Run command on host under its login shell (tcsh);
    return (returncode, stdout, stderr), with returncode None on timeout
    or when ssh fails to start or to connect (ssh exit status 255). stdin
    is /dev/null: ssh reading a terminal's stdin stops a background caller
    (SIGTTIN)."""
    try:
        res = subprocess.run(['ssh'] + SSH_OPTIONS + [host, command],
                             stdin=subprocess.DEVNULL,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, '', 'ssh timed out after {0} s'.format(timeout)
    except OSError as e:
        return None, '', 'ssh failed to start: {0}'.format(e)
    if res.returncode == 255:
        return None, res.stdout, res.stderr
    return res.returncode, res.stdout, res.stderr


def _bash_command(lines):
    """One bash -c command for the script lines, quoted for the remote
    tcsh (single line: tcsh rejects newlines inside quotes)."""
    return 'bash -c {0}'.format(shlex.quote('; '.join(lines)))


def _slurm_state(raw):
    if raw in SLURM_RUNNING:
        return 'running'
    if raw in SLURM_QUEUED:
        return 'queued'
    return 'other:{0}'.format(raw)


def _sge_state(raw):
    if raw in SGE_RUNNING:
        return 'running'
    # e.g. qw, hqw, hRwq; error, deletion and suspended states excluded
    if 'q' in raw and not set(raw) & set('EdsST'):
        return 'queued'
    return 'other:{0}'.format(raw)


def query_jobs(cluster, scheduler, user):
    """{jobid: (jobname, state)} for the user's jobs on the cluster, or None
    when the head node does not answer or the query fails."""
    if scheduler == 'slurm':
        cmd = "squeue -h -u {0} -o '%i|%j|%T'".format(user)
    elif scheduler == 'sge':
        cmd = 'qstat -u {0} -xml'.format(user)
    else:
        raise ValueError('unknown scheduler: {0}'.format(scheduler))
    rc, out, _ = run_ssh(cluster, cmd)
    if rc != 0:
        return None
    jobs = {}
    if scheduler == 'slurm':
        for line in out.splitlines():
            parts = line.strip().split('|')
            if len(parts) != 3:
                continue
            jobid, name, raw = parts
            jobs[jobid] = (name, _slurm_state(raw))
        return jobs
    try:
        root = ET.fromstring(out)
    except ET.ParseError:
        return None
    for elem in root.iter('job_list'):
        jobid = elem.findtext('JB_job_number', '').strip()
        if not jobid:
            continue
        # array tasks share a job number: one entry per task
        tasks = elem.findtext('tasks', '').strip()
        if tasks:
            jobid = '{0}.{1}'.format(jobid, tasks)
        name = elem.findtext('JB_name', '').strip()
        raw = elem.findtext('state', '').strip()
        jobs[jobid] = (name, _sge_state(raw))
    return jobs


def count_jobs(jobs):
    """(n_running, n_queued) over a query_jobs result."""
    states = [s for _, s in jobs.values()]
    return states.count('running'), states.count('queued')


def sge_mem_per_slot(clustmem, nslots):
    """clustmem (an SGE memory value: number and optional unit K, M, G or T,
    upper case 1024-based, lower case 1000-based, none bytes) divided by
    nslots, in the same unit when whole, else in the next smaller unit,
    rounded up. Raises ValueError."""
    m = re.match(r'^\s*([0-9]*\.?[0-9]+)\s*([kKmMgGtT]?)\s*$', str(clustmem))
    if not m:
        raise ValueError('invalid clustmem: {0}'.format(clustmem))
    unit = m.group(2)
    base = 1000 if unit.islower() else 1024
    i = SGE_MEM_UNITS.index(unit.upper())
    val = float(m.group(1)) / nslots
    if abs(val - round(val)) > 1e-9 and i > 0:
        val *= base
        i -= 1
    if abs(val - round(val)) < 1e-9:
        val = round(val)
    u = SGE_MEM_UNITS[i]
    return '{0}{1}'.format(int(math.ceil(val)),
                           u.lower() if unit.islower() else u)


def scale_mem(clustmem, factor):
    """clustmem (number and optional unit K, M, G or T, either case) times
    factor, in the same unit, rounded up. Raises ValueError."""
    m = re.match(r'^\s*([0-9]*\.?[0-9]+)\s*([kKmMgGtT]?)\s*$', str(clustmem))
    if not m:
        raise ValueError('invalid clustmem: {0}'.format(clustmem))
    val = float(m.group(1)) * factor
    if abs(val - round(val)) < 1e-9:
        val = round(val)
    return '{0}{1}'.format(int(math.ceil(val)), m.group(2))


def submit_lines(scheduler, batchdir, job_names, clustmem, clustcpu=1,
                 cstr='qcshjob', batchname=None):
    """The remote script lines that submit the jobs, as qjobs does; on SGE,
    with clustcpu > 1, -pe smp and -binding linear for clustcpu slots and
    h_vmem per slot (clustmem / clustcpu). clustmem: one value for all jobs,
    or {job_name: value}. Raises ValueError."""
    clustcpu = int(clustcpu or 1)
    pe = ''
    if scheduler == 'sge' and clustcpu > 1:
        pe = ' -pe smp {0} -binding linear:{0}'.format(clustcpu)
    if not batchname:
        batchname = os.path.basename(os.path.normpath(batchdir))
    lines = ['cd {0}'.format(shlex.quote(batchdir)), 'mkdir -p pbsout']
    for job in job_names:
        name = '{0}_{1}'.format(batchname, job)
        out = 'pbsout/{0}.out'.format(job)
        err = 'pbsout/{0}.err'.format(job)
        wrap = 'echo "{0} {1} {2}"'.format(cstr, batchdir, job)
        mem = clustmem.get(job) if isinstance(clustmem, dict) else clustmem
        if scheduler == 'slurm':
            sub = ('{0} | xargs -I {{}} sbatch --wrap="{{}}" --mem={1} '
                   '--job-name={2} --output={3} --error={4} --nodes=1 '
                   '--ntasks-per-node=1 --cpus-per-task={5}').format(
                       wrap, mem, name, out, err, clustcpu)
        elif scheduler == 'sge':
            vmem = sge_mem_per_slot(mem, clustcpu) if pe else mem
            sub = ('{0} | qsub -V -hard -l mem_free={1},h_vmem={2}{3} -cwd '
                   '-N {4} -o {5} -e {6}').format(wrap, mem, vmem, pe,
                                                  name, out, err)
        else:
            raise ValueError('unknown scheduler: {0}'.format(scheduler))
        lines += ['echo {0} {1}'.format(MARKER, job),
                  'rm -f {0} {1}'.format(out, err),
                  '{0} 2>&1'.format(sub),
                  'sleep 0.25']
    return lines


def parse_submit_output(out):
    """{job_name: jobid} from submit output, by the marker line before each
    job's submit."""
    ids = {}
    job = None
    for line in out.splitlines():
        line = line.strip()
        if line.startswith(MARKER + ' '):
            job = line[len(MARKER) + 1:]
            continue
        if job is None:
            continue
        m = (re.match(r'Submitted batch job (\d+)', line) or
             re.match(r'Your job (\d+) \(".*"\) has been submitted', line))
        if m:
            ids[job] = m.group(1)
            job = None
    return ids


def submit_jobs(cluster, scheduler, batchdir, job_names, clustmem, clustcpu=1,
                cstr='qcshjob', batchname=None, dry_run=False):
    """Submit batchdir/<job>.csh for each job in one ssh call (clustmem: one
    value, or {job_name: value}); return {job_name: jobid} for the jobs whose
    submit reported an id.
    With dry_run, print the remote command and return {}."""
    if not job_names:
        return {}
    lines = submit_lines(scheduler, batchdir, job_names, clustmem, clustcpu,
                         cstr, batchname)
    if dry_run:
        print('ssh {0} bash -c, with the lines:'.format(cluster))
        for line in lines:
            print('  {0}'.format(line))
        return {}
    timeout = 120 + len(job_names)
    rc, out, err = run_ssh(cluster, _bash_command(lines), timeout=timeout)
    if rc is None:
        print('WARNING: submit to {0} failed: {1}'.format(
            cluster, err.strip()), file=sys.stderr)
    ids = parse_submit_output(out)
    missing = [j for j in job_names if j not in ids]
    if missing:
        print('WARNING: no job id for {0} of {1} jobs on {2}'.format(
            len(missing), len(job_names), cluster), file=sys.stderr)
    return ids


def cancel_jobs(cluster, scheduler, jobids, dry_run=False):
    """Cancel the jobs (scancel / qdel) in one ssh call; return True on
    success. With dry_run, print the remote command only."""
    jobids = [str(j) for j in jobids]
    if not jobids:
        return True
    bad = [j for j in jobids if not re.match(r'^\d+(\.\d+)?$', j)]
    if bad:
        raise ValueError('invalid job ids: {0}'.format(' '.join(bad)))
    if scheduler == 'slurm':
        cmd = 'scancel {0}'.format(' '.join(jobids))
    elif scheduler == 'sge':
        cmd = 'qdel {0}'.format(' '.join(jobids))
    else:
        raise ValueError('unknown scheduler: {0}'.format(scheduler))
    if dry_run:
        print('ssh {0} {1}'.format(cluster, cmd))
        return True
    rc, _, err = run_ssh(cluster, cmd)
    if rc != 0:
        print('WARNING: cancel on {0} failed: {1}'.format(
            cluster, err.strip()), file=sys.stderr)
        return False
    return True


def cluster_from_host(host):
    """The cluster of a compute node host name, or None."""
    for pattern, cluster in HOST_CLUSTERS:
        if re.match(pattern, host or ''):
            return cluster
    return None


def _seconds(text):
    """Seconds from [D-][HH:]MM:SS[.sss] (sacct) or N[.n][s] (qacct);
    None if empty or unparsed."""
    text = (text or '').strip()
    if not text:
        return None
    try:
        if ':' not in text:
            return float(text.rstrip('s'))
        days = 0
        if '-' in text:
            d, text = text.split('-', 1)
            days = int(d)
        parts = [float(p) for p in text.split(':')]
    except ValueError:
        return None
    while len(parts) < 3:
        parts.insert(0, 0.0)
    return days * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def _gb(text):
    """GB from a size such as 8489416K, 107367.17M (sacct) or 11.689GB
    (qacct); None if empty or unparsed."""
    m = re.match(r'^([0-9.]+)\s*([KMGTP]?)B?$', (text or '').strip())
    if not m:
        return None
    return float(m.group(1)) * SIZE_UNITS[m.group(2)] / GB


def _usage_record(elapsed, cpu, mem_gb, mem_type, read_gb, write_gb, node,
                  state=None, exit_status=None, failed=None, oom=None):
    rec = dict.fromkeys(USAGE_FIELDS)
    rec.update({'mem_gb': mem_gb, 'mem_type': mem_type, 'read_gb': read_gb,
                'write_gb': write_gb, 'node': node, 'state': state,
                'exit_status': exit_status, 'failed': failed, 'oom': oom})
    if elapsed is not None:
        rec['elapsed_h'] = elapsed / 3600.0
    if cpu is not None:
        rec['cpu_h'] = cpu / 3600.0
    if elapsed and cpu is not None:
        rec['cores'] = cpu / elapsed
    return rec


def _slurm_usage(cluster, jobids):
    """job_usage for slurm: the job row's Elapsed, State and NodeList, the
    .batch row's TotalCPU, MaxRSS, MaxDiskRead and MaxDiskWrite (empty on
    the job row)."""
    cmd = 'sacct -j {0} -P -o {1}'.format(','.join(jobids),
                                          ','.join(SACCT_FIELDS))
    rc, out, _ = run_ssh(cluster, cmd)
    if rc != 0:
        return None
    rows = {}
    for line in out.splitlines():
        parts = line.strip().split('|')
        if len(parts) != len(SACCT_FIELDS) or parts[0] == 'JobID':
            continue
        row = dict(zip(SACCT_FIELDS, parts))
        # an array task is shown as <array>_<task>, its JobIDRaw is its id
        rows[row['JobID']] = row
        rows[row['JobIDRaw']] = row
    usage = {}
    for jobid in jobids:
        job = rows.get(jobid)
        if job is None:
            continue
        batch = rows.get(job['JobIDRaw'] + '.batch', {})
        # the job row or its .batch step may carry OUT_OF_MEMORY
        oom = None
        if any((s or '').startswith('OUT_OF_MEMORY')
               for s in (job['State'], batch.get('State'))):
            oom = 'slurm OUT_OF_MEMORY'
        usage[jobid] = _usage_record(
            _seconds(job['Elapsed']), _seconds(batch.get('TotalCPU')),
            _gb(batch.get('MaxRSS')), 'rss', _gb(batch.get('MaxDiskRead')),
            _gb(batch.get('MaxDiskWrite')), job['NodeList'],
            state=job['State'], oom=oom)
    return usage


def _parse_qacct(text):
    """The last record of qacct -j output, as {field: value}."""
    rec = {}
    for line in text.splitlines():
        if line.startswith('====='):
            rec = {}
            continue
        m = re.match(r'^(\w+)\s+(.*?)\s*$', line)
        if m:
            rec[m.group(1)] = m.group(2)
    return rec


def _sge_usage(cluster, jobids):
    """job_usage for SGE: qacct -j per job in one ssh; memory is maxvmem,
    write is ru_oublock (512-byte blocks), no read."""
    lines = []
    for jobid in jobids:
        job, _, task = jobid.partition('.')
        lines += ['echo {0} {1}'.format(MARKER, jobid),
                  'qacct -j {0}{1} 2>&1'.format(
                      job, ' -t {0}'.format(task) if task else '')]
    lines.append('exit 0')
    rc, out, _ = run_ssh(cluster, _bash_command(lines),
                         timeout=120 + 5 * len(jobids))
    if rc != 0:
        return None
    chunks = {}
    jobid = None
    for line in out.splitlines():
        if line.startswith(MARKER + ' '):
            jobid = line[len(MARKER) + 1:].strip()
            chunks[jobid] = []
        elif jobid is not None:
            chunks[jobid].append(line)
    usage = {}
    for jobid, chunk in chunks.items():
        rec = _parse_qacct('\n'.join(chunk))
        if rec.get('jobnumber') != jobid.partition('.')[0]:
            continue
        try:
            write_gb = int(rec.get('ru_oublock', '')) * 512 / GB
        except ValueError:
            write_gb = None
        # 137: killed (SIGKILL), as when the job passes its h_vmem
        exit_status = rec.get('exit_status')
        oom = 'SGE exit_status 137' if (exit_status or '').split(' ')[0] == '137' \
            else None
        usage[jobid] = _usage_record(
            _seconds(rec.get('ru_wallclock')), _seconds(rec.get('cpu')),
            _gb(rec.get('maxvmem')), 'vmem', None, write_gb,
            rec.get('hostname'), exit_status=exit_status,
            failed=rec.get('failed'), oom=oom)
    return usage


def job_usage(cluster, scheduler, jobids):
    """{jobid: usage} from the scheduler's accounting, one ssh call; None
    when the head node does not answer or the query fails; jobs it does not
    know are absent. usage holds elapsed_h, cpu_h, cores (cpu / elapsed),
    mem_gb with mem_type ('rss', slurm MaxRSS; 'vmem', SGE maxvmem), read_gb
    (slurm only), write_gb, node, state (slurm) or exit_status and failed
    (SGE), and oom (why the scheduler ended the job for its memory: 'slurm
    OUT_OF_MEMORY' or 'SGE exit_status 137', else None); a value accounting
    does not give is None."""
    jobids = sorted(set(str(j) for j in jobids if str(j).strip()))
    if not jobids:
        return {}
    bad = [j for j in jobids if not re.match(r'^\d+([._]\d+)?$', j)]
    if bad:
        raise ValueError('invalid job ids: {0}'.format(' '.join(bad)))
    if scheduler == 'slurm':
        return _slurm_usage(cluster, jobids)
    if scheduler == 'sge':
        return _sge_usage(cluster, jobids)
    raise ValueError('unknown scheduler: {0}'.format(scheduler))


def parse_args():
    parser = argparse.ArgumentParser(
        description='Query cluster jobs, or show the commands that would '
                    'submit or cancel them.')
    parser.add_argument('--cluster', required=True,
                        help='Head node name (e.g. brightcluster1)')
    parser.add_argument('--scheduler', required=True, choices=SCHEDULERS,
                        help='Scheduler of the cluster')
    parser.add_argument('--user', default=getpass.getuser(),
                        help='User whose jobs are queried (default: current user)')
    parser.add_argument('--query', action='store_true',
                        help='Print the jobs (id, state, name) and the counts')
    parser.add_argument('--usage', action='store_true',
                        help='Print the accounting of --jobids (memory, CPU, '
                             'disk I/O)')
    parser.add_argument('--submit_dry_run', action='store_true',
                        help='Print the remote command that would submit --jobs')
    parser.add_argument('--cancel_dry_run', action='store_true',
                        help='Print the remote command that would cancel --jobids')
    parser.add_argument('--batchdir', default=None,
                        help='Batch directory holding <job>.csh')
    parser.add_argument('--jobs', nargs='+', default=None,
                        help='Job names (without .csh)')
    parser.add_argument('--jobids', nargs='+', default=None,
                        help='Job ids to cancel, or to report with --usage')
    parser.add_argument('--clustmem', default='20G',
                        help='Memory per job (default: 20G)')
    parser.add_argument('--clustcpu', type=int, default=1,
                        help='CPUs per job (SGE: slots, with -pe smp and '
                             '-binding when > 1) (default: 1)')
    parser.add_argument('--cstr', default='qcshjob',
                        help='Job wrapper command (default: qcshjob)')
    parser.add_argument('--batchname', default=None,
                        help='Job name prefix (default: basename of batchdir)')
    return parser.parse_args()


def run(args):
    if not (args.query or args.usage or args.submit_dry_run or
            args.cancel_dry_run):
        errmsg = ('one of --query, --usage, --submit_dry_run, '
                  '--cancel_dry_run required')
        print('ERROR: {0}'.format(errmsg), file=sys.stderr)
        return 1, errmsg
    if args.query:
        jobs = query_jobs(args.cluster, args.scheduler, args.user)
        if jobs is None:
            errmsg = 'query failed on {0}'.format(args.cluster)
            print('ERROR: {0}'.format(errmsg), file=sys.stderr)
            return 1, errmsg
        for jobid in sorted(jobs, key=lambda j: (len(j.split('.')[0]), j)):
            name, state = jobs[jobid]
            print('{0}\t{1}\t{2}'.format(jobid, state, name))
        nrun, nqueue = count_jobs(jobs)
        print('{0}: {1} jobs, {2} running, {3} queued, {4} other'.format(
            args.cluster, len(jobs), nrun, nqueue, len(jobs) - nrun - nqueue))
    if args.usage:
        if not args.jobids:
            errmsg = '--usage requires --jobids'
            print('ERROR: {0}'.format(errmsg), file=sys.stderr)
            return 1, errmsg
        usage = job_usage(args.cluster, args.scheduler, args.jobids)
        if usage is None:
            errmsg = 'accounting query failed on {0}'.format(args.cluster)
            print('ERROR: {0}'.format(errmsg), file=sys.stderr)
            return 1, errmsg
        print('\t'.join(['jobid'] + USAGE_FIELDS))
        for jobid in sorted(usage):
            vals = [usage[jobid][k] for k in USAGE_FIELDS]
            print('\t'.join([jobid] + ['' if v is None else
                                        '{0:.3f}'.format(v) if isinstance(v, float)
                                        else str(v) for v in vals]))
        missing = [j for j in args.jobids if j not in usage]
        if missing:
            print('not in accounting: {0}'.format(' '.join(missing)))
    if args.submit_dry_run:
        if not args.batchdir or not args.jobs:
            errmsg = '--submit_dry_run requires --batchdir and --jobs'
            print('ERROR: {0}'.format(errmsg), file=sys.stderr)
            return 1, errmsg
        submit_jobs(args.cluster, args.scheduler, args.batchdir, args.jobs,
                    args.clustmem, args.clustcpu, args.cstr, args.batchname,
                    dry_run=True)
    if args.cancel_dry_run:
        if not args.jobids:
            errmsg = '--cancel_dry_run requires --jobids'
            print('ERROR: {0}'.format(errmsg), file=sys.stderr)
            return 1, errmsg
        cancel_jobs(args.cluster, args.scheduler, args.jobids, dry_run=True)
    return 0, None


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
