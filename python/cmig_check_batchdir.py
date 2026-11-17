#!/usr/bin/env python3
"""
Report on a batch job directory: job accounting, run times and error analysis.

Derived from cmig_check_batch_errors.py, which it may replace. The job list comes
from scriptlist.txt rather than from the .err files present, so a job that never
started is reported rather than missed.

Memory, CPU and disk I/O per job come from the scheduler's accounting
(cmig_cluster_utils.job_usage), with the job ID and host name read from the
header the job wrapper prints at the top of each .out.

Given --raw_rootdir and the output roots, each session of the batch is also
checked for steps left stale (cmig_mkjobs_bids_pipeline.check_session), as a
job that runs only some steps leaves the steps after them; these are listed
as warnings.

Created:  09/08/2026 by Don Hagler
Prev mod: 09/11/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Last mod: 11/10/2026 by Don Hagler
"""

import argparse
import collections
import os
import re
import statistics
import sys
import time
from datetime import datetime

import cmig_cluster_utils as clu

# plots are optional: the tool is useful without them, and matplotlib is not
# present everywhere MMPS is. same guard the other MMPS python scripts use
try:
    import matplotlib
    matplotlib.use('Agg')  # non-interactive backend
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


N_TAIL = 10

# the wrappers qcshjob and qmatjob echo a bare 'finished' with no period, so
# the MATLAB 'finished\.' pattern never matches; match the bare word
FINISHED_RE = re.compile(r'^finished\s*$')

STATE_FINISHED = 'finished'
STATE_STARTED = 'started'
STATE_NEVER_STARTED = 'never_started'

# cmig_write_status prints these four lines, in this order, and they are read
# by fixed offset from the outdir anchor, so nothing may be inserted between
# them. anything added to the block prints ABOVE outdir instead, which is
# where module and startdate go
STATUS_PREFIX = 'cmig_write_status: '
STATUS_FIELDS = ('outdir', 'procdate', 'errcode', 'errmsg')
# the lines above the anchor, which a block may print in any order and may
# omit; they are scanned for rather than read at fixed offsets
STATUS_HEADER_FIELDS = ('startdate', 'module')
PROCDATE_FMT = '%Y-%m-%dT%H:%M:%S.%f'

# a module that succeeds follows its block with these two lines; one that
# fails prints neither, and names itself only in the .err
MODULE_DONE_RE = re.compile(r'^(\w+): finished\s*$')

# the pipeline banners every module it runs, with the task where there is one,
# so a banner names a block better than anything else available
BANNER_RE = re.compile(r'^([A-Za-z_]\w*(?: \([^)]*\))?)$')
BANNER_RULE_RE = re.compile(r'^=+$')
ERR_MODULE_RE = re.compile(r'^(\w+)(?:\((\w+)\))?:')

UNATTRIBUTED = '(unattributed)'

# mmil_check_batchdir truncates an errmsg at this length in its tsv
MAX_ERR_LENGTH = 3000
TSV_COLUMNS = ('jobname', 'state', 'errored', 'category', 'start', 'finish',
               'runtime_hours', 'cluster', 'jobid', 'mem_gb', 'mem_type',
               'cores', 'cpu_hours', 'read_gb', 'write_gb', 'errmsg')

# qcshjob prints these at the top of the .out, on slurm and SGE alike
OUT_HEADER_RE = re.compile(r'^(HOSTNAME|JOB_ID|JOB_NAME) =\s*(\S*)')
N_HEADER_LINES = 20

# memory as each scheduler's accounting reports it
MEM_LABELS = {'rss': 'resident (MaxRSS)', 'vmem': 'virtual (maxvmem)'}

FIG_SIZE = (4, 4)
FIG_DPI = 300
N_BINS = 20

BASIS_MARKER = 'marker'
BASIS_MTIME = 'mtime'

# a rename lands the marker's ctime a fraction of a second after the job's last
# write to the .out. a wider gap means something other than the mv moved the
# ctime, a copy being the usual cause, so the .out and .err mtimes are used
CTIME_TOLERANCE = 60.0

STATE_LABELS = {
    STATE_FINISHED: 'finished',
    STATE_STARTED: 'started but not finished',
    STATE_NEVER_STARTED: 'never started',
}

# stderr output that does not indicate job failure; a .err holding only these
# is treated as empty, so the job is judged from the .out like any other
BENIGN_ERR_PATTERNS = [
    r'do_ypcall',                          # NIS lookup noise from the batch shell
    r'clnt_call: RPC:',
    r'Warning: no access to tty',
    r'Thus no job control in this shell',
    r'^Fontconfig warning',                # graphics library config parse warning
]

CATEGORY_ORDER = [
    'success',
    'success_benign',
    'started',
    'never_started',
    'status_missing',
    'matlab_logical',
    'mat_missing',
    'fsurf_error',
    'mgh_read',
    'input_error',
    'proc_error',
    'slurm_cancel',
    'unknown',
]

CATEGORY_LABELS = {
    'success':        'Successful completion',
    'success_benign': 'Successful completion (benign .err output)',
    'started':        'Started but not finished',
    'never_started':  'Never started (no .out or .err)',
    'status_missing': 'Status file not found (pipeline step failed)',
    'matlab_logical': 'MATLAB logical conditional error (select_scan)',
    'mat_missing':    'Missing .mat file (atlas registration failure)',
    'fsurf_error':    'FreeSurfer recon-all failure',
    'mgh_read':       'MGH file read error (corrupted/incomplete file)',
    'input_error':    'No usable input data',
    'proc_error':     'Processing module error (errcode=1)',
    'slurm_cancel':   'Job cancelled by SLURM scheduler',
    'unknown':        'Unknown error',
}


def strip_benign(err_text):
    """Drop stderr lines that do not indicate job failure."""
    keep = []
    for line in err_text.splitlines():
        if not line.strip():
            continue
        if any(re.search(pat, line) for pat in BENIGN_ERR_PATTERNS):
            continue
        keep.append(line)
    return '\n'.join(keep)


def categorize(err_text, state, blocks=None):
    # a job that never ran has no .err to read
    if state == STATE_NEVER_STARTED:
        return 'never_started'
    # a .err holding nothing but benign output is not evidence of failure, so
    # the job is judged by its blocks and then its state; success_benign is
    # kept separate from success so the noise stays visible in the summary
    if not strip_benign(err_text).strip():
        # a module can report its failure on stdout alone, as the input checks
        #   in cmig_fsurf_bids.py do, so a silent .err is not proof of success
        #   while a block carries a non-zero errcode
        failed = [b for b in (blocks or []) if b.get('errcode', '0') != '0']
        if failed:
            errmsg = failed[0].get('errmsg') or ''
            if ('input not found' in errmsg or 'no usable' in errmsg
                    or 'no data found' in errmsg):
                return 'input_error'
            return 'proc_error'
        if state == STATE_STARTED:
            return 'started'
        if err_text.strip():
            return 'success_benign'
        return 'success'
    # a job with a real .err is reported by its error whatever its state, so a
    # job the scheduler killed reads as cancelled rather than as started
    if 'CANCELLED' in err_text:
        return 'slurm_cancel'
    if 'fs_load_mgh' in err_text or ('tried to read' in err_text and 'actually read' in err_text):
        return 'mgh_read'
    if 'nonLogicalConditional' in err_text or 'Operands to the logical' in err_text:
        return 'matlab_logical'
    if 'couldNotReadFile' in err_text or ('Unable to find file' in err_text and '.mat' in err_text):
        return 'mat_missing'
    if 'FreeSurfer command ended in error' in err_text:
        return 'fsurf_error'
    if 'errcode=1' in err_text:
        # 'no usable' from cmig_proc_bids_anat, 'no data found' from cmig_rawqc_bids
        if 'no usable' in err_text or 'no data found' in err_text:
            return 'input_error'
        return 'proc_error'
    if 'status file not found' in err_text:
        return 'status_missing'
    return 'unknown'


def extract_context(out_text, category):
    ctx = {}
    if not out_text:
        return ctx

    # proc outdir
    m = re.search(r'cmig_proc_bids_anat: outdir = (\S+)', out_text)
    if m:
        ctx['proc_outdir'] = m.group(1)

    # FreeSurfer subject dir (from "NOTE: complete FreeSurfer recon dir:" line)
    m = re.search(r'NOTE: complete FreeSurfer recon dir:\n(\S+)', out_text)
    if m:
        ctx['fsurf_dir'] = m.group(1)

    # FreeSurfer subject dir (from standalone path line after SUBJECTS_DIR)
    if 'fsurf_dir' not in ctx:
        m = re.search(r'\n(/\S+/freesurfer/sub-[^\s/]+)\s*\n', out_text)
        if m:
            ctx['fsurf_dir'] = m.group(1)

    # Truncated NIFTI path (FreeSurfer mri_convert failure)
    m = re.search(r'error: niiRead\(\): error reading from (\S+)', out_text)
    if m:
        ctx['bad_nifti'] = m.group(1)

    # Analysis ERROR-level missing file
    m = re.search(r'cmig_analyze_bids_anat: ERROR: missing file (\S+)', out_text)
    if m:
        ctx['analysis_error'] = m.group(1)

    # Processing module errmsg (from .out)
    m = re.search(r'errmsg = (.+)', out_text)
    if m and m.group(1).strip():
        ctx['proc_errmsg'] = m.group(1).strip()

    # Last N non-empty lines of .out
    all_lines = out_text.splitlines()
    tail = [l for l in all_lines if l.strip()][-N_TAIL:]
    if tail:
        ctx['out_tail'] = tail

    return ctx


def err_excerpt(err_text, category):
    """Return .err lines relevant to the root cause."""
    lines = err_text.strip().splitlines()
    # For categories where "status file not found" is a consequence, strip it
    if category != 'status_missing':
        lines = [l for l in lines if 'status file not found' not in l]
    while lines and not lines[-1].strip():
        lines.pop()
    return '\n'.join(lines)


def err_signature(excerpt):
    """MATLAB error identifier and innermost stack frame from a .err excerpt.

    One category can hold unrelated failures: 'status file not found' is written
    by the python wrapper for any MATLAB error in a pipeline step, so a single
    22-job category held three separate bugs. Grouping on the identifier and the
    top 'Error in' frame separates them without a hardcoded category per bug.
    """
    ident = ''
    match = re.search(r'^([A-Za-z]\w*(?::\w+)+)\s*$', excerpt, re.MULTILINE)
    if match:
        ident = match.group(1)
    frame = ''
    match = re.search(r'^Error in (\S+)', excerpt, re.MULTILINE)
    if match:
        frame = match.group(1)
    sig = ' '.join(p for p in (ident, frame) if p)
    return sig or 'no MATLAB error in .err'


def group_by_signature(jobs):
    """Group (job_name, excerpt, ctx) tuples by err_signature, first seen first."""
    groups = []
    index = {}
    for job in jobs:
        sig = err_signature(job[1])
        if sig not in index:
            index[sig] = len(groups)
            groups.append((sig, []))
        groups[index[sig]][1].append(job)
    return groups


def format_context(ctx):
    lines = []
    if 'proc_outdir' in ctx:
        lines.append(f"  proc outdir:     {ctx['proc_outdir']}")
    if 'fsurf_dir' in ctx:
        lines.append(f"  fsurf dir:       {ctx['fsurf_dir']}")
    if 'bad_nifti' in ctx:
        lines.append(f"  bad NIFTI:       {ctx['bad_nifti']}")
    if 'analysis_error' in ctx:
        lines.append(f"  analysis error:  {ctx['analysis_error']}")
    if 'proc_errmsg' in ctx:
        lines.append(f"  errmsg:          {ctx['proc_errmsg']}")
    if 'out_tail' in ctx:
        lines.append(f"  last {N_TAIL} lines of .out:")
        for line in ctx['out_tail']:
            lines.append(f"    {line}")
    return '\n'.join(lines)


def read_file(path):
    try:
        with open(path) as f:
            return f.read().replace('\x00', '')
    except OSError:
        return ''


def read_scriptlist(batchdir):
    """Job names the batch was supposed to run, from the master scriptlist.

    A multi-cluster run also leaves scriptlist_batch<N>_<clust>.txt subsets
    beside the master, but the master holds every job, so only it is read.
    """
    path = os.path.join(batchdir, 'scriptlist.txt')
    if not os.path.isfile(path):
        print(f'WARNING: scriptlist not found: {path}', file=sys.stderr)
        return []
    jobs = []
    for line in read_file(path).splitlines():
        name = line.strip()
        if name:
            jobs.append(name)
    return jobs


def resolve_job_files(pbsout_dir, job_name, entries):
    """The .out and .err paths for one job, None where either is absent.

    Some jobs write {job_name}_{date_time}.out instead of {job_name}.out; the
    most recent of those wins and the .err is re-derived from the same stem so
    the pair stays matched.
    """
    out_path = os.path.join(pbsout_dir, job_name + '.out')
    err_path = os.path.join(pbsout_dir, job_name + '.err')
    if not os.path.exists(out_path):
        suffixed = [e for e in entries
                    if e.startswith(job_name + '_') and e.endswith('.out')]
        if suffixed:
            newest = max(suffixed, key=lambda e: os.path.getmtime(
                os.path.join(pbsout_dir, e)))
            out_path = os.path.join(pbsout_dir, newest)
            err_path = out_path[:-4] + '.err'
    return (out_path if os.path.exists(out_path) else None,
            err_path if os.path.exists(err_path) else None)


def out_header(out_path):
    """{'HOSTNAME': ..., 'JOB_ID': ..., 'JOB_NAME': ...} from the top of
    the .out, empty where the .out lacks the header."""
    header = {}
    if out_path is None:
        return header
    try:
        with open(out_path, errors='replace') as f:
            for _ in range(N_HEADER_LINES):
                line = f.readline()
                if not line:
                    break
                match = OUT_HEADER_RE.match(line)
                if match and match.group(2):
                    header[match.group(1)] = match.group(2)
    except OSError:
        pass
    return header


def markers_present(pbsout_dir):
    """Whether the batchdir carries any marker files at all.

    A batchdir with .out files but no markers was run outside qcshjob and
    qmatjob, and then the last line of the .out is the only finished signal.
    """
    if not os.path.isdir(pbsout_dir):
        return False
    for entry in os.listdir(pbsout_dir):
        if entry.startswith('.finished_') or entry.startswith('.isrunning_'):
            return True
    return False


def out_finished(out_path):
    """Whether the .out ends with the wrappers' bare 'finished' line."""
    if out_path is None:
        return False
    lines = [l for l in read_file(out_path).splitlines() if l.strip()]
    return bool(lines) and bool(FINISHED_RE.match(lines[-1]))


def job_state(pbsout_dir, job_name, out_path, err_path, has_markers):
    """finished, started but not finished, or never started.

    Marker existence is the signal: the wrappers touch .isrunning_ when the
    job starts and rename it to .finished_ only after the job body returns.
    A job killed by the scheduler never reaches the rename, so its stale
    .isrunning_ reads the same as a live one.
    """
    if has_markers:
        if os.path.exists(os.path.join(pbsout_dir, '.finished_' + job_name)):
            return STATE_FINISHED
        if os.path.exists(os.path.join(pbsout_dir, '.isrunning_' + job_name)):
            return STATE_STARTED
    if out_path is None and err_path is None:
        return STATE_NEVER_STARTED
    # files but no marker: judge from the .out, as for a marker-less batchdir
    return STATE_FINISHED if out_finished(out_path) else STATE_STARTED


def job_errored(err_text):
    """Whether the .err holds anything that indicates failure.

    A .err carrying only benign stderr is not evidence of failure, so it is
    stripped before the emptiness test. State and error are independent: the
    wrappers finish a job whatever its body did, so failure comes from here.
    """
    return bool(strip_benign(err_text).strip())


def job_times(pbsout_dir, job_name, state, out_path, err_path):
    """Start and finish for one job as epoch seconds, with the basis used.

    The wrappers touch .isrunning_ when the job starts and rename it to
    .finished_ when it ends, and rename preserves mtime, so a .finished_
    marker carries the start in its mtime and the finish in its ctime. A
    started job has no finish yet, so its elapsed time runs to now.

    Falls back to the .err and .out mtimes where there is no marker, or where
    the marker's ctime has been moved by something other than the rename. That
    fallback collapses to zero for a job that wrote to stderr as it died,
    since the .err mtime is the last stderr write rather than the start.
    """
    if state == STATE_NEVER_STARTED:
        return None, None, None

    if state == STATE_STARTED:
        run_marker = os.path.join(pbsout_dir, '.isrunning_' + job_name)
        if os.path.exists(run_marker):
            return os.path.getmtime(run_marker), time.time(), BASIS_MARKER
        # no marker: the .err is created at job start and stays untouched
        # unless the job writes to stderr
        if err_path:
            return os.path.getmtime(err_path), time.time(), BASIS_MTIME
        return None, None, None

    fin_marker = os.path.join(pbsout_dir, '.finished_' + job_name)
    if os.path.exists(fin_marker):
        st = os.stat(fin_marker)
        if st.st_ctime >= st.st_mtime:
            if out_path is None:
                return st.st_mtime, st.st_ctime, BASIS_MARKER
            gap = abs(st.st_ctime - os.path.getmtime(out_path))
            if gap <= CTIME_TOLERANCE:
                return st.st_mtime, st.st_ctime, BASIS_MARKER

    if out_path and err_path:
        return (os.path.getmtime(err_path), os.path.getmtime(out_path),
                BASIS_MTIME)
    return None, None, None


def summarize(values):
    """Count, min, max, median and mean of a list of run times, None if empty."""
    if not values:
        return None
    return {
        'n': len(values),
        'min': min(values),
        'max': max(values),
        'median': statistics.median(values),
        'mean': statistics.mean(values),
    }


def format_runtimes(records):
    """The run-time statistics for finished and for started jobs, in hours.

    Errored jobs are included: the marker basis measures them correctly, and
    it was only the .out and .err mtime basis collapsing on them that made
    mmil_check_batchdir leave them out.
    """
    lines = []
    for state in (STATE_FINISHED, STATE_STARTED):
        hours = [(r['finish'] - r['start']) / 3600.0 for r in records
                 if r['state'] == state and r['start'] is not None]
        stats = summarize(hours)
        if stats is None:
            continue
        lines.append('')
        lines.append('run-time (hours) for %s jobs:' % STATE_LABELS[state])
        lines.append('  # timed: %d' % stats['n'])
        for key in ('min', 'max', 'median', 'mean'):
            lines.append('  %s: %.4f' % (key, stats[key]))
    fallback = sum(1 for r in records if r['basis'] == BASIS_MTIME)
    if fallback:
        lines.append('')
        lines.append('  NOTE: %d job(s) timed from the .out and .err mtimes '
                     'rather than the marker' % fallback)
    return '\n'.join(lines)


def resolve_batchdir(batchdir, root_batchdir):
    """A full path, or a bare batch name resolved against the root.

    qjobs resolves a batch name the same way, so either form works here.
    """
    if os.path.isdir(batchdir):
        return batchdir
    candidate = os.path.join(os.path.expanduser(root_batchdir), batchdir)
    if os.path.isdir(candidate):
        return candidate
    return batchdir


def tsv_errmsg(err_text):
    """The .err as one tsv field: benign lines dropped, whitespace escaped."""
    text = strip_benign(err_text).strip()
    if not text:
        return ''
    for char, name in (('\t', '<TAB> '), ('\n', '<NL> '), ('\r', '<CR> ')):
        text = text.replace(char, name)
    text = text.replace('{', '').replace('}', '')
    return text[:MAX_ERR_LENGTH]


def stamp_time(epoch):
    """A procdate-shaped timestamp, or empty where there is none."""
    if epoch is None:
        return ''
    return time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(epoch))


def fmt_value(value, fmt):
    """A number for the tsv, or empty where there is none."""
    return '' if value is None else fmt % value


def write_jobs_tsv(path, records):
    """One row per scriptlist job, the machine-readable form of the accounting.

    This supersedes mmil_check_batchdir's errors.tsv, which held only the jobs
    with a non-empty .err; filtering on the errored column reproduces it.
    """
    with open(path, 'w') as f:
        f.write('\t'.join(TSV_COLUMNS) + '\n')
        for record in records:
            runtime = ''
            if record['start'] is not None and record['finish'] is not None:
                runtime = '%.4f' % ((record['finish'] - record['start'])
                                    / 3600.0)
            usage = record.get('usage') or {}
            f.write('\t'.join((
                record['name'],
                record['state'],
                '1' if record['errored'] else '0',
                record.get('category', ''),
                stamp_time(record['start']),
                stamp_time(record['finish']),
                runtime,
                record.get('cluster', ''),
                record.get('jobid', ''),
                fmt_value(usage.get('mem_gb'), '%.3f'),
                usage.get('mem_type') or '',
                fmt_value(usage.get('cores'), '%.3f'),
                fmt_value(usage.get('cpu_h'), '%.4f'),
                fmt_value(usage.get('read_gb'), '%.3f'),
                fmt_value(usage.get('write_gb'), '%.3f'),
                record.get('errmsg', ''),
            )) + '\n')


def write_summary(path, text):
    """The counts, run times and module table, as printed."""
    with open(path, 'w') as f:
        f.write(text + '\n')


def save_fig(fig, path):
    """Write one figure and release it."""
    fig.savefig(path, dpi=FIG_DPI, format='png', bbox_inches='tight')
    plt.close(fig)


def barh_plot(path, values, names, xlabel, title, xmax=None):
    """A horizontal bar chart, in place of the MATLAB's mmil_barh."""
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    positions = list(range(len(values)))
    ax.barh(positions, values)
    ax.set_yticks(positions)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    if xmax:
        ax.set_xlim(0, xmax)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    save_fig(fig, path)


def histo_plot(path, values, title):
    """The distribution of run times, which the four summary numbers cannot show."""
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.hist(values, bins=N_BINS)
    ax.set_xlabel('run time (hours)')
    ax.set_ylabel('# of jobs')
    ax.set_title(title)
    save_fig(fig, path)


def scatter_plot(path, starts, times, title):
    """Run time against start time, which shows contention across the batch."""
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.plot(starts, times, 'bo')
    ax.set_xlabel('start time (hours)')
    ax.set_ylabel('run time (hours)')
    ax.set_title(title)
    save_fig(fig, path)


def make_plots(records, outdir, batchname):
    """The seven plots mmil_check_batchdir makes, as PNGs.

    Skipped where matplotlib is missing. The started-job plots carry _started
    in their names rather than the MATLAB's _running, matching the category
    rename.
    """
    if not MATPLOTLIB_AVAILABLE:
        print('WARNING: matplotlib not available, skipping plots',
              file=sys.stderr)
        return
    counts = {state: 0 for state in STATE_LABELS}
    for record in records:
        counts[record['state']] += 1
    n_jobs = len(records)
    if not n_jobs:
        return
    title = batchname.replace('_', ' ')

    barh_plot(os.path.join(outdir, batchname + '_summary.png'),
              [n_jobs, counts[STATE_FINISHED], counts[STATE_STARTED],
               counts[STATE_NEVER_STARTED],
               sum(1 for r in records if r['errored'])],
              ['total', 'finished', 'started', 'never started', 'errors'],
              '# of jobs', 'summary of %s jobs' % title, xmax=n_jobs)

    # the earliest start anchors the x axis, as mmil_check_batchdir does
    starts = [r['start'] for r in records if r['start'] is not None]
    origin = min(starts) if starts else 0.0

    for state, tag in ((STATE_FINISHED, 'finished'), (STATE_STARTED, 'started')):
        rows = [r for r in records
                if r['state'] == state and r['start'] is not None]
        if not rows:
            continue
        hours = [(r['finish'] - r['start']) / 3600.0 for r in rows]
        offsets = [(r['start'] - origin) / 3600.0 for r in rows]
        label = '%s (%d of %d)' % (STATE_LABELS[state], len(rows), n_jobs)
        histo_plot(os.path.join(outdir, '%s_histo_%s.png' % (batchname, tag)),
                   hours, label)
        scatter_plot(os.path.join(outdir, '%s_time_%s.png' % (batchname, tag)),
                     offsets, hours, label)
        stats = summarize(hours)
        barh_plot(os.path.join(outdir, '%s_summary_%s.png' % (batchname, tag)),
                  [stats['min'], stats['max'], stats['median'], stats['mean']],
                  ['min', 'max', 'median', 'mean'],
                  'run time (hours)', label, xmax=stats['max'])


def usage_histo_plot(path, rows, key, xlabel, title):
    """One usage measure for the finished jobs, stacked by cluster; None
    values (read on SGE) are left out."""
    data, labels = [], []
    for cluster, usage in rows.items():
        values = [u[key] for u in usage if u[key] is not None]
        if values:
            data.append(values)
            label = cluster
            if key == 'mem_gb':
                types = sorted(set(u['mem_type'] for u in usage))
                label = '%s (%s)' % (cluster, ', '.join(types))
            labels.append(label)
    if not data:
        return
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.hist(data, bins=N_BINS, stacked=True, label=labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('# of jobs')
    ax.set_title(title)
    ax.legend(fontsize=7)
    save_fig(fig, path)


def make_usage_plots(records, outdir, batchname):
    """Histograms of memory, cores, write GB and read GB per finished job."""
    if not MATPLOTLIB_AVAILABLE:
        return
    rows = usage_rows(records)
    if not rows:
        return
    stem = os.path.join(outdir, batchname)
    usage_histo_plot(stem + '_usage_mem.png', rows, 'mem_gb', 'memory (GB)',
                     'peak memory (rss: slurm, vmem: SGE)')
    usage_histo_plot(stem + '_usage_cores.png', rows, 'cores',
                     'cores (CPU time / elapsed)', 'cores per job')
    usage_histo_plot(stem + '_usage_write.png', rows, 'write_gb',
                     'written (GB)', 'disk write per job')
    usage_histo_plot(stem + '_usage_read.png', rows, 'read_gb', 'read (GB)',
                     'disk read per job (slurm only)')


def module_colors(names):
    """A stable module to colour mapping, so every module plot agrees."""
    cmap = plt.get_cmap('tab20')
    return {name: cmap(i % 20) for i, name in enumerate(names)}


def module_box_plot(path, stats, colors):
    """Per-module run times as boxes on a log axis.

    Durations span three orders of magnitude, so a linear axis would flatten
    everything but the dominant step. This shows the spread the table's four
    numbers cannot.
    """
    names, data = [], []
    for name, entry in stats.items():
        hours = [t / 3600.0 for t in entry['times'] if t > 0]
        if hours:
            names.append(name)
            data.append(hours)
    if not data:
        return
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(names) + 2))
    boxes = ax.boxplot(data, labels=names, vert=False, patch_artist=True)
    for patch, name in zip(boxes['boxes'], names):
        patch.set_facecolor(colors[name])
    ax.set_xscale('log')
    ax.set_xlabel('run time (hours, log scale)')
    ax.set_title('run time by module')
    save_fig(fig, path)


def module_timeline_plot(path, job_segments, origin, colors):
    """One row per job, modules as segments, sorted by total duration.

    The only view that shows per-job structure: which module the slow jobs are
    slow in, and where a failure cut a job short.
    """
    rows = []
    for name, segments in job_segments:
        timed = [g for g in segments
                 if g['start'] is not None and g['finish'] is not None]
        if timed:
            rows.append((timed[-1]['finish'] - timed[0]['start'], timed))
    if not rows:
        return
    rows.sort(key=lambda r: r[0])
    height = min(0.12 * len(rows) + 1.5, 24)
    fig, ax = plt.subplots(figsize=(9, height))
    for y, (_, timed) in enumerate(rows):
        for segment in timed:
            start = (segment['start'] - origin) / 3600.0
            width = (segment['finish'] - segment['start']) / 3600.0
            ax.broken_barh([(start, width)], (y - 0.4, 0.8),
                           facecolors=colors[segment['module']])
    ax.set_xlabel('hours since the first job started')
    ax.set_ylabel('jobs, sorted by duration')
    ax.set_ylim(-1, len(rows))
    ax.set_title('module timeline per job')
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[n]) for n in colors]
    # outside the axes, so it does not cover the longest jobs
    ax.legend(handles, list(colors), loc='lower center',
              bbox_to_anchor=(0.5, 1.02), ncol=2, fontsize=7)
    save_fig(fig, path)


def module_stack_plot(path, stats, n_jobs, colors):
    """What a typical job is made of: median duration times runs per job."""
    left = 0.0
    fig, ax = plt.subplots(figsize=(8, 2.5))
    for name, entry in stats.items():
        hours = [t / 3600.0 for t in entry['times']]
        if not hours or not n_jobs:
            continue
        share = statistics.median(hours) * len(hours) / float(n_jobs)
        ax.barh([0], [share], left=left, color=colors[name], label=name)
        left += share
    if left == 0.0:
        plt.close(fig)
        return
    ax.set_yticks([])
    ax.set_xlabel('hours')
    ax.set_title('composition of a typical job')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol=2,
              fontsize=7)
    save_fig(fig, path)


def module_time_plot(path, job_segments, origin, colors):
    """Module duration against start time, which attributes a batch-wide
    slowdown to a particular step."""
    fig, ax = plt.subplots(figsize=(7, 5))
    seen = set()
    for _, segments in job_segments:
        for segment in segments:
            if segment['start'] is None or segment['finish'] is None:
                continue
            name = segment['module']
            ax.plot((segment['start'] - origin) / 3600.0,
                    (segment['finish'] - segment['start']) / 3600.0,
                    'o', markersize=3, color=colors[name],
                    label=name if name not in seen else '')
            seen.add(name)
    if not seen:
        plt.close(fig)
        return
    ax.set_xlabel('start time (hours)')
    ax.set_ylabel('run time (hours)')
    ax.set_title('module run time against start time')
    ax.legend(fontsize=7, markerscale=2)
    save_fig(fig, path)


def module_count_plot(path, stats, colors, key, xlabel, title):
    """Runs or errors per module, as a horizontal bar."""
    names = list(stats)
    values = [stats[n]['runs'] if key == 'runs' else stats[n]['errors']
              for n in names]
    if not names or not any(values):
        return
    fig, ax = plt.subplots(figsize=(8, 0.4 * len(names) + 1.5))
    positions = list(range(len(names)))
    ax.barh(positions, values, color=[colors[n] for n in names])
    ax.set_yticks(positions)
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    save_fig(fig, path)


def module_histo_plot(path, stats, colors):
    """A run-time histogram per module, on a shared figure."""
    panels = [(n, [t / 3600.0 for t in e['times']])
              for n, e in stats.items() if e['times']]
    if not panels:
        return
    fig, axes = plt.subplots(len(panels), 1, figsize=(6, 2 * len(panels)),
                             squeeze=False)
    for ax, (name, hours) in zip(axes[:, 0], panels):
        ax.hist(hours, bins=N_BINS, color=colors[name])
        ax.set_ylabel('# runs')
        ax.set_title(name, fontsize=9)
    axes[-1, 0].set_xlabel('run time (hours)')
    fig.tight_layout()
    save_fig(fig, path)


def make_module_plots(job_segments, stats, records, outdir, batchname):
    """The plots the per-module timing makes possible, as PNGs."""
    if not MATPLOTLIB_AVAILABLE or not stats:
        return
    starts = [r['start'] for r in records if r['start'] is not None]
    origin = min(starts) if starts else 0.0
    colors = module_colors(list(stats))
    stem = os.path.join(outdir, batchname)
    module_box_plot(stem + '_module_box.png', stats, colors)
    module_timeline_plot(stem + '_module_timeline.png', job_segments, origin,
                         colors)
    module_stack_plot(stem + '_module_stack.png', stats, len(records), colors)
    module_time_plot(stem + '_module_time.png', job_segments, origin, colors)
    module_count_plot(stem + '_module_runs.png', stats, colors, 'runs',
                      '# of runs', 'runs per module')
    module_count_plot(stem + '_module_errors.png', stats, colors, 'errors',
                      '# of errors', 'errors per module')
    module_histo_plot(stem + '_module_histo.png', stats, colors)


def parse_procdate(text):
    """Epoch seconds from a procdate stamp, None if it will not parse.

    startdate is passed here too and a block may omit it, so None is a normal
    input rather than a malformed stamp.
    """
    try:
        return datetime.strptime(text, PROCDATE_FMT).timestamp()
    except (TypeError, ValueError):
        return None


def module_from_pair(lines, index):
    """The module named by the finished pair after a status block, or None."""
    if index >= len(lines):
        return None
    match = MODULE_DONE_RE.match(lines[index])
    if not match:
        return None
    return match.group(1)


def parse_status_blocks(out_text):
    """The cmig_write_status blocks in the .out, in the order they were written.

    procdate is stamped as the block is written, so it is the moment the module
    finished; startdate, where the module supplied one, is its entry time, so a
    duration needs no differencing against the module before it.

    A module that succeeded is named by the finished pair that follows its
    block, as before. One that failed prints no pair, and is now named by the
    block's own module line where there is one, so the .err has to supply the
    name only for a .out written before cmig_write_status printed it. The two
    give different forms, cmig_proc_bids_anat against desc-proc, so the pair
    stays primary and the module line is kept separately as status_module.
    """
    blocks = []
    lines = out_text.splitlines()
    for i, line in enumerate(lines):
        if not line.startswith(STATUS_PREFIX + 'outdir = '):
            continue
        fields = {}
        for offset, key in enumerate(STATUS_FIELDS):
            prefix = STATUS_PREFIX + key + ' = '
            if i + offset < len(lines) and lines[i + offset].startswith(prefix):
                fields[key] = lines[i + offset].split(' = ', 1)[1]
        if 'procdate' not in fields:
            continue
        # the lines above the anchor, scanned backward rather than read at
        #   fixed offsets, since any of them may be absent: a module that
        #   passes no startdate leaves module directly above outdir
        for j in range(i - 1, max(i - 1 - len(STATUS_HEADER_FIELDS), -1), -1):
            for key in STATUS_HEADER_FIELDS:
                prefix = STATUS_PREFIX + key + ' = '
                if lines[j].startswith(prefix):
                    fields[key] = lines[j].split(' = ', 1)[1]
                    break
            else:
                break
        fields['line'] = i
        # the finished pair stays the primary name, since it gives the module
        #   function name where the block's own module line gives the status
        #   entity, desc-proc against cmig_proc_bids_anat; the module line is
        #   the fallback, which is what a failed block now has where it had
        #   nothing. both are kept, so a consumer can take either
        fields['status_module'] = fields.get('module')
        fields['module'] = module_from_pair(lines, i + len(STATUS_FIELDS)) \
            or fields.get('module')
        blocks.append(fields)
    return blocks


def parse_banners(out_text):
    """The module banners in the .out, as (line index, module name).

    A banner is a module name fenced by rules of = characters. It marks a
    module starting, which the status blocks do not.
    """
    banners = []
    lines = out_text.splitlines()
    for i in range(1, len(lines) - 1):
        if (BANNER_RULE_RE.match(lines[i - 1])
                and BANNER_RULE_RE.match(lines[i + 1])
                and BANNER_RE.match(lines[i])):
            banners.append((i, lines[i]))
    return banners


def err_module(err_text):
    """The module the .err reports as failing, with its task where named.

    A module that fails prints no finished pair, so the .err is the only place
    it names itself, as '<module>: errcode=1' and '<module>(<task>): <message>'.
    """
    name = None
    task = None
    for line in strip_benign(err_text).splitlines():
        match = ERR_MODULE_RE.match(line)
        if not match:
            continue
        if name is None:
            name = match.group(1)
        if match.group(2) and match.group(1) == name:
            task = match.group(2)
    if name and task:
        return '%s(%s)' % (name, task)
    return name


def offset_from(when, job_start):
    """Seconds from the job start, None where either stamp is missing.

    A module that was skipped has no duration but still has a place in the
    job, which is the boundary it was reached at; order_modules needs that.
    """
    if when is None or job_start is None:
        return None
    return when - job_start


def job_module_segments(blocks, banners, job_start, err_text):
    """One timed segment per status block: module, start, finish, errored.

    A segment runs from the previous block's procdate, or from the job start
    for the first, since each procdate is a completion. Where a banner inside a
    segment names a module other than the block's own, that module wrote no
    status block of its own and gets an untimed segment, so every module is
    named once and none is folded into another.
    """
    segments = []
    failed_name = None
    previous = job_start
    prev_line = -1
    for block in blocks:
        finish = parse_procdate(block['procdate'])
        errored = block.get('errcode', '0') != '0'
        name = block['module']
        if name is None and errored:
            if failed_name is None:
                failed_name = err_module(err_text)
            name = failed_name
        # the banner is the best name available: it covers a module that
        # succeeded and one that failed alike, and it carries the task. the
        # last banner in a segment is the one this block belongs to
        in_segment = [n for (i, n) in banners if prev_line < i < block['line']]
        if in_segment:
            # every banner but the last names a module that wrote no block of
            #   its own, typically one cmig_write_status skipped because its
            #   status file was already there. it gets its own segment, with no
            #   times to give it, rather than having its name folded into the
            #   module that follows
            for earlier in in_segment[:-1]:
                segments.append({'module': earlier, 'start': None,
                                 'finish': None, 'errored': False,
                                 'offset': offset_from(previous, job_start)})
            name = in_segment[-1]
        if name is None:
            name = UNATTRIBUTED
        # startdate is the module's own entry time, so whatever ran between the
        #   previous block and this one stays out of this module's duration
        start = parse_procdate(block.get('startdate'))
        if start is None:
            start = previous
        segments.append({'module': name, 'start': start, 'finish': finish,
                         'errored': errored,
                         'offset': offset_from(start, job_start)})
        if finish is not None:
            previous = finish
        prev_line = block['line']
    return segments


def add_module_times(stats, segments):
    """Fold one job's segments into the per-module statistics.

    Every segment's errcode counts: a failed module does not end the job, so
    the last block can report success for a job that contains a failure.

    runs counts every segment and times only the ones that carry both stamps,
    so a module that was skipped is still reported as having run. position is
    where the segment sits in its own job, which order_modules reduces to the
    running order.
    """
    for position, segment in enumerate(segments):
        entry = stats.setdefault(segment['module'],
                                 {'times': [], 'errors': 0, 'runs': 0,
                                  'positions': [], 'offsets': []})
        entry['runs'] += 1
        entry['positions'].append(position)
        if segment.get('offset') is not None:
            entry['offsets'].append(segment['offset'])
        if segment['errored']:
            entry['errors'] += 1
        if segment['start'] is not None and segment['finish'] is not None:
            entry['times'].append(segment['finish'] - segment['start'])


def collect_module_aliases(blocks, aliases):
    """Learn the status entity to module name map from the blocks that carry
    both, so a block that carries only the entity can be named the same way.

    A module that succeeded is named by the finished pair that follows its
    block, cmig_fsurf_bids, where the block's own module line gives the status
    entity, desc-fsurf. One that failed prints no pair and has only the entity,
    which would report the same module under two names.
    """
    for block in blocks:
        entity = block.get('status_module')
        name = block.get('module')
        if entity and name and entity != name:
            aliases.setdefault(entity, set()).add(name)


def resolve_module_aliases(stats, job_segments, aliases):
    """Rename the entity-named modules, merging their statistics.

    Only an entity that names one module across the whole batch is renamed:
    desc-proc is written by cmig_proc_bids_anat and cmig_proc_bids_func alike,
    and a job carrying both leaves it ambiguous.
    """
    resolved = {entity: names.pop() for entity, names in aliases.items()
                if len(names) == 1}
    if not resolved:
        return stats
    for _, segments in job_segments:
        for segment in segments:
            segment['module'] = resolved.get(segment['module'],
                                             segment['module'])
    merged = collections.OrderedDict()
    for name, entry in stats.items():
        key = resolved.get(name, name)
        if key not in merged:
            merged[key] = entry
            continue
        target = merged[key]
        target['times'].extend(entry['times'])
        target['positions'].extend(entry['positions'])
        target['offsets'].extend(entry['offsets'])
        target['errors'] += entry['errors']
        target['runs'] += entry['runs']
    return merged


def order_modules(stats, job_segments):
    """The modules in the order they run.

    stats is built as the jobs are read, so its own order is first encounter
    across the batch: a job that fails early, or one whose task list differs
    from the first job's, inserts its modules out of sequence.

    Within a job the segments are in order, so every consecutive pair is a
    precedence fact, and a topological sort over those facts orders modules
    that no single job ever ran together. A median start offset alone is not
    enough: a job that analyzed one task reaches postqc hours earlier than a
    job that analyzed four, so postqc's median lands among the tasks. The
    offset is used only to break ties among the modules that are ready at the
    same point, and to break a cycle where jobs disagree on the order.
    """
    successors = {name: set() for name in stats}
    indegree = {name: 0 for name in stats}
    edges = set()
    for _, segments in job_segments:
        names = [segment['module'] for segment in segments]
        for first, second in zip(names, names[1:]):
            if first == second or (first, second) in edges:
                continue
            if first not in successors or second not in indegree:
                continue
            edges.add((first, second))
            successors[first].add(second)
            indegree[second] += 1

    def chronological(name):
        entry = stats[name]
        offset = (statistics.median(entry['offsets']) if entry['offsets']
                  else float('inf'))
        return (offset, statistics.median(entry['positions']), name)

    ordered = collections.OrderedDict()
    remaining = set(stats)
    while remaining:
        ready = [name for name in remaining if indegree[name] == 0]
        if not ready:
            # the jobs disagree on the order, so drop what holds back the
            #   earliest module and carry on rather than stopping
            ready = [min(remaining, key=chronological)]
        pick = min(ready, key=chronological)
        ordered[pick] = stats[pick]
        remaining.discard(pick)
        for name in successors[pick]:
            if name in remaining:
                indegree[name] = max(0, indegree[name] - 1)
    return ordered


def format_module_times(stats):
    """The per-module run times in hours, one row per module.

    A module that runs more than once in a job, as the analysis step does per
    task, contributes one sample per run. # runs counts every time the module
    was reached and # timed only the runs that stamped both a startdate and a
    procdate, so a module skipped because its output was already there is
    reported with no times of its own.
    """
    if not stats:
        return ''
    width = max(len(name) for name in stats)
    width = max(width, len('module'))
    lines = ['', 'module run-time (hours), from the status messages:',
             '  %-*s %7s %8s %8s %8s %8s %8s %7s'
             % (width, 'module', '# runs', '# timed', 'min', 'max', 'median',
                'mean', '# err')]
    for name, entry in stats.items():
        hours = [t / 3600.0 for t in entry['times']]
        summary = summarize(hours)
        if summary is None:
            lines.append('  %-*s %7d %8d %8s %8s %8s %8s %7d'
                         % (width, name, entry['runs'], 0, '-', '-', '-', '-',
                            entry['errors']))
            continue
        lines.append('  %-*s %7d %8d %8.4f %8.4f %8.4f %8.4f %7d'
                     % (width, name, entry['runs'], summary['n'],
                        summary['min'], summary['max'], summary['median'],
                        summary['mean'], entry['errors']))
    return '\n'.join(lines)


def account_jobs(batchdir):
    """One accounting record per job in the scriptlist.

    Each record carries the job's name, its resolved .out and .err paths, its
    state and whether its .err indicates failure. The .out text is read and
    discarded inside the loop rather than kept, since a single .out runs to
    thousands of lines.
    """
    pbsout_dir = os.path.join(batchdir, 'pbsout')

    # a batch that was built and never submitted has no pbsout/ at all, which
    # is a state to report rather than a reason to quit
    if os.path.isdir(pbsout_dir):
        entries = os.listdir(pbsout_dir)
    else:
        entries = []
        print(f'WARNING: directory not found: {pbsout_dir}', file=sys.stderr)

    has_markers = markers_present(pbsout_dir)
    records = []
    for job_name in read_scriptlist(batchdir):
        out_path, err_path = resolve_job_files(pbsout_dir, job_name, entries)
        err_text = read_file(err_path) if err_path else ''
        state = job_state(pbsout_dir, job_name, out_path, err_path,
                          has_markers)
        start, finish, basis = job_times(pbsout_dir, job_name, state,
                                         out_path, err_path)
        header = out_header(out_path)
        records.append({
            'name': job_name,
            'out': out_path,
            'err': err_path,
            'state': state,
            'errored': job_errored(err_text),
            'start': start,
            'finish': finish,
            'basis': basis,
            'host': header.get('HOSTNAME', ''),
            'jobid': header.get('JOB_ID', ''),
            'cluster': clu.cluster_from_host(header.get('HOSTNAME')) or '',
            'usage': None,
        })
    return records


def add_usage(records):
    """Each record's accounting, one query per cluster; a cluster that
    does not answer leaves its records without usage."""
    by_cluster = collections.OrderedDict()
    for record in records:
        if record['cluster'] and record['jobid']:
            by_cluster.setdefault(record['cluster'], []).append(record)
    for cluster, rows in by_cluster.items():
        scheduler = clu.CLUSTER_SCHEDULERS.get(cluster)
        if scheduler is None:
            continue
        try:
            usage = clu.job_usage(cluster, scheduler,
                                  [r['jobid'] for r in rows])
        except ValueError as e:
            print(f'WARNING: {cluster}: {e}', file=sys.stderr)
            continue
        if usage is None:
            print(f'WARNING: no accounting from {cluster}', file=sys.stderr)
            continue
        for record in rows:
            record['usage'] = usage.get(record['jobid'])


def quartiles(values):
    """25th, 50th and 75th percentiles, interpolated between ranks."""
    xs = sorted(values)

    def pct(frac):
        pos = frac * (len(xs) - 1)
        low = int(pos)
        high = min(low + 1, len(xs) - 1)
        return xs[low] + (xs[high] - xs[low]) * (pos - low)
    return pct(0.25), pct(0.5), pct(0.75)


def usage_rows(records):
    """{cluster: [usage]} for the finished jobs with accounting."""
    rows = collections.OrderedDict()
    for record in sorted(records, key=lambda r: r['cluster']):
        if record['state'] == STATE_FINISHED and record['usage']:
            rows.setdefault(record['cluster'], []).append(record['usage'])
    return rows


def format_usage(records):
    """Memory, cores, CPU time and disk I/O per cluster, median (IQR),
    for the finished jobs with accounting."""
    rows = usage_rows(records)
    if not rows:
        return ''
    lines = ['', 'memory, CPU and disk I/O per finished job, from the '
             'scheduler accounting, median (IQR):']
    for cluster, usage in rows.items():
        scheduler = clu.CLUSTER_SCHEDULERS.get(cluster, '')
        lines.append(f'  {cluster} ({scheduler}): {len(usage)} jobs')
        items = [('elapsed hours', 'elapsed_h'), ('cores', 'cores'),
                 ('CPU hours', 'cpu_h')]
        for mem_type in sorted(set(u['mem_type'] for u in usage)):
            items.append(('memory GB, %s' % MEM_LABELS.get(mem_type, mem_type),
                          ('mem_gb', mem_type)))
        items += [('read GB', 'read_gb'), ('write GB', 'write_gb')]
        for label, key in items:
            if isinstance(key, tuple):
                values = [u[key[0]] for u in usage
                          if u['mem_type'] == key[1] and u[key[0]] is not None]
            else:
                values = [u[key] for u in usage if u[key] is not None]
            # read_gb is None on SGE
            if not values:
                continue
            q1, median, q3 = quartiles(values)
            lines.append('    %-32s %9.2f (%.2f-%.2f)' % (label, median, q1, q3))
        writes = [u['write_gb'] for u in usage if u['write_gb'] is not None]
        reads = [u['read_gb'] for u in usage if u['read_gb'] is not None]
        if writes:
            lines.append('    %-32s %9.1f' % ('total write GB', sum(writes)))
        if reads:
            lines.append('    %-32s %9.1f' % ('total read GB', sum(reads)))
    n_missing = sum(1 for r in records
                    if r['state'] == STATE_FINISHED and not r['usage'])
    if n_missing:
        lines.append('  # finished jobs without accounting: %d' % n_missing)
    return '\n'.join(lines)


def format_accounting(batchdir, records):
    """The job counts, in mmil_check_batchdir's order.

    A job is counted once by state; the error count crosses all three, since
    the wrappers finish a job whatever its body did.
    """
    counts = {state: 0 for state in STATE_LABELS}
    for record in records:
        counts[record['state']] += 1
    lines = [f'summary of jobs for {os.path.basename(batchdir.rstrip("/"))}:',
             f'  # of jobs: {len(records)}']
    for state in (STATE_FINISHED, STATE_STARTED, STATE_NEVER_STARTED):
        lines.append(f'  # {STATE_LABELS[state]}: {counts[state]}')
    lines.append('  # with errors: %d'
                 % sum(1 for r in records if r['errored']))
    return '\n'.join(lines)


STALE_ROOTS = ('raw_rootdir', 'proc_rootdir', 'fsurf_rootdir',
               'analysis_rootdir', 'postqc_rootdir')


def batch_sessions(records):
    """(participant_id, session_id) of the batch's jobs, from their names
    (job_NNN_sub-X_ses-Y), in job order."""
    sessions = []
    for record in records:
        m = re.search(r'(sub-[^_]+)_(ses-[^_]+)', record['name'])
        if m and m.groups() not in sessions:
            sessions.append(m.groups())
    return sessions


def format_stale(records, dirs):
    """Warnings for the batch's sessions left stale: the full-pipeline check
    of cmig_mkjobs_bids_pipeline (a step whose input ran after it, whose
    events copies or scans.tsv selection changed, or downstream of such a
    step), as a job that runs only some steps leaves the steps after them."""
    import cmig_mkjobs_bids_pipeline as mkjobs
    sessions = batch_sessions(records)
    lines = []
    for pid, sid in sessions:
        status, why = mkjobs.check_session('cmig_pipeline_bids', pid, sid,
                                           dirs, {})
        if status == 'stale':
            lines.append(f'  WARNING: {pid} {sid}: {why or "stale"}')
    return '\n'.join([f'sessions left stale: {len(lines)} of {len(sessions)}']
                     + lines)


def run(args):
    batchdir = resolve_batchdir(args.batchdir, args.root_batchdir)
    if not os.path.isdir(batchdir):
        sys.exit(f'Error: batchdir not found: {args.batchdir}')

    records = account_jobs(batchdir)
    if args.usageflag:
        add_usage(records)
    summary = [format_accounting(batchdir, records)]
    runtimes = format_runtimes(records)
    if runtimes:
        summary.append(runtimes)
    usage = format_usage(records)
    if usage:
        summary.append(usage)
    dirs = dict((k, getattr(args, k)) for k in STALE_ROOTS if getattr(args, k))
    if dirs.get('raw_rootdir') and len(dirs) > 1:
        summary.append(format_stale(records, dirs))

    jobs_by_category = {cat: [] for cat in CATEGORY_ORDER}
    module_times = collections.OrderedDict()
    job_segments = []
    module_aliases = {}

    for record in records:
        err_text = read_file(record['err']) if record['err'] else ''
        out_text = read_file(record['out']) if record['out'] else ''

        # a job that is not a cmig_bids module writes no status blocks, and
        # then the overall timing and the .err are all there is
        blocks = parse_status_blocks(out_text)
        collect_module_aliases(blocks, module_aliases)
        segments = job_module_segments(blocks, parse_banners(out_text),
                                       record['start'], err_text)
        add_module_times(module_times, segments)
        if segments:
            job_segments.append((record['name'], segments))

        cat = categorize(err_text, record['state'], blocks)
        record['category'] = cat
        record['errmsg'] = tsv_errmsg(err_text)
        excerpt = err_excerpt(err_text, cat)
        ctx = extract_context(out_text, cat)
        jobs_by_category[cat].append((record['name'], excerpt, ctx))

    # the table and every module plot read this dict in order, so ordering it
    #   once here puts all of them in running order
    module_times = resolve_module_aliases(module_times, job_segments,
                                          module_aliases)
    module_times = order_modules(module_times, job_segments)
    modules = format_module_times(module_times)
    if modules:
        summary.append(modules)
    summary_text = '\n'.join(summary)
    print(summary_text)

    batchdir_name = os.path.basename(batchdir.rstrip('/'))
    if not args.outdir:
        return
    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.join(args.outdir, batchdir_name)
    write_summary(stem + '_summary.txt', summary_text)
    write_jobs_tsv(stem + '_jobs.tsv', records)
    if args.plotflag:
        make_plots(records, args.outdir, batchdir_name)
        make_module_plots(job_segments, module_times, records, args.outdir,
                          batchdir_name)
        make_usage_plots(records, args.outdir, batchdir_name)
    outpath = os.path.join(args.outdir, f'{batchdir_name}_error_summary.txt')

    total = len(records)

    # signature sub-rows share the count column, so they set the width too
    widths = [len(CATEGORY_LABELS[c]) for c in CATEGORY_ORDER]
    sig_groups = {}
    for cat in CATEGORY_ORDER:
        if cat == 'success':
            continue
        groups = group_by_signature(jobs_by_category[cat])
        sig_groups[cat] = groups
        if len(groups) > 1:
            widths.extend(len(sig) + 2 for sig, _ in groups)
    label_width = max(widths)

    with open(outpath, 'w') as f:
        f.write(f'Error Summary: {batchdir_name}\n')
        f.write(f'{total} jobs total\n\n')

        # Count table
        f.write(f'{"Category":<{label_width}}  Count\n')
        f.write('-' * (label_width + 8) + '\n')
        for cat in CATEGORY_ORDER:
            count = len(jobs_by_category[cat])
            if count == 0:
                continue
            f.write(f'{CATEGORY_LABELS[cat]:<{label_width}}  {count:>5}\n')
            groups = sig_groups.get(cat, [])
            if len(groups) > 1:
                for sig, sig_jobs in groups:
                    f.write(f'{"  " + sig:<{label_width}}  {len(sig_jobs):>5}\n')
        f.write('\n')

        # Per-category detail sections (skip success)
        for cat in CATEGORY_ORDER:
            if cat == 'success':
                continue
            jobs = jobs_by_category[cat]
            if not jobs:
                continue

            sep = '=' * 55
            f.write(f'{sep}\n')
            f.write(f'{CATEGORY_LABELS[cat]} ({len(jobs)} jobs)\n')
            f.write(f'{sep}\n')

            groups = sig_groups.get(cat, [])
            for sig, sig_jobs in groups:
                # one block per root cause, so a category holding several
                # unrelated failures does not read as a single problem
                if len(groups) > 1:
                    f.write(f'-- {sig} ({len(sig_jobs)} jobs)\n')

                f.write('Jobs:\n')
                for job_name, _, _ in sig_jobs:
                    f.write(f'  {job_name}\n')
                f.write('\n')

                job_name, excerpt, ctx = sig_jobs[0]
                ctx_str = format_context(ctx)
                if excerpt.strip() or ctx_str:
                    f.write(f'Example: {job_name}\n')
                    if excerpt.strip():
                        f.write('  .err:\n')
                        for line in excerpt.splitlines():
                            f.write(f'    {line}\n')
                    if ctx_str:
                        f.write('  .out:\n')
                        f.write(ctx_str + '\n')
                    f.write('\n')

    print(f'wrote {outpath}')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Report on a batch job directory.'
    )
    parser.add_argument('--batchdir', required=True,
                        help='batch job directory holding pbsout/, either a '
                             'full path or a name under --root_batchdir')
    parser.add_argument('--root_batchdir',
                        default=os.path.join('~', 'batchdirs'),
                        help='root to resolve a bare --batchdir name against '
                             '(default ~/batchdirs)')
    parser.add_argument('--outdir',
                        help='output directory; if omitted, results are '
                             'printed and no files are written')
    parser.add_argument('--plotflag', type=int, default=1, choices=(0, 1),
                        help='whether to make plots (default 1)')
    parser.add_argument('--usageflag', type=int, default=1, choices=(0, 1),
                        help='whether to query the scheduler accounting '
                             '(ssh to the head nodes) for memory, CPU and '
                             'disk I/O per job (default 1)')
    # with raw_rootdir and the output roots, the batch's sessions are checked
    #   for steps left stale (e.g. by a job that ran only some steps)
    for key in STALE_ROOTS:
        parser.add_argument('--' + key,
                            help=key + ' of the batch, for the check of the '
                                 'steps left stale (with --raw_rootdir and '
                                 'at least one other root)')
    return parser.parse_args()


def main():
    run(parse_args())


if __name__ == '__main__':
    main()
