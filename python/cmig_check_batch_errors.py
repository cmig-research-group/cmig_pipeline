#!/usr/bin/env python3
"""
Analyze .err/.out files from a CMIG batch job directory and write an error summary.

Created:  05/18/2026 by Don Hagler
Prev mod: 05/31/2026 by Don Hagler
Prev mod: 06/02/2026 by Don Hagler
Prev mod: 08/19/2026 by Don Hagler
Last mod: 10/19/2026 by Don Hagler
"""

import argparse
import os
import re
import sys


N_TAIL = 10

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
    'uncounted_error',
    'running',
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
    'uncounted_error': 'Pipeline exit 0 with a module failure it does not count (errcode=1)',
    'running':        'Still running or incomplete (no errcode in .out)',
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


def pipeline_result(out_text):
    """(finished, errcode) from a pipeline job's closing lines,
    '{pipeline}: errcode = N' and '{pipeline}: finished'; (False, None) for a
    job of a single module."""
    errcodes = re.findall(r'^cmig_pipeline_bids\w*: errcode = (\S+)\s*$',
                          out_text, re.MULTILINE)
    finished = bool(re.search(r'^cmig_pipeline_bids\w*: finished\s*$',
                              out_text, re.MULTILINE))
    return finished, (errcodes[-1] if errcodes else None)


def categorize(err_text, out_text=''):
    # a pipeline that finished with errcode 0 succeeded, but can carry a
    #   module failure it does not count (an analyze func task while another
    #   task succeeded), seen in the .err as that module's errcode=1
    finished, errcode = pipeline_result(out_text)
    if finished and errcode == '0':
        if 'errcode=1' in err_text:
            return 'uncounted_error'
        if err_text.strip():
            return 'success_benign'
        return 'success'
    # a .err holding nothing but benign output is not evidence of failure, so
    # the job is judged from the .out; kept in its own category rather than
    # merged into success so the noise stays visible in the summary
    if not strip_benign(err_text).strip():
        finished = bool(re.search(r'^finished\s*$', out_text, re.MULTILINE))
        errcodes = re.findall(r'errcode = (\S+)', out_text)
        if finished and errcodes and errcodes[-1] == '0':
            if err_text.strip():
                return 'success_benign'
            return 'success'
        return 'running'
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


def run(args):
    pbsout_dir = os.path.join(args.batchdir, 'pbsout')
    if not os.path.isdir(pbsout_dir):
        sys.exit(f'Error: directory not found: {pbsout_dir}')

    err_files = sorted(f for f in os.listdir(pbsout_dir) if f.endswith('.err'))
    if not err_files:
        sys.exit(f'Error: no .err files found in {pbsout_dir}')

    jobs_by_category = {cat: [] for cat in CATEGORY_ORDER}

    for err_fname in err_files:
        job_name = err_fname[:-4]
        err_text = read_file(os.path.join(pbsout_dir, err_fname))
        out_text = read_file(os.path.join(pbsout_dir, job_name + '.out'))

        cat = categorize(err_text, out_text)
        excerpt = err_excerpt(err_text, cat)
        ctx = extract_context(out_text, cat)
        jobs_by_category[cat].append((job_name, excerpt, ctx))

    os.makedirs(args.outdir, exist_ok=True)
    batchdir_name = os.path.basename(args.batchdir.rstrip('/'))
    outpath = os.path.join(args.outdir, f'{batchdir_name}_error_summary.txt')

    total = sum(len(v) for v in jobs_by_category.values())

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
                f.write(f'Example: {job_name}\n')
                if excerpt.strip():
                    f.write('  .err:\n')
                    for line in excerpt.splitlines():
                        f.write(f'    {line}\n')
                ctx_str = format_context(ctx)
                if ctx_str:
                    f.write('  .out:\n')
                    f.write(ctx_str + '\n')
                f.write('\n')

    print(f'wrote {outpath}')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Analyze batch job .err/.out files and write an error summary.'
    )
    parser.add_argument('--batchdir', required=True,
                        help='Full path to the batch job directory (contains pbsout/)')
    parser.add_argument('--outdir', required=True,
                        help='Output directory for the summary file')
    return parser.parse_args()


def main():
    run(parse_args())


if __name__ == '__main__':
    main()
