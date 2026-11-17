#!/usr/bin/env python3
"""
BIDS app front end of the cmig bids container (run.py is a symbolic link to
this script): runs the pipelines, or one module or pipeline by name, for the
sessions of a BIDS dataset, and the group-level steps across them. The
MATLAB modules run through cmig_run_bids_module.py, and so through the
single compiled app cmig_run_bids_module (compiled_flag is always 1).

Usage:
  cmig_run_bids_app.py bids_dir output_dir {participant,group} [options]
  cmig_run_bids_app.py --module NAME -h     the options of NAME
  cmig_run_bids_app.py -h | -v

participant level: for each session, the rawqc pipeline
  (cmig_pipeline_bids_rawqc), then the full pipeline (cmig_pipeline_bids),
  which runs even if rawqc failed; the session fails if either fails. With
  --module NAME, that module or pipeline alone.
group level: cmig_compile_bids_tabdat (the tabulated data of the rawqc,
  analysis and postqc roots, into --tabulated_dir), then
  cmig_concat_bids_sessions for each datatype (the analysis root, into
  --concat_dir); a datatype with no data is skipped with a warning unless
  named with --datatype. Both run with --forceflag 1 unless --forceflag 0 is
  given. With --module cmig_concat_bids_sessions or cmig_compile_bids_tabdat,
  that alone. --module bids_validator (only on request) runs bids-validator
  on bids_dir and on each per-session output root given by its option (not
  tabulated or concatenated), its JSON report
  into output_dir/logs/; --skip_bids_validator is accepted and ignored.

Sessions: --participant_label and --session_label take one or more labels,
with or without the sub-/ses- prefix (repeats dropped); without them, every
sub-* in bids_dir and every ses-* of each. bids_dir must have a session
level (sub-X/ses-Y/). A step that does not read rawdata (fsurf, analyze,
postqc) finds its sessions in the proc root instead. Each session runs on
its own, a failure not stopping the rest; the exit status is 1 if any
session (or group step) failed. At the group level the labels select the
sessions combined; without them, every session found in the output roots.

Roots: each output root defaults to a subdirectory of output_dir, and is
passed to the steps that use it under the step's own option name (the root
a step writes as its --output_rootdir); bids_dir is every step's
--raw_rootdir. A root a step only reads is passed only where it exists,
unless the step requires it; a root given that is only read must exist.
Protocol compliance (cmig_pc_bids_abcd) writes into --pc_rootdir, default
the rawqc root. dataset_description.json is written into each per-session
output root the run writes (where absent; not into the group roots
tabulated and concatenated), and the ABCD T1 atlas into the proc root.

Step options: every option run.py does not take itself is passed to each
step to be run whose parser accepts it (run.py --module NAME -h lists
them); an option no such step accepts is an error, before anything runs.
--participant_id, --session_id, --raw_rootdir, --output_rootdir,
--compiled_flag and the group steps' --outdir, --sublist and --seslist are
set by run.py and cannot be given. A value that begins with - is given as
--option=value (e.g. --fsurf_options=-all).

--fs_license names the FreeSurfer license file (license.txt), a path the
container can see (bind it in); else FS_LICENSE. The participant level
stops before anything runs if neither names a readable file. --work_dir
sets SCRATCHROOT and MCR_CACHE_ROOT (work_dir/mcr_cache) for every step.

Versions: every status file records the image tag (CMIG_IMAGE_TAG) as
app_version, CMIG_IMAGE_SHA256 (the .sif's sha256, which mkjobs and the
rounds set in the job; n/a where not set) as app_sha256, and the
package versions of CMIG_VERSIONS_FILE (versions.json); the steps inherit
them from run.py's environment. dataset_description.json GeneratedBy lists
the same packages.

Created:  10/06/2026 by Don Hagler
Prev mod: 10/06/2026 by Don Hagler
Prev mod: 10/07/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/13/2026 by Don Hagler
"""

import argparse
import csv
import importlib
import glob
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback

import cmig_run_bids_module

LEVELS = ['participant', 'group']
PIPELINES = ['cmig_pipeline_bids_rawqc', 'cmig_pipeline_bids']
# the python modules of a session (the MATLAB ones are those of
#   cmig_run_bids_module)
PY_MODULES = ['cmig_proc_bids_dwi', 'cmig_fsurf_bids', 'cmig_analyze_bids_behav',
              'cmig_pc_bids_abcd', 'cmig_check_event_pcqc']
CONCAT = 'cmig_concat_bids_sessions'
TABDAT = 'cmig_compile_bids_tabdat'
VALIDATOR = 'bids_validator'
VALIDATOR_EXE = 'bids-validator'
# the modules run across sessions, at the group level
GROUP_MODULES = [CONCAT, TABDAT, VALIDATOR]

# root option, and its default subdirectory of output_dir
ROOTS = [('rawqc_rootdir', 'cmig_rawqc'), ('proc_rootdir', 'cmig_proc'),
         ('fsurf_rootdir', 'freesurfer'), ('analysis_rootdir', 'cmig_analysis'),
         ('postqc_rootdir', 'cmig_postqc'), ('tabulated_dir', 'tabulated'),
         ('concat_dir', 'concatenated')]
SESSION_ROOTS = ['rawqc_rootdir', 'proc_rootdir', 'fsurf_rootdir',
                 'analysis_rootdir', 'postqc_rootdir']
# the roots of cmig_write_bids_dataset_description's table, and those the
#   validator checks (the group roots are not laid out as BIDS datasets)
DD_ROOTS = SESSION_ROOTS + ['pc_rootdir']
# root labels, for the validator reports and the Name fallback
ROOT_LABELS = dict(ROOTS, pc_rootdir='cmig_pc')
# the roots of the per-session tabulated data the group level compiles
TABDAT_ROOTS = ['rawqc_rootdir', 'analysis_rootdir', 'postqc_rootdir']

# the root each step writes, by step (a rawqc alias resolves to
#   cmig_rawqc_bids); a pipeline missing here writes several, by name
OUTPUT_ROOT = {
    'cmig_pipeline_bids_rawqc': 'rawqc_rootdir',
    'cmig_rawqc_bids': 'rawqc_rootdir',
    'cmig_pc_bids_abcd': 'pc_rootdir',
    'cmig_check_event_pcqc': 'rawqc_rootdir',
    'cmig_proc_bids_anat': 'proc_rootdir',
    'cmig_proc_bids_func': 'proc_rootdir',
    'cmig_proc_bids_dwi': 'proc_rootdir',
    'cmig_fsurf_bids': 'fsurf_rootdir',
    'cmig_analyze_bids_anat': 'analysis_rootdir',
    'cmig_analyze_bids_func': 'analysis_rootdir',
    'cmig_analyze_bids_dwi': 'analysis_rootdir',
    'cmig_analyze_bids_behav': 'analysis_rootdir',
    'cmig_postqc_bids_anat': 'postqc_rootdir',
    'cmig_postqc_bids_func': 'postqc_rootdir',
    'cmig_postqc_bids_dwi': 'postqc_rootdir',
    'cmig_postqc_bids_fsurf': 'postqc_rootdir',
    'cmig_postqc_bids_fmap': 'postqc_rootdir',
    CONCAT: 'concat_dir',
    TABDAT: 'tabulated_dir',
}
WRITTEN_ROOTS = {'cmig_pipeline_bids': ['proc_rootdir', 'fsurf_rootdir',
                                        'analysis_rootdir', 'postqc_rootdir'],
                 'cmig_pipeline_bids_rawqc': ['rawqc_rootdir', 'pc_rootdir'],
                 VALIDATOR: []}
# the root a step that does not read rawdata finds its sessions in (as
#   mkjobs scans); the others scan bids_dir
SCAN_ROOT = dict((n, 'proc_rootdir') for n in [
    'cmig_fsurf_bids', 'cmig_analyze_bids_anat', 'cmig_analyze_bids_func',
    'cmig_analyze_bids_dwi', 'cmig_postqc_bids_anat', 'cmig_postqc_bids_func',
    'cmig_postqc_bids_dwi', 'cmig_postqc_bids_fsurf', 'cmig_postqc_bids_fmap'])

# options set by run.py, never given
RESERVED = ['participant_id', 'session_id', 'raw_rootdir', 'output_rootdir',
            'compiled_flag', 'outdir', 'sublist', 'seslist']

# per concat datatype: the analysis subdirectory and the per-session files
#   it reads (any one is data); one table, shared with the module's check
CONCAT_DATA = cmig_run_bids_module.CONCAT_DATA
# the visits csv's place in a step's arguments until it is written
VISITS_CSV = '<visits csv>'

class Terminated(BaseException):
    """SIGTERM (a scheduler's kill), raised past run_step so the finally
    blocks run."""


def on_sigterm(signum, frame):
    raise Terminated()


def list_names():
    """Return the modules and pipelines run.py runs by name."""
    return sorted(set(cmig_run_bids_module.list_names()) | set(PY_MODULES)
                  | set(PIPELINES) | set([TABDAT, VALIDATOR]))


def is_matlab(name):
    return name in cmig_run_bids_module.list_names()


def base_name(name):
    """The module name with a rawqc alias resolved."""
    return cmig_run_bids_module.resolve(name)[0] if is_matlab(name) else name


def written_roots(name):
    """The root options a step writes."""
    if name in WRITTEN_ROOTS:
        return WRITTEN_ROOTS[name]
    return [OUTPUT_ROOT[base_name(name)]]


def step_parser(name, argv=()):
    """The parser of a step's arguments, named as run.py shows it."""
    if is_matlab(name):
        parser = cmig_run_bids_module.build_parser(name, argv)
    else:
        parser = importlib.import_module(name).build_parser()
    parser.prog = 'run.py --module ' + name
    return parser


def parser_options(parser):
    """The option names (without --) a parser accepts."""
    return set(s[2:] for s in parser._option_string_actions if s.startswith('--'))


def required_options(parser):
    return set(a.dest for a in parser._actions if a.required and a.option_strings)


def build_parser():
    prog = os.path.basename(sys.argv[0])
    parser = argparse.ArgumentParser(
        prog=prog, add_help=False, allow_abbrev=False,
        usage='%(prog)s bids_dir output_dir {participant,group} [options]')
    parser.add_argument('bids_dir', help='BIDS raw data directory (read only)')
    parser.add_argument('output_dir',
                        help='output directory, holding the output roots '
                             'unless they are given')
    parser.add_argument('analysis_level', choices=LEVELS,
                        help='participant: per session; group: across sessions')
    sel = parser.add_argument_group('selection')
    sel.add_argument('--participant_label', nargs='+', default=None,
                     help='participant labels, with or without sub- '
                          '(default: every sub-* in bids_dir)')
    sel.add_argument('--session_label', nargs='+', default=None,
                     help='session labels, with or without ses- (default: '
                          'every ses-* of each participant)')
    sel.add_argument('--module', default=None,
                     help='run this module or pipeline alone (see names '
                          'below); at the group level, only '
                          + ', '.join(GROUP_MODULES))
    sel.add_argument('--datatype', nargs='+', default=None,
                     help='group level: the cmig_concat_bids_sessions '
                          'datatypes (default: all, those with no data '
                          'skipped)')
    roots = parser.add_argument_group('output roots')
    for key, sub in ROOTS:
        roots.add_argument('--' + key, default=None,
                           help='(default: output_dir/{0})'.format(sub))
    roots.add_argument('--pc_rootdir', default=None,
                       help='protocol compliance (default: the rawqc root)')
    other = parser.add_argument_group('other')
    other.add_argument('--fs_license', default=None,
                       help='FreeSurfer license file, set as FS_LICENSE '
                            '(default: $FS_LICENSE); the participant level '
                            'requires one of them')
    other.add_argument('--work_dir', default=None,
                       help='temporary files: sets SCRATCHROOT=DIR and '
                            'MCR_CACHE_ROOT=DIR/mcr_cache for every step, '
                            'created if missing (default: the image\'s, /tmp)')
    other.add_argument('--skip_bids_validator', action='store_true',
                       help='accepted and ignored: the validator runs only '
                            'as the group-level --module ' + VALIDATOR)
    other.add_argument('--dataset_name', default=None,
                       help='Name in the dataset_description.json written into '
                            'each per-session output root, suffixed with the root\'s '
                            'label (default: the Name of bids_dir\'s '
                            'dataset_description.json, else its basename)')
    other.add_argument('--dry_flag', type=int, default=0, choices=[0, 1],
                       help='1: print each step\'s arguments per session, '
                            'run nothing and write nothing (default: 0)')
    other.add_argument('-v', '--version', action='store_true',
                       help='print the cmig_pipeline version and exit')
    other.add_argument('-h', '--help', action='store_true',
                       help='this help; with --module NAME, the options of NAME')
    return parser


def version():
    return 'cmig_pipeline {0}'.format(
        os.environ.get('CMIG_PIPELINE_VERSION') or 'unknown (CMIG_PIPELINE_VERSION not set)')


def usage(parser):
    prog = os.path.basename(sys.argv[0])
    text = __doc__.split('Created:')[0].rstrip()
    print(text.replace('cmig_run_bids_app.py', prog) + '\n')
    print(parser.format_help().rstrip())
    print('\nnames (--module):\n  ' + '\n  '.join(list_names()))


def check_module(name, level):
    """The --module name, checked against the level; exits on error."""
    if name.endswith('.py'):
        name = name[:-3]
    if name not in list_names():
        sys.exit('Error: unknown module: {0}; run.py -h lists them'.format(name))
    if level == 'group' and name not in GROUP_MODULES:
        sys.exit('Error: {0} is a participant-level module; the group level '
                 'runs only {1}'.format(name, ', '.join(GROUP_MODULES)))
    if level == 'participant' and name in GROUP_MODULES:
        sys.exit('Error: {0} runs across sessions, at the group level'.format(name))
    return name


def strip_labels(labels, prefix):
    """The labels without the prefix, repeats dropped, in the order given."""
    if labels is None:
        return None
    out = []
    for lab in labels:
        lab = lab[len(prefix):] if lab.startswith(prefix) else lab
        if lab not in out:
            out.append(lab)
    return out


def subdirs(dirname, prefix):
    try:
        return sorted(d for d in os.listdir(dirname) if d.startswith(prefix)
                      and os.path.isdir(os.path.join(dirname, d)))
    except OSError:
        return []


def find_sessions(rootdir, participant_labels, session_labels, strict):
    """[(pid, sid)] of the sub-*/ses-* directories of rootdir that the labels
    select. strict: a participant label without a directory is an error."""
    if participant_labels is None:
        pids = subdirs(rootdir, 'sub-')
    else:
        pids = []
        for lab in participant_labels:
            pid = 'sub-' + lab
            if os.path.isdir(os.path.join(rootdir, pid)):
                pids.append(pid)
            elif strict:
                sys.exit('Error: participant not found: {0}'.format(
                    os.path.join(rootdir, pid)))
    sessions = []
    nmissing = 0
    what = 'selected ses-*' if session_labels else 'ses-*'
    for pid in pids:
        sids = subdirs(os.path.join(rootdir, pid), 'ses-')
        if session_labels is not None:
            sids = [s for s in sids if s[4:] in session_labels]
        if not sids and strict:
            if participant_labels is None:
                nmissing += 1
            else:
                print('WARNING: {0}: no {1} directory in {2}'.format(
                    pid, what, rootdir))
        sessions += [(pid, sid) for sid in sids]
    if nmissing:
        # one line for a scanned root, not one per participant
        sel = ' ({0})'.format(' '.join('ses-' + s for s in session_labels)) \
            if session_labels else ''
        print('WARNING: {0} of {1} participants have no {2} directory{3} in '
              '{4}'.format(nmissing, len(pids), what, sel, rootdir))
    return sessions


def check_session_level(bids_dir):
    """Exit if bids_dir has sub-* directories but no sub-*/ses-*."""
    pids = subdirs(bids_dir, 'sub-')
    if pids and not any(subdirs(os.path.join(bids_dir, p), 'ses-') for p in pids):
        sys.exit('Error: {0} has no ses-* directories; cmig_pipeline requires '
                 'a session level (sub-X/ses-Y/)'.format(bids_dir))


def option_groups(argv):
    """[(name, tokens)] of the options given as --name [values] or
    --name=value, in order; exits on a value with no option."""
    groups = []
    for tok in argv:
        if tok.startswith('--') and len(tok) > 2:
            groups.append((tok[2:].split('=', 1)[0], [tok]))
        elif groups:
            groups[-1][1].append(tok)
        else:
            sys.exit('Error: unexpected argument: {0}'.format(tok))
    return groups


def route_options(groups, parsers):
    """Check that each option given is accepted by a step to be run, and is
    not one run.py sets; exits on error."""
    errors = []
    for name, _ in groups:
        if name in RESERVED:
            errors.append('--{0} is set by run.py'.format(name))
        elif not any(name in opts for opts in parsers.values()):
            errors.append('--{0} is accepted by none of the steps to be run '
                          '({1})'.format(name, ', '.join(parsers)))
    if errors:
        sys.exit('Error: ' + '\nError: '.join(errors))


def as_tokens(val):
    if isinstance(val, (list, tuple)):
        return [str(v) for v in val]
    return [str(val)]


def step_argv(opts, required, fixed, defaults, groups, positional=()):
    """A step's argument list: run.py's values for the options it takes,
    the defaults of the required ones not given, then the options given that
    it accepts."""
    argv = list(positional)
    given = set(name for name, _ in groups)
    for key, val in fixed.items():
        if key in opts and val is not None:
            argv += ['--' + key] + as_tokens(val)
    for key, val in defaults.items():
        if key in opts and key in required and key not in given and key not in fixed:
            argv += ['--' + key] + as_tokens(val)
    for name, toks in groups:
        if name in opts and name not in fixed:
            argv += toks
    return argv


def session_fixed(name, roots, bids_dir, pid=None, sid=None, absent=()):
    """The values run.py sets for a step, by option name; the roots in absent
    are left off."""
    vals = {'raw_rootdir': bids_dir, 'compiled_flag': 1}
    for key in SESSION_ROOTS + ['pc_rootdir']:
        if key not in absent:
            vals[key] = roots[key]
    if pid:
        vals['participant_id'] = pid
        vals['session_id'] = sid
    out = OUTPUT_ROOT.get(base_name(name))
    if out:
        vals['output_rootdir'] = roots[out]
    return vals


def parse_step(name, parser, argv):
    """The parsed arguments of a step: run_module kwargs for a MATLAB module,
    else the namespace its run() takes; argparse exits on an error."""
    args = parser.parse_args(argv)
    if is_matlab(name):
        return cmig_run_bids_module.parsed_kwargs(name, args)
    return args


def run_step(name, parsed):
    """Run one step; return (errcode, errmsg). A step that raises or exits
    is a failure, not the end of the run."""
    print('\n' + '#' * 60)
    print('run.py: {0}'.format(name))
    print('#' * 60, flush=True)
    # run_module reports a MATLAB module's errcode and errmsg itself (on
    #   stderr when it failed); not when it raised or exited
    reported = is_matlab(name)
    try:
        if is_matlab(name):
            result = cmig_run_bids_module.run_module(name, **parsed)
        else:
            result = importlib.import_module(name).run(parsed)
        if isinstance(result, tuple):
            errcode, errmsg = result[0], result[1]
        else:
            errcode, errmsg = result, None
    except SystemExit as e:
        reported = False
        errcode = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
        errmsg = None if isinstance(e.code, int) or e.code is None else str(e.code)
    except Exception as e:
        reported = False
        traceback.print_exc()
        errcode, errmsg = 1, '{0}: {1}'.format(type(e).__name__, e)
    errcode = 1 if errcode else 0
    sys.stdout.flush()
    if not reported:
        # the closing lines cmig_check_batch_errors.py judges a pipeline job
        #   by; a failure also on stderr, where it looks for one
        print('{0}: errcode = {1}'.format(name, errcode))
        if errcode:
            print('{0}: ERROR: {1}'.format(name, errmsg or 'errcode 1'),
                  file=sys.stderr, flush=True)
        elif errmsg:
            print('{0}: warning: {1}'.format(name, errmsg))
    print('{0}: finished'.format(name), flush=True)
    return errcode, errmsg


def dataset_name(bids_dir):
    """bids_dir's Name without a trailing "Raw data" or root label (bids_dir
    is the input root in an image job of a step that does not read rawdata),
    else bids_dir's basename."""
    try:
        with open(os.path.join(bids_dir, 'dataset_description.json')) as f:
            name = json.load(f).get('Name') or ''
        labels = ['Raw data'] + sorted(ROOT_LABELS.values())
        name = re.sub(r'(^|\s+-\s+)({0})$'.format('|'.join(map(re.escape, labels))),
                      '', name.strip(), flags=re.IGNORECASE).strip()
        if name:
            return name
    except (OSError, ValueError, AttributeError):
        pass
    return os.path.basename(os.path.normpath(bids_dir))


def written_keys(names, roots):
    """The root options the steps write; the pc root is the rawqc root
    unless given apart from it."""
    keys = []
    for n in names:
        for k in written_roots(n):
            if k == 'pc_rootdir' and roots[k] == roots['rawqc_rootdir']:
                k = 'rawqc_rootdir'
            if k not in keys:
                keys.append(k)
    return keys


def write_root_files(keys, roots, bids_dir, name, FS_version=None):
    """dataset_description.json into each per-session root written, and the
    atlas into the proc root; a failure is a warning."""
    import cmig_write_bids_dataset_description as dd
    # none into the group roots (tabulated, concatenated)
    keys = [k for k in keys if k in DD_ROOTS]
    for key in keys:
        msg = dd.check_root_name(roots[key], 'fsurf' if key == 'fsurf_rootdir'
                                 else 'cmig')
        if msg:
            print('WARNING: {0}'.format(msg))
    try:
        if keys:
            kwargs = {}
            if FS_version is not None:
                kwargs['FS_version'] = FS_version
            dd.write_dataset_descriptions(
                name, write_roots=keys, raw_rootdir=bids_dir,
                **dict([(k, roots[k]) for k in DD_ROOTS], **kwargs))
    except Exception as e:
        print('WARNING: dataset_description.json not written: {0}'.format(e))
    if 'proc_rootdir' in keys:
        try:
            import cmig_write_bids_atlas
            cmig_write_bids_atlas.write_bids_atlas(roots['proc_rootdir'])
        except Exception as e:
            print('WARNING: atlas not written into {0}: {1}'.format(
                roots['proc_rootdir'], e))


def read_visits(fname):
    """The (participant_id, session_id) pairs of an id_events csv."""
    with open(fname) as f:
        return set((r['participant_id'].strip(), r['session_id'].strip())
                   for r in csv.DictReader(f))


def given_values(groups):
    """{name: [values]} of the step options given."""
    return dict((n, ([t[0].split('=', 1)[1]] if '=' in t[0] else []) + t[1:])
                for n, t in groups)


def input_roots(names, parsers, roots, given_roots):
    """{step: the roots to leave off}: a root a step only reads is passed
    where it exists; exits if a given one, or one the step requires, does
    not exist."""
    errors = []
    absent = {}
    for n in names:
        opts = parser_options(parsers[n])
        req = required_options(parsers[n])
        written = set(roots[k] for k in written_roots(n))
        absent[n] = set()
        for key in SESSION_ROOTS + ['pc_rootdir']:
            if key not in opts or roots[key] in written or os.path.isdir(roots[key]):
                continue
            if key in given_roots:
                msg = '--{0} not found: {1}'.format(key, roots[key])
            elif key in req:
                msg = '--{0} not found: {1} ({2} reads it)'.format(key, roots[key], n)
            else:
                absent[n].add(key)
                continue
            if msg not in errors:
                errors.append(msg)
    if errors:
        sys.exit('Error: ' + '\nError: '.join(errors))
    return absent


def participant_jobs(args, roots, groups, names, given_roots):
    """[(label, [(step, argv, parsed)])] per session, every step parsed."""
    pids = strip_labels(args.participant_label, 'sub-')
    sids = strip_labels(args.session_label, 'ses-')
    parsers = dict((n, step_parser(n, sum([t for _, t in groups], []))) for n in names)
    route_options(groups, dict((n, parser_options(p)) for n, p in parsers.items()))
    absent = input_roots(names, parsers, roots, given_roots)
    # a step that does not read rawdata finds its sessions in its input root
    scan_key = SCAN_ROOT.get(base_name(names[0])) if len(names) == 1 else None
    scan_dir = roots[scan_key] if scan_key else args.bids_dir
    sessions = find_sessions(scan_dir, pids, sids, strict=True)
    if not sessions:
        sys.exit('Error: no sessions selected in {0}'.format(scan_dir))
    defaults = {'manualqc_rootdir': roots['rawqc_rootdir']}
    jobs = []
    for pid, sid in sessions:
        steps = []
        for n in names:
            p = parsers[n]
            argv = step_argv(parser_options(p), required_options(p),
                             session_fixed(n, roots, args.bids_dir, pid, sid,
                                           absent[n]),
                             defaults, groups)
            steps.append((n, argv, parse_step(n, p, argv)))
        jobs.append(('{0} {1}'.format(pid, sid), steps))
    return jobs


def has_concat_data(rootdir, sessions, datatypes):
    """The datatypes with a per-session file of theirs in one of the sessions
    (stops at the first session that has each)."""
    left = list(datatypes)
    found = []
    for pid, sid in sessions:
        for d in list(left):
            subdir, patterns = CONCAT_DATA[d]
            indir = os.path.join(rootdir, pid, sid, subdir)
            if any(glob.glob(os.path.join(indir, '{0}_{1}{2}'.format(pid, sid, pat)))
                   for pat in patterns):
                found.append(d)
                left.remove(d)
        if not left:
            break
    return found


def selected_sessions(rootdir, pids, sids):
    if not os.path.isdir(rootdir):
        return []
    return find_sessions(rootdir, pids, sids, strict=False)


def validator_jobs(args, roots, groups, given_roots):
    """[('group', [('_validator', label, dir, report)])]: bids_dir and each
    per-session output root given."""
    if groups:
        sys.exit('Error: --module {0} takes no step options: {1}'.format(
            VALIDATOR, ' '.join('--' + n for n, _ in groups)))
    stamp = time.strftime('%y%m%d%H%M')
    logdir = os.path.join(args.output_dir, 'logs')
    targets = [('rawdata', args.bids_dir)]
    for key in DD_ROOTS:
        if key in given_roots and roots[key] not in [t[1] for t in targets]:
            if not os.path.isdir(roots[key]):
                sys.exit('Error: --{0} not found: {1}'.format(key, roots[key]))
            targets.append((ROOT_LABELS[key], roots[key]))
    steps = [('_validator', label, d, os.path.join(
        logdir, 'bids_validator_{0}_{1}.json'.format(label, stamp)))
        for label, d in targets]
    return [('group', steps)]


def group_jobs(args, roots, groups, names, given_roots):
    """[(label, [(step, argv, parsed)])] of the group steps, every step
    parsed."""
    if names == [VALIDATOR]:
        return validator_jobs(args, roots, groups, given_roots)
    import cmig_mkjobs_bids_concat
    datatypes = args.datatype or sorted(cmig_mkjobs_bids_concat.DATATYPE_CONFIG)
    unknown = [d for d in datatypes if d not in cmig_mkjobs_bids_concat.DATATYPE_CONFIG]
    if unknown:
        sys.exit('Error: unknown datatype: {0}'.format(' '.join(unknown)))
    pids = strip_labels(args.participant_label, 'sub-')
    sids = strip_labels(args.session_label, 'ses-')
    labels = bool(pids or sids)
    parsers = dict((n, step_parser(n, sum([t for _, t in groups], []))) for n in names)
    route_options(groups, dict((n, parser_options(p)) for n, p in parsers.items()))
    given = given_values(groups)
    # a root given that is only read must exist
    read_keys = (TABDAT_ROOTS if TABDAT in names else []) + \
        (['analysis_rootdir'] if CONCAT in names else [])
    for key in read_keys:
        if key in given_roots and not os.path.isdir(roots[key]):
            sys.exit('Error: --{0} not found: {1}'.format(key, roots[key]))
    # forceflag 1 unless given: a rerun rewrites the tables and concat files
    force = {} if 'forceflag' in given else {'forceflag': 1}
    if TABDAT in names and labels and not any(
            selected_sessions(roots[k], pids, sids) for k in TABDAT_ROOTS):
        sys.exit('Error: no selected sessions in {0}, {1} or {2}'.format(
            *[roots[k] for k in TABDAT_ROOTS]))
    steps = []
    if TABDAT in names:
        p = parsers[TABDAT]
        fixed = dict(force, outdir=roots['tabulated_dir'],
                     sublist=['sub-' + x for x in pids] if pids else None,
                     seslist=['ses-' + x for x in sids] if sids else None)
        defaults = {'phenotype_dir': os.path.join(args.bids_dir, 'phenotype')}
        if 'phenotype_dir' not in given and not os.path.isdir(defaults['phenotype_dir']):
            sys.exit('Error: no phenotype directory in bids_dir ({0}); give '
                     '--phenotype_dir, the release\'s per-table json '
                     'directory'.format(defaults['phenotype_dir']))
        rootdirs = [roots[k] for k in TABDAT_ROOTS if os.path.isdir(roots[k])]
        if not rootdirs:
            sys.exit('Error: none of the roots of tabulated data exists: {0}'.format(
                ' '.join(roots[k] for k in TABDAT_ROOTS)))
        argv = step_argv(parser_options(p), required_options(p), fixed, defaults,
                         groups, positional=rootdirs)
        steps.append((TABDAT, argv, parse_step(TABDAT, p, argv)))
    if CONCAT in names:
        steps += concat_steps(args, roots, groups, parsers[CONCAT], given, force,
                              datatypes, pids, sids, names == [CONCAT])
    return [('group', steps)]


def concat_steps(args, roots, groups, p, given, force, datatypes, pids, sids,
                 alone):
    """The concat steps, one per datatype with data; alone: --module concat,
    where no selected session is an error."""
    adir = roots['analysis_rootdir']
    fixed = dict(session_fixed(CONCAT, roots, args.bids_dir), **force)
    fixed['seslist'] = ['ses-' + x for x in sids] if sids else None
    visits = set(selected_sessions(adir, pids, sids))
    if pids and 'id_events_file' in given:
        visits &= read_visits(given['id_events_file'][0])
    if (pids or sids) and not visits:
        msg = 'no selected sessions in {0}'.format(adir)
        if alone:
            sys.exit('Error: ' + msg)
        print('WARNING: {0}; concat skipped'.format(msg))
        return []
    steps = []
    if pids:
        # the module selects participants only through an id_events csv: one
        #   listing the selected sessions of the analysis root (those of a
        #   given --id_events_file too) is written for the run
        fixed['id_events_file'] = VISITS_CSV
        steps.append(('_visits', None, sorted(visits)))
    if not args.datatype:
        found = has_concat_data(adir, sorted(visits), datatypes)
        for d in datatypes:
            if d not in found:
                print('WARNING: no {0} data in {1}; concat {0} skipped'.format(d, adir))
        datatypes = [d for d in datatypes if d in found]
    if not datatypes:
        return []
    for datatype in datatypes:
        fixed['datatype'] = datatype
        argv = step_argv(parser_options(p), required_options(p), fixed, {}, groups)
        steps.append((CONCAT, argv, parse_step(CONCAT, p, argv)))
    return steps


def scratch_dir():
    """SCRATCHROOT where it exists, else python's ($TMPDIR, /tmp)."""
    d = os.environ.get('SCRATCHROOT')
    return d if d and os.path.isdir(d) else tempfile.gettempdir()


def write_visits(visits):
    """Write the visits csv to a new temporary file; return its name."""
    fd, fname = tempfile.mkstemp(prefix='run_py_visits_', suffix='.csv',
                                 dir=scratch_dir())
    with os.fdopen(fd, 'w') as f:
        f.write('participant_id,session_id\n')
        for pid, sid in visits:
            f.write('{0},{1}\n'.format(pid, sid))
    return fname


def run_validator(label, dirname, fname):
    """bids-validator --json on one dataset, stdout into fname; return its
    exit code."""
    print('\n' + '#' * 60)
    print('run.py: {0} {1}'.format(VALIDATOR, dirname))
    print('#' * 60, flush=True)
    os.makedirs(os.path.dirname(fname), exist_ok=True)
    env = dict(os.environ)
    # deno's cache, where a container's HOME may not be writable
    env.setdefault('DENO_DIR', os.path.join(scratch_dir(), 'deno_dir'))
    cmd = [shutil.which(VALIDATOR_EXE), '--json', dirname]
    print(' '.join(cmd), flush=True)
    with open(fname, 'w') as f:
        errcode = subprocess.call(cmd, stdout=f, env=env)
    name = '{0} {1}'.format(VALIDATOR, label)
    print('{0}: report {1}'.format(name, fname))
    print('{0}: errcode = {1}'.format(name, errcode))
    if errcode:
        print('{0}: ERROR: {1} exited {2}; see {3}'.format(
            name, VALIDATOR_EXE, errcode, fname), file=sys.stderr, flush=True)
    print('{0}: finished'.format(name), flush=True)
    return errcode


def print_dry(jobs, level, roots):
    for label, steps in jobs:
        print('{0}:'.format(label))
        for step in steps:
            if step[0] == '_visits':
                print('  (writes a temporary visits csv into {0}: {1} '
                      'sessions)'.format(scratch_dir(), len(step[2])))
            elif step[0] == '_validator':
                print('  {0} --json {1} > {2}'.format(VALIDATOR_EXE, step[2], step[3]))
            else:
                print('  {0} {1}'.format(step[0], ' '.join(step[1])))
    print('{0} {1}, roots:'.format(len(jobs), 'sessions' if level == 'participant'
                                    else 'group run'))
    for key in [k for k, _ in ROOTS] + ['pc_rootdir']:
        print('  --{0} {1}'.format(key, roots[key]))


def check_fs_license(fs_license, level):
    """Set FS_LICENSE from --fs_license, else keep it; the participant level
    exits without a readable license file."""
    lic = fs_license or os.environ.get('FS_LICENSE')
    ok = bool(lic) and os.path.isfile(lic) and os.access(lic, os.R_OK)
    if fs_license and not ok:
        sys.exit('Error: FreeSurfer license file not found or not readable: '
                 '{0}'.format(fs_license))
    if level == 'participant' and not ok:
        sys.exit('Error: a FreeSurfer license is required: --fs_license FILE '
                 'or FS_LICENSE')
    if ok:
        os.environ['FS_LICENSE'] = os.path.abspath(lic)


def set_work_dir(work_dir, dry_flag):
    """SCRATCHROOT and MCR_CACHE_ROOT for every step, created (not with
    dry_flag)."""
    work_dir = os.path.abspath(work_dir)
    os.environ['SCRATCHROOT'] = work_dir
    os.environ['MCR_CACHE_ROOT'] = os.path.join(work_dir, 'mcr_cache')
    if not dry_flag:
        os.makedirs(os.environ['MCR_CACHE_ROOT'], exist_ok=True)


def run_jobs(jobs):
    """Run the steps; return [(label, [failed steps])] and the validator's
    first nonzero exit code (0 if none)."""
    failed = []
    vcode = 0
    visits_csv = None
    try:
        for label, steps in jobs:
            errs = []
            for step in steps:
                if step[0] == '_visits':
                    visits_csv = write_visits(step[2])
                    continue
                if step[0] == '_validator':
                    errcode = run_validator(*step[1:])
                    if errcode:
                        vcode = vcode or errcode
                        errs.append('{0} ({1})'.format(VALIDATOR, step[1]))
                    continue
                name, _, parsed = step
                if is_matlab(name) and parsed.get('id_events_file') == VISITS_CSV:
                    parsed['id_events_file'] = visits_csv
                errcode, errmsg = run_step(name, parsed)
                if errcode:
                    errs.append(name if name != CONCAT
                                else '{0} ({1})'.format(name, parsed['datatype']))
            if errs:
                failed.append((label, errs))
    finally:
        if visits_csv and os.path.isfile(visits_csv):
            os.remove(visits_csv)
    return failed, vcode


def main():
    parser = build_parser()
    argv = sys.argv[1:]
    # the options of a step, or run.py's help, before the positionals are
    #   required
    if not argv or '-h' in argv or '--help' in argv:
        pre = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
        pre.add_argument('--module', default=None)
        known, rest = pre.parse_known_args(argv)
        if known.module:
            name = known.module[:-3] if known.module.endswith('.py') else known.module
            if name not in list_names():
                sys.exit('Error: unknown module: {0}'.format(name))
            if name == VALIDATOR:
                print('run.py bids_dir output_dir group --module {0} '
                      '[--<root> DIR ...]\n\nruns {1} --json on bids_dir and on '
                      'each output root given by its option; each report goes '
                      'to output_dir/logs/{0}_<root label>_<YYMMDDHHMM>.json; '
                      'the exit status is the validator\'s (the first nonzero '
                      'one). No step options.'.format(VALIDATOR, VALIDATOR_EXE))
                sys.exit(0)
            step_parser(name, [a for a in rest if a not in ('-h', '--help')]).print_help()
            sys.exit(0)
        usage(parser)
        sys.exit(0 if argv else 1)
    if '-v' in argv or '--version' in argv:
        print(version())
        sys.exit(0)
    if any(a.startswith('-') for a in argv[:3]):
        sys.exit('Error: bids_dir, output_dir and analysis_level come first '
                 '(run.py -h)')
    args, rest = parser.parse_known_args(argv)

    args.bids_dir = os.path.abspath(args.bids_dir)
    args.output_dir = os.path.abspath(args.output_dir)
    if not os.path.isdir(args.bids_dir):
        sys.exit('Error: bids_dir not found: {0}'.format(args.bids_dir))
    roots = {}
    given_roots = set()
    for key, sub in ROOTS:
        val = getattr(args, key)
        if val:
            given_roots.add(key)
        roots[key] = os.path.abspath(val) if val else os.path.join(args.output_dir, sub)
    if args.pc_rootdir:
        given_roots.add('pc_rootdir')
    roots['pc_rootdir'] = (os.path.abspath(args.pc_rootdir) if args.pc_rootdir
                           else roots['rawqc_rootdir'])

    level = args.analysis_level
    if args.datatype and level != 'group':
        sys.exit('Error: --datatype is a group-level option')
    if args.module:
        names = [check_module(args.module, level)]
    elif level == 'participant':
        names = list(PIPELINES)
    else:
        names = [TABDAT, CONCAT]
    if names != [VALIDATOR]:
        check_session_level(args.bids_dir)
    # fail now, not hours later at the first FreeSurfer tool
    check_fs_license(args.fs_license, level)
    if args.work_dir:
        set_work_dir(args.work_dir, args.dry_flag)

    groups = option_groups(rest)
    if level == 'participant':
        jobs = participant_jobs(args, roots, groups, names, given_roots)
    else:
        jobs = group_jobs(args, roots, groups, names, given_roots)

    if args.dry_flag:
        print_dry(jobs, level, roots)
        sys.exit(0)
    if names == [VALIDATOR] and not shutil.which(VALIDATOR_EXE):
        sys.exit('Error: {0} not found on PATH (the cmig image has it)'.format(
            VALIDATOR_EXE))

    FS_version = given_values(groups).get('FS_version')
    write_root_files(written_keys(names, roots), roots, args.bids_dir,
                     args.dataset_name or dataset_name(args.bids_dir),
                     FS_version[0] if FS_version else None)

    signal.signal(signal.SIGTERM, on_sigterm)
    try:
        failed, vcode = run_jobs(jobs)
    except Terminated:
        sys.exit('Error: run.py terminated (SIGTERM)')

    print('\nrun.py: {0} of {1} {2} failed'.format(
        len(failed), len(jobs), 'sessions' if level == 'participant' else 'group runs'))
    for label, errs in failed:
        print('  {0}: {1}'.format(label, ', '.join(errs)))
    if vcode:
        sys.exit(vcode)
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
