#!/usr/bin/env python3
"""
Generic python wrapper for the cmig bids MATLAB modules: runs any of them,
by name, with --compiled_flag 1 (the default) through the single compiled app
cmig_run_bids_module, or with 0 as MATLAB code (a script run with
run_matlab_script, cmig_pipeline on the MATLAB path; $CMIG_MATLAB_ADDPATH, dirs
separated by ':', adds those dirs first, e.g. working copies). With --reproducible_flag
1 (the default) the module runs with OpenBLAS and MKL settings and a fixed
MATLAB thread count that give the same results on every CPU type and host.

Usage:
  cmig_run_bids_module.py <module> [module arguments] [runtime arguments]
  cmig_run_bids_module.py <module> --help

The module arguments are the module's positional inputs, given as named
options (e.g. --participant_id), and the options documented in the module's
MATLAB help. An option not given is not passed, so the module's own default
applies. List values are given space-separated and passed to MATLAB as a cell
array or a numeric vector.

From python, run_module(module, **kwargs) does the same and returns
(errcode, errmsg).

Created:  10/13/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/14/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 10/23/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/03/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Prev mod: 11/13/2026 by Don Hagler
Prev mod: 11/14/2026 by Don Hagler
Last mod: 11/15/2026 by Don Hagler
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

import cmig_sing_utils
import cmig_write_scripts

DISPATCH_APPNAME = 'cmig_run_bids_module'
DEFAULT_APPVER = 'mmps_dev'
DEFAULT_APPDIR = '/usr/pubsw/packages/MMPS/apps'
DEFAULT_PUBSDIR = '/usr/pubsw'
DEFAULT_MRTDIR = '/usr/pubsw/packages/matlab/R2023a'
DEFAULT_MRTDIR_TARGET = '/opt/matlabruntime/R2023a'
RUNTIME_DEFAULTS = {
    'compiled_flag': 1,
    'reproducible_flag': 1,
    'appver': DEFAULT_APPVER,
    'appdir': DEFAULT_APPDIR,
    'pubsdir': DEFAULT_PUBSDIR,
    'mrtdir': DEFAULT_MRTDIR,
    'mrtdir_target': DEFAULT_MRTDIR_TARGET,
}

REPRODUCIBLE_HELP = ('1 = the same results on every CPU type and host: sets '
                     'OPENBLAS_CORETYPE=HASWELL, OPENBLAS_NUM_THREADS=1, '
                     'MKL_CBWR=AVX2 and CMIG_MATLAB_NTHREADS=1 (MATLAB '
                     'threads) where not already set (needs a CPU with '
                     'AVX2); 0 = the inherited environment, for a CPU where 1 '
                     'fails (default: %(default)s)')

# how an option's value is passed to MATLAB: as given, as a numeric vector
#   [1 2], as a cell array of char {'a' 'b'}, or as given with none meaning []
NUM = 'num'
CELL = 'cell'
NONE_EMPTY = 'none_empty'

STATUS = '{output_rootdir}/{participant_id}/{session_id}/{modality}/{participant_id}_{session_id}_{entities}desc-{type}_status.tsv'
POSTQC_ROOTS = ['participant_id', 'session_id', 'proc_rootdir', 'fsurf_rootdir',
                'output_rootdir']
POSTQC_OPTIONS = {'cleanup_flag': '', 'verbose': '', 'forceflag': ''}

# per module: the compiled MATLAB function (app), its positional inputs in
#   order, those that may be left empty (optional), the options documented in
#   its help with how their values are passed, the options a call must give
#   (required), the root dirs it only reads (indirs), the positional inputs
#   that take several values, passed as a cell array (multi), and its status
#   file (modality and type, or a function of the arguments)
MODULES = {
    'cmig_proc_bids_anat': {
        'positional': ['participant_id', 'session_id', 'raw_rootdir', 'output_rootdir'],
        'options': {'select_flag': '', 'atlas_flag': '', 'native_flag': '',
                    'gradunwarp_flag': '', 'cleanup_flag': '', 'forceflag': ''},
        'indirs': ['raw_rootdir'],
        'status': ('anat', 'proc'),
    },
    'cmig_proc_bids_func': {
        'positional': ['participant_id', 'session_id', 'raw_rootdir', 'output_rootdir'],
        'options': {'B0unwarp_flag': '', 'gradunwarp_flag': '', 'native_flag': '',
                    'nvoxels': NUM, 'resolution': NUM, 'std_orient': NONE_EMPTY,
                    'regT1flag': '', 'cleanup_flag': '', 'forceflag': ''},
        'indirs': ['raw_rootdir'],
        'status': ('func', 'proc'),
    },
    'cmig_analyze_bids_anat': {
        'positional': POSTQC_ROOTS,
        'options': {'cleanup_flag': '', 'forceflag': ''},
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': ('anat', 'analyze'),
    },
    'cmig_analyze_bids_func': {
        'positional': POSTQC_ROOTS,
        'options': {
            'task': '', 'deconv_flag': '', 'stim_times_flag': '', 'block_gap': '',
            'minfrac': '', 'max_dur_std': '', 'stim_times_model': '',
            'stim_times_model_mid': '', 'stim_times_model_sst': '',
            'stim_times_model_nback': '', 'censor_thresh_deconv': '',
            'censor_thresh_corr': '', 'contig_thresh_corr': '', 'contig_thresh': '',
            'min_numTRs': '', 'min_ntpoints': '', 'prefilt_flag': '',
            'mc_resp_flag': '', 'resp_low': '', 'resp_high': '', 'skipTRs': '',
            **{f'ntpoints_{t}': '' for t in ['rest', 'mid', 'nback', 'sst']},
            'TR': '', 'pydir': '', 'ignore_delay_flag': '',
            'exclude_noncompliant_flag': '', 'overshoot_thresh': '',
            'max_delay': '', 'short_flag': '', 'norm_flag': '', 'norm_thresh': '',
            'despike_flag': '', 'baseline_flag': '', 'njobs': '', 'outliers_flag': '',
            'outliers_thresh': '', 'prep_tseries_flag': '', 'extract_ico_order': NUM,
            'corr_ico_order': NUM, 'trim_minutes': NUM, 'rest_trim_minutes': NUM,
            'trim_rand_flag': '', 'cleanup_flag': '', 'forceflag': ''},
        # the status file is named by task
        'required': ['task'],
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': lambda a: ('func', 'analyze', f'task-{a["task"]}_'),
    },
    'cmig_analyze_bids_dwi': {
        'positional': POSTQC_ROOTS,
        'options': {
            'dti_flag': '', 'rsi_flag': '', 'fiber_flag': '', 'aseg_flag': '',
            'aseg_disp_flag': '', 'aseg_dispfact': '', 'aseg_disp_suffix': '',
            'aseg_disp_roicodes': NUM, 'aseg_dispvec': NUM, 'xcg_flag': '',
            'thresh_prob': '', 'inner_dti_flag': '', 'inner_bval_max': '',
            'dt_outfix_inner': '', 'gwcsurf_outfix_inner': '', 'dt_outfix_full': '',
            'rsi_outfix': '', 'dt_measlist_full': CELL, 'dt_scalefacts_full': NUM,
            'reregister_flag': '', 'wmbc_flag': '', 'wmbc_vol_target': '',
            'wmbc_image_range': NUM, 'wmbc_wm_codes': NUM, 'wmbc_gm_codes': NUM,
            'wmbc_excl_codes': NUM, 'rsi_measlist': CELL, 'motion_radius': '',
            'fibers': NUM, 'trans_fibers': NUM, 'fiber_atlasdir': '', 'atlasdir': '',
            'tabdir': '', 'verbose': '', 'cleanup_flag': '', 'forceflag': ''},
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': ('dwi', 'analyze'),
    },
    'cmig_rawqc_bids': {
        'positional': ['participant_id', 'session_id', 'raw_rootdir', 'output_rootdir'],
        'options': {'modality': '', 'cleanup_flag': '', 'export_flag': '', 'verbose': '',
                    'forceflag': '', 'motion_radius': '',
                    'method': '', 'fMRI_minTRs': '', 'afni_flag': '',
                    'censor_reg_flag': ''},
        # the status file is in the modality subdirectory
        'required': ['modality'],
        'indirs': ['raw_rootdir'],
        'status': lambda a: (a['modality'], 'rawqc', ''),
    },
    'cmig_postqc_bids_anat': {
        'positional': POSTQC_ROOTS,
        'optional': ['fsurf_rootdir'],
        'options': POSTQC_OPTIONS,
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': ('anat', 'postqc'),
    },
    'cmig_postqc_bids_func': {
        'positional': POSTQC_ROOTS,
        'optional': ['fsurf_rootdir'],
        'options': dict(POSTQC_OPTIONS, skipTRs='', skipTRs_GE_DV25='', skipTRs_GE='',
                        skipTRs_Philips='', skipTRs_Siemens='', fMRI_topval_prctile=''),
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': ('func', 'postqc'),
    },
    'cmig_postqc_bids_dwi': {
        'positional': POSTQC_ROOTS,
        'optional': ['fsurf_rootdir'],
        'options': dict(POSTQC_OPTIONS, min_ndirs=''),
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': ('dwi', 'postqc'),
    },
    'cmig_postqc_bids_fsurf': {
        'positional': POSTQC_ROOTS,
        'options': POSTQC_OPTIONS,
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        'status': ('fsurf', 'postqc'),
    },
    'cmig_postqc_bids_fmap': {
        'positional': POSTQC_ROOTS,
        'optional': ['fsurf_rootdir'],
        'options': dict(POSTQC_OPTIONS, acq_types=CELL),
        'indirs': ['proc_rootdir', 'fsurf_rootdir'],
        # one acquisition type is named in the status file
        'status': lambda a: ('fmap', 'postqc',
                             f'acq-{_as_list(a["acq_types"])[0]}_'
                             if a.get('acq_types') is not None
                             and len(_as_list(a['acq_types'])) == 1 else ''),
    },
    'cmig_concat_bids_sessions': {
        'positional': ['datatype', 'analysis_rootdir', 'output_rootdir'],
        'options': dict(
            {'seslist': CELL, 'id_events_file': '',
             'measlist': CELL, 'layerlist': CELL, 'surflist': CELL,
             'smoothlist': NUM, 'hemilist': CELL, 'tasklist': CELL, 'atlaslist': CELL,
             'ico_order': '', 'vol_info_flag': '', 'forceflag': ''},
            **{f'{n}_{t}': '' for n in ['nruns', 'ntpoints']
               for t in ['rest', 'mid', 'nback', 'sst']}),
        # sessions, a struct array, cannot be given on the command line
        'excluded': ['sessions'],
        'indirs': ['analysis_rootdir'],
        # several analysis roots, their sessions combined
        'multi': ['analysis_rootdir'],
        # a cohort-level function: no per-session status file; its invocation
        #   status is read by check_concat_output
        'status': None,
    },
}

# the names of the earlier per-modality rawqc wrappers, each one module with
#   its modality set
ALIASES = {f'cmig_rawqc_bids_{m}': ('cmig_rawqc_bids', {'modality': m})
           for m in ['anat', 'func', 'dwi', 'fmap']}


def resolve(name):
    """Return the module name and the options an alias fixes."""
    if name.endswith('.py'):
        name = name[:-3]
    if name in ALIASES:
        return ALIASES[name]
    if name not in MODULES:
        raise ValueError(f'unknown module: {name}; expected one of: '
                         + ', '.join(list_names()))
    return name, {}


# option descriptions and defaults from the MATLAB help, written beside the
#   compiled app by cmig_build_run_bids_module.py, since the MATLAB code is not
#   installed with it
HELP_NAME = 'cmig_run_bids_module_help.json'


def load_help(appdir, appver):
    """Return the option help of the compiled app, or {} if it has none."""
    fname = os.path.join(
        cmig_sing_utils._compiled_subdir(appdir, DISPATCH_APPNAME, appver), HELP_NAME)
    try:
        with open(fname) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def list_names():
    return sorted(list(MODULES) + list(ALIASES))


def _as_list(val):
    if isinstance(val, (list, tuple)):
        return list(val)
    return str(val).split()


def encode(val, kind):
    """Return the string MATLAB's cmig_parseval turns into the intended value."""
    if kind == NUM:
        return '[' + ' '.join(str(v) for v in _as_list(val)) + ']'
    if kind == CELL:
        # each element quoted: {ses-00A ses-02A} is not a cell of char
        return '{' + ' '.join(f"'{v}'" for v in _as_list(val)) + '}'
    if kind == NONE_EMPTY and str(val).lower() == 'none':
        return '[]'
    if isinstance(val, bool):
        return str(int(val))
    return str(val)


def typed_value(val, kind):
    """Return val as the python value cmig_write_scripts writes as a MATLAB literal."""
    if kind == NUM:
        return [_number(v) for v in _as_list(val)]
    if kind == CELL:
        return [str(v) for v in _as_list(val)]
    if kind == NONE_EMPTY and str(val).lower() == 'none':
        return None
    return _number(val)


def _number(val):
    """Return val as an int or float if it reads as one, else as given."""
    if isinstance(val, bool):
        return int(val)
    if isinstance(val, (int, float)):
        return val
    for conv in (int, float):
        try:
            return conv(str(val))
        except ValueError:
            pass
    return val


def status_file(module, kwargs):
    """Return the module's status file path, or None if it writes none."""
    spec = MODULES[module]['status']
    if spec is None:
        return None
    if callable(spec):
        modality, stype, entities = spec(kwargs)
    else:
        (modality, stype), entities = spec, ''
    return STATUS.format(output_rootdir=kwargs['output_rootdir'],
                         participant_id=kwargs['participant_id'],
                         session_id=kwargs['session_id'],
                         modality=modality, type=stype, entities=entities)


# per concat datatype: the analysis subdirectory and the per-session files
#   it reads, any one of which is data (as cmig_concat_bids_sessions
#   check_datatype); cmig_run_bids_app.py has_concat_data uses it
CONCAT_DATA = {
    'smri_measures': ('anat', ['*_hemi-*_atlas-ico*_desc-fwhm*.mat']),
    'smri_surfaces': ('anat', ['*_hemi-*_atlas-ico*_white.mat',
                               '*_hemi-*_atlas-ico*_pial.mat']),
    'dmri_dti': ('dwi', ['*_desc-fwhm*_dti*.mat']),
    'dmri_isdti': ('dwi', ['*_desc-fwhm*_isdti*.mat']),
    'dmri_rsi': ('dwi', ['*_desc-fwhm*_rsi*.mat']),
    'tfmri_stats': ('func', ['*_hemi-*_atlas-ico*_desc-*_beta.mat']),
    'fmri_corr': ('func', ['*_task-*_atlas-*_corrmat.mat']),
    'fmri_tseries': ('func', ['*_task-*_atlas-*_tseries.mat']),
}


# the concat datatypes cmig_mkjobs_bids_concat.py expands by measure and
#   smoothing, whose invocation status is then named by both
CONCAT_EXPANDABLE = ['smri_measures', 'dmri_dti', 'dmri_isdti', 'dmri_rsi']


def concat_status_name(kwargs):
    """The invocation status file name of cmig_concat_bids_sessions, as its
    status_stem names it: concat_{datatype}_status.tsv, or
    concat_{datatype}_{meas}_sm{smooth}_status.tsv for an expandable datatype
    given one measure and one smoothing step."""
    datatype = kwargs['datatype']
    measlist = kwargs.get('measlist')
    smoothlist = kwargs.get('smoothlist')
    if isinstance(measlist, str):
        measlist = [measlist]
    if smoothlist is not None and not isinstance(smoothlist, (list, tuple)):
        smoothlist = [smoothlist]
    if (datatype in CONCAT_EXPANDABLE and measlist and smoothlist
            and len(measlist) == 1 and len(smoothlist) == 1):
        # %g, as the module formats it
        stem = 'concat_{0}_{1}_sm{2:g}'.format(datatype, measlist[0],
                                               float(smoothlist[0]))
    else:
        stem = 'concat_{0}'.format(datatype)
    return stem + '_status.tsv'


def concat_status_file(kwargs):
    """The invocation status file, in output_rootdir (the module writes it
    there, since its outputs go to the release 7.0 layout folders, two of
    them for fmri_corr and fmri_tseries). None if it is not there."""
    fname = os.path.join(kwargs['output_rootdir'], concat_status_name(kwargs))
    if not os.path.isfile(fname):
        return None
    return fname


def check_concat_output(kwargs, app_name=None, app_version=None):
    """Report the result of cmig_concat_bids_sessions and return (errcode, errmsg).

    Read from its invocation status file (concat_status_file), as check_status
    reads a session's, with the version columns appended. A no-data skip has
    errcode 0 and its errmsg, "no {datatype} data in ...".
    """
    module = 'cmig_concat_bids_sessions'
    fname = concat_status_file(kwargs)
    if fname is None:
        errmsg = 'missing status file ' + os.path.join(
            kwargs['output_rootdir'], concat_status_name(kwargs))
        print(module + ': errcode=1', file=sys.stderr)
        print(module + ': errmsg=' + errmsg, file=sys.stderr)
        return 1, errmsg
    return cmig_sing_utils.check_status(fname, module, app_name=app_name,
                                        app_version=app_version)


def run_module(name, **kwargs):
    """Run one module; return (errcode, errmsg).

    kwargs are the module's positional inputs and options, by name, and the
    runtime settings (compiled_flag, reproducible_flag, appver, appdir,
    pubsdir, mrtdir, mrtdir_target). Options set to None are not passed.
    compiled_flag: 0 = MATLAB code, 1 = compiled app (default)
    reproducible_flag: 1 = set OPENBLAS_CORETYPE, OPENBLAS_NUM_THREADS,
      MKL_CBWR and CMIG_MATLAB_NTHREADS where not already set (default),
      0 = the inherited environment
    """
    module, fixed = resolve(name)
    spec = MODULES[module]
    runtime = dict(RUNTIME_DEFAULTS)
    for key in list(kwargs):
        if key in RUNTIME_DEFAULTS:
            val = kwargs.pop(key)
            if val is not None:
                runtime[key] = val
    kwargs.update(fixed)
    unknown = [k for k in kwargs
               if k not in spec['positional'] and k not in spec['options']]
    if unknown:
        raise ValueError(f'{module}: unknown option(s): {", ".join(unknown)}')
    optional = spec.get('optional', [])
    for key in spec['positional']:
        if kwargs.get(key) is None:
            if key not in optional:
                raise ValueError(f'{module}: missing input: {key}')
            kwargs[key] = ''
    for key in spec.get('required', []):
        if kwargs.get(key) is None:
            raise ValueError(f'{module}: missing option: {key}')

    multi = spec.get('multi', [])
    for key in spec['indirs']:
        dirs = kwargs[key] if isinstance(kwargs[key], (list, tuple)) else [kwargs[key]]
        for d in dirs:
            if d and not os.path.isdir(d):
                sys.exit(f'Error: directory not found: {d}')
    os.makedirs(kwargs['output_rootdir'], exist_ok=True)

    # a list given for a multi input goes as a cell array, a single value as is
    def is_list(key):
        return key in multi and isinstance(kwargs[key], (list, tuple))
    container_args = [encode(kwargs[key], CELL) if is_list(key) else str(kwargs[key])
                      for key in spec['positional']]
    for key, kind in spec['options'].items():
        val = kwargs.get(key)
        if val is None:
            continue
        container_args += [key, encode(val, kind)]

    appver = runtime['appver']
    compiled_flag = int(runtime['compiled_flag'])
    env = None
    if int(runtime['reproducible_flag']):
        env = cmig_sing_utils.reproducible_env()
    if compiled_flag == 0:
        # the module's parameters are parsed from strings only when deployed,
        #   so the call is written as a MATLAB script with typed values and
        #   run with run_matlab_script; through cmig_run_bids_module.m, as the
        #   compiled app is, so the same thread count applies
        app = 'matlab'
        # no app runs, so the status gets n/a for app_version
        app_version = None
        margs = [typed_value(kwargs[key], CELL) if is_list(key) else kwargs[key]
                 for key in spec['positional']]
        parms = {key: typed_value(kwargs[key], kind)
                 for key, kind in spec['options'].items()
                 if kwargs.get(key) is not None}
        fd, fname_m = tempfile.mkstemp(prefix='cmig_run_bids_module_', suffix='.m')
        os.close(fd)
        try:
            # CMIG_MATLAB_ADDPATH (dirs separated by ':'): put working copies
            #   ahead of startup's addpath calls, which precede MATLABPATH
            dirs = [d for d in os.environ.get('CMIG_MATLAB_ADDPATH', '').split(':') if d]
            # addpath prepends, so the first dir listed is added last
            header = '\n'.join("addpath('{0}');".format(d) for d in reversed(dirs)) or None
            cmig_write_scripts.cmig_write_matlab_script(
                fname_m, DISPATCH_APPNAME, [module] + margs, tags=list(parms),
                parms=parms, header=header)
            subprocess.run(['run_matlab_script', fname_m], env=env)
        finally:
            os.remove(fname_m)
    else:
        app = DISPATCH_APPNAME
        app_version = appver
        cmig_sing_utils.check_compiled(runtime['appdir'], app, appver)
        cmig_sing_utils.run_compiled(runtime['appdir'], app, runtime['mrtdir'],
                                     [module] + container_args, appver, env=env)

    fname = status_file(module, kwargs)
    if fname is None:
        return check_concat_output(kwargs, app_name=app,
                                   app_version=app_version)
    return cmig_sing_utils.check_status(fname, module, app_name=app,
                                        app_version=app_version)


def build_parser(name, argv=()):
    """Return the parser of module name's arguments; argv, the arguments,
    only names the compiled app whose option help it shows. Raises
    ValueError for an unknown name."""
    module, fixed = resolve(name)
    spec = MODULES[module]
    # the help is read from the compiled app the runtime arguments name
    pre = argparse.ArgumentParser(add_help=False)
    for key in ['appdir', 'appver']:
        pre.add_argument(f'--{key}', default=RUNTIME_DEFAULTS[key])
    known, _ = pre.parse_known_args(list(argv))
    parser = argparse.ArgumentParser(
        prog=f'{os.path.basename(sys.argv[0])} {name}',
        usage='%(prog)s --<input> value ... [--<option> value ...] [-h]',
        description=f'Run {module}. Options not given take the module\'s '
                    f'own defaults; see its MATLAB help for their meaning.')
    inputs = parser.add_argument_group('inputs')
    optional = spec.get('optional', [])
    multi = spec.get('multi', [])
    for key in spec['positional']:
        if key in multi:
            inputs.add_argument(f'--{key}', required=key not in optional, nargs='+',
                                help='one or more')
        else:
            inputs.add_argument(f'--{key}', required=key not in optional)
    options = parser.add_argument_group('options (from the MATLAB help)')
    doc = load_help(known.appdir, known.appver).get(module, {}).get('options', {})
    for key, kind in spec['options'].items():
        if key in fixed:
            continue
        entry = doc.get(key)
        text = None
        if entry:
            # argparse expands %, so a literal one is doubled
            text = entry['description'].replace('%', '%%')
            if entry.get('default'):
                text += ' (default: ' + entry['default'].replace('%', '%%') + ')'
        if kind in (NUM, CELL):
            options.add_argument(f'--{key}', nargs='+', help=text)
        else:
            options.add_argument(f'--{key}', help=text,
                                 required=key in spec.get('required', []))
    runtime = parser.add_argument_group('runtime')
    runtime.add_argument('--compiled_flag', type=int, choices=[0, 1],
                         default=RUNTIME_DEFAULTS['compiled_flag'],
                         help='0 = MATLAB code, 1 = compiled app '
                              '(default: %(default)s)')
    runtime.add_argument('--reproducible_flag', type=int, choices=[0, 1],
                         default=RUNTIME_DEFAULTS['reproducible_flag'],
                         help=REPRODUCIBLE_HELP)
    for key in ['appver', 'appdir', 'pubsdir', 'mrtdir', 'mrtdir_target']:
        runtime.add_argument(f'--{key}', default=RUNTIME_DEFAULTS[key],
                             help='(default: %(default)s)')
    return parser


def parsed_kwargs(name, args):
    """Return the run_module kwargs of module name's parsed arguments."""
    kwargs = vars(args)
    # a single value of a multi input is passed as before, not as a cell array
    for key in MODULES[resolve(name)[0]].get('multi', []):
        val = kwargs.get(key)
        if isinstance(val, list) and len(val) == 1:
            kwargs[key] = val[0]
    return kwargs


def parse_args(argv=None):
    """Parse '<module> [arguments]'; return (module name, kwargs)."""
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ('-h', '--help'):
        print(__doc__.split('Created:')[0].rstrip())
        print('\nmodules:\n  ' + '\n  '.join(list_names()))
        sys.exit(0 if argv else 1)
    name = argv[0]
    try:
        parser = build_parser(name, argv[1:])
    except ValueError as e:
        sys.exit(f'Error: {e}')
    return name, parsed_kwargs(name, parser.parse_args(argv[1:]))


def main():
    name, kwargs = parse_args()
    errcode, _ = run_module(name, **kwargs)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
