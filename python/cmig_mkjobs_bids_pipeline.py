#!/usr/bin/env python3
"""
Generate cluster batch job scripts for CMIG BIDS pipeline and module scripts.

Created:  05/01/2026 by Don Hagler
Prev mod: 05/27/2026 by Don Hagler
Prev mod: 06/11/2026 by Don Hagler
Prev mod: 06/17/2026 by Don Hagler
Prev mod: 06/22/2026 by Don Hagler
Prev mod: 06/24/2026 by Don Hagler
Prev mod: 07/24/2026 by Don Hagler
Prev mod: 08/10/2026 by Don Hagler
Prev mod: 08/13/2026 by Don Hagler
Prev mod: 08/27/2026 by Don Hagler
Prev mod: 09/29/2026 by Don Hagler
Prev mod: 10/02/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/02/2026 by Don Hagler
Prev mod: 11/03/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/09/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Prev mod: 11/13/2026 by Don Hagler
Last mod: 11/14/2026 by Don Hagler
"""

import argparse
import datetime
import glob
import os
import re
import shlex
import sys

import cmig_get_bids_data
import cmig_pipeline_bids
import cmig_run_bids_module
import cmig_sing_utils
import cmig_write_bids_atlas
import cmig_write_bids_dataset_description
import cmig_write_status
import cmig_select_bids_func_runs

DEFAULT_BATCHDIR = os.path.join(os.path.expanduser('~'), 'batchdirs')
DEFAULT_APPVER = 'mmps_dev'
DEFAULT_APPDIR = '/usr/pubsw/packages/MMPS/apps'
DEFAULT_PUBSDIR = '/usr/pubsw'
DEFAULT_MRTDIR = '/usr/pubsw/packages/matlab/R2023a'
DEFAULT_MRTDIR_TARGET = '/opt/matlabruntime/R2023a'
DEFAULT_PYDIR = os.path.dirname(os.path.abspath(__file__))
# the python the job scripts run, e.g. the csh/cmig_python3 wrapper (the
#   python singularity image) in place of the host's python3
DEFAULT_PYTHON_CMD = 'python3'
DEFAULT_FS_VERSION = 820
# tasks whose E-Prime runs the eprime-scan check covers
DELAY_TASKS = ('mid', 'sst', 'nback')
IMAGE_ONLY_MSG = '{0} needs the dwi tools, in the image only: give --image (opts image)'
# image jobs: the license's path in the container, the job options that are
#   host settings, and the run.py root option of each dir key
IMAGE_FS_LICENSE = '/opt/fs_license.txt'
IMAGE_HOST_OPTS = ('compiled_flag', 'pydir')
IMAGE_ROOTS = {'rawqc_rootdir': 'rawqc_rootdir', 'proc_rootdir': 'proc_rootdir',
               'fsurf_rootdir': 'fsurf_rootdir', 'analysis_rootdir': 'analysis_rootdir',
               'postqc_rootdir': 'postqc_rootdir', 'pc_rootdir': 'pc_rootdir',
               'manualqc_rootdir': 'manualqc_rootdir',
               'proc_anat_rootdir': 'proc_anat_rootdir'}


# ---------------------------------------------------------------------------
# Script configuration
# ---------------------------------------------------------------------------
# Each entry has:
#   scan_dir_arg        : key in dirs dict used to discover {pid}/{sid} pairs
#   required_subdirs    : modality subdirectories that must exist under {scan_dir}/{pid}/{sid}/
#                         before a job is created (e.g. ['anat'], ['anat', 'func'])
#   status_files        : ordered list of (dir_key, subpath_template) — one per module;
#                         for pipelines, in execution order (used to detect prior failures)
#   dir_args            : dir keys passed as CLI flags to the script
#   optional_dir_args   : dir keys passed only when given, the script applying
#                         its own default otherwise
#   optional_output_dir_args : those of optional_dir_args the script writes to
#                         (created and described when given)
#   output_rootdir_from : (single-module scripts only) dir key whose value is passed
#                         as --output_rootdir to the script
#   extra_args          : non-dir, non-infra CLI flags the script accepts
#   image_arg_names     : optional {extra_args key: the script's own option name}
#   image_only          : if 1, the script runs only in the image (--image)
#   mrt_flag            : if 1, include appver/appdir/pubsdir/mrtdir/mrtdir_target
#   events_status       : optional (dir_key, subpath with {mod} and {task}): the
#                         status files are those of the tasks with events in rawdata
#                         func/ or beh/, one per such directory ({mod}), a group
#                         beside status_files; with no status_files (the
#                         standalone behavioral step) a session needs a task
#                         events file there, not required_subdirs, and the job
#                         gets --task and --ignore_delay_flag for its tasks.
#                         A status is stale where its events changed after it
#                         ran, or (func/) its recorded ignore_delay_flag
#                         differs from the REDCap flags'
#   events_newer        : optional (dir_key, status subpath): stale where any of
#                         the session's rawdata events changed after that
#                         status's run started
#   no_delete           : if 1, errflag deletes nothing (the script reruns its failures)
#   modality_status     : optional (dir_key, subpath with {mod}, [mods]): one status
#                         file per modality whose rawdata {pid}/{sid}/{mod} exists,
#                         ahead of status_files, so errflag 1 deletes only the
#                         failed modality's directory
#   stale_check         : optional {'output': (dir_key, subpath), 'inputs':[]}.
#                         A session is stale where an input ran after the
#                         output did, comparing each input's status procdate
#                         against the output's status startdate. A dir_key
#                         absent from dirs falls back to the output's, as the
#                         script's own default does
#   stale_steps         : optional list of STALE_STEPS names the script runs
#                         ('pipeline': those cmig_pipeline_bids runs for the
#                         session); each is stale where an input step ran
#                         after it started, and the job deletes its outputs

SCRIPT_CONFIG = {
    'cmig_proc_bids_anat': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('proc_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-proc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'proc_rootdir'],
        'output_rootdir_from': 'proc_rootdir',
        'extra_args': ['gradunwarp_flag', 'cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['proc_anat'],
    },
    'cmig_proc_bids_func': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['func'],
        'status_files': [
            ('proc_rootdir', '{pid}/{sid}/func/{pid}_{sid}_desc-proc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'proc_rootdir'],
        'output_rootdir_from': 'proc_rootdir',
        'extra_args': ['gradunwarp_flag', 'cleanup_flag', 'std_orient', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['proc_func'],
    },
    # the dwi tools are in the image only
    'cmig_proc_bids_dwi': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['dwi'],
        'status_files': [
            ('proc_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-proc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'proc_rootdir'],
        'optional_dir_args': ['proc_anat_rootdir'],
        'output_rootdir_from': 'proc_rootdir',
        'extra_args': ['dwi_n_cpus', 'cleanup_flag'],
        'image_arg_names': {'dwi_n_cpus': 'n_cpus'},
        'mrt_flag': 0,
        'image_only': 1,
        'stale_steps': ['proc_dwi'],
    },
    'cmig_fsurf_bids': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('fsurf_rootdir', '{pid}_{sid}/{pid}_{sid}_desc-fsurf_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir'],
        'output_rootdir_from': 'fsurf_rootdir',
        'extra_args': ['T2_flag', 'FS_version'],
        'mrt_flag': 0,
        'stale_steps': ['fsurf'],
    },
    'cmig_analyze_bids_anat': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('analysis_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-analyze_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'analysis_rootdir'],
        'output_rootdir_from': 'analysis_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag'],
        'mrt_flag': 1,
        'stale_steps': ['analyze_anat'],
    },
    'cmig_analyze_bids_dwi': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['dwi'],
        'status_files': [
            ('analysis_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-analyze_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'analysis_rootdir'],
        'output_rootdir_from': 'analysis_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['analyze_dwi'],
    },
    'cmig_analyze_bids_func': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['func'],
        'status_files': [
            ('analysis_rootdir', '{pid}/{sid}/func/{pid}_{sid}_task-{task}_desc-analyze_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'analysis_rootdir'],
        'output_rootdir_from': 'analysis_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag', 'task', 'stim_times_model',
                       'stim_times_model_mid', 'stim_times_model_sst',
                       'stim_times_model_nback', 'stim_times_flag',
                       'baseline_flag', 'skipTRs', 'despike_flag',
                       'ignore_delay_flag', 'exclude_noncompliant_flag', 'pydir'],
        'mrt_flag': 1,
        'stale_steps': ['analyze_func'],
    },
    'cmig_pipeline_bids': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('proc_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-proc_status.tsv'),
            ('fsurf_rootdir', '{pid}_{sid}/{pid}_{sid}_desc-fsurf_status.tsv'),
            ('analysis_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-analyze_status.tsv'),
            ('postqc_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-postqc_status.tsv'),
            ('postqc_rootdir', '{pid}/{sid}/fsurf/{pid}_{sid}_desc-postqc_status.tsv'),
        ],
        # expected besides for a session with task fMRI data, unless func_flag 0
        'func_status_files': [
            ('proc_rootdir', '{pid}/{sid}/func/{pid}_{sid}_desc-proc_status.tsv'),
            ('analysis_rootdir', '{pid}/{sid}/func/{pid}_{sid}_task-{task}_desc-analyze_status.tsv'),
            ('postqc_rootdir', '{pid}/{sid}/func/{pid}_{sid}_desc-postqc_status.tsv'),
            ('postqc_rootdir', '{pid}/{sid}/fmap/{pid}_{sid}_acq-func_desc-postqc_status.tsv'),
        ],
        # expected besides for a session with dwi data in rawdata, unless
        #   dwi_flag 0; analyze_dwi only with a complete recon
        'dwi_status_files': [
            ('proc_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-proc_status.tsv'),
            ('analysis_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-analyze_status.tsv'),
            ('postqc_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-postqc_status.tsv'),
            ('postqc_rootdir', '{pid}/{sid}/fmap/{pid}_{sid}_acq-dwi_desc-postqc_status.tsv'),
        ],
        'events_status': ('analysis_rootdir',
                          '{pid}/{sid}/{mod}/{pid}_{sid}_task-{task}_desc-behav_status.tsv'),
        'dir_args': ['raw_rootdir', 'proc_rootdir', 'fsurf_rootdir', 'analysis_rootdir', 'postqc_rootdir'],
        'extra_args': ['gradunwarp_flag', 'cleanup_flag', 'std_orient', 'func_flag', 'dwi_flag', 'dwi_n_cpus', 'analyze_anat_flag', 'task', 'compiled_flag', 'reproducible_flag',
                       'ignore_delay_tasks', 'exclude_noncompliant_flag', 'FS_version', 'pydir'],
        'mrt_flag': 1,
        'stale_steps': 'pipeline',
    },
    'cmig_postqc_bids_fmap': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['fmap'],
        'status_files': [
            ('postqc_rootdir', '{pid}/{sid}/fmap/{pid}_{sid}_desc-postqc_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'postqc_rootdir'],
        'output_rootdir_from': 'postqc_rootdir',
        'extra_args': ['acq_types', 'cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['postqc_fmap_func', 'postqc_fmap_dwi'],
    },
    'cmig_postqc_bids_anat': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('postqc_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-postqc_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'postqc_rootdir'],
        'output_rootdir_from': 'postqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['postqc_anat'],
    },
    'cmig_postqc_bids_dwi': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['dwi'],
        'status_files': [
            ('postqc_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-postqc_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'postqc_rootdir'],
        'output_rootdir_from': 'postqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['postqc_dwi'],
    },
    'cmig_postqc_bids_func': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['func'],
        'status_files': [
            ('postqc_rootdir', '{pid}/{sid}/func/{pid}_{sid}_desc-postqc_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'postqc_rootdir'],
        'output_rootdir_from': 'postqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['postqc_func'],
    },
    'cmig_postqc_bids_fsurf': {
        'scan_dir_arg': 'proc_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('postqc_rootdir', '{pid}/{sid}/fsurf/{pid}_{sid}_desc-postqc_status.tsv'),
        ],
        'dir_args': ['proc_rootdir', 'fsurf_rootdir', 'postqc_rootdir'],
        'output_rootdir_from': 'postqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
        'stale_steps': ['postqc_fsurf'],
    },
    'cmig_rawqc_bids_anat': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('rawqc_rootdir', '{pid}/{sid}/anat/{pid}_{sid}_desc-rawqc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'rawqc_rootdir'],
        'output_rootdir_from': 'rawqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
    },
    'cmig_rawqc_bids_dwi': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['dwi'],
        'status_files': [
            ('rawqc_rootdir', '{pid}/{sid}/dwi/{pid}_{sid}_desc-rawqc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'rawqc_rootdir'],
        'output_rootdir_from': 'rawqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
    },
    'cmig_rawqc_bids_fmap': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['fmap'],
        'status_files': [
            ('rawqc_rootdir', '{pid}/{sid}/fmap/{pid}_{sid}_desc-rawqc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'rawqc_rootdir'],
        'output_rootdir_from': 'rawqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
    },
    'cmig_rawqc_bids_func': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['func'],
        'status_files': [
            ('rawqc_rootdir', '{pid}/{sid}/func/{pid}_{sid}_desc-rawqc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'rawqc_rootdir'],
        'output_rootdir_from': 'rawqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag'],
        'mrt_flag': 1,
    },
    'cmig_pipeline_bids_rawqc': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['anat'],
        # the modalities the pipeline runs (those in rawdata), then its own
        #   session status, last, for pc and pcqc and a job that did not finish
        'modality_status': ('rawqc_rootdir',
                            '{pid}/{sid}/{mod}/{pid}_{sid}_desc-rawqc_status.tsv',
                            ['anat', 'func', 'dwi', 'fmap']),
        'status_files': [
            ('rawqc_rootdir', '{pid}/{sid}/{pid}_{sid}_desc-rawqc_status.tsv'),
        ],
        'dir_args': ['raw_rootdir', 'rawqc_rootdir'],
        'optional_dir_args': ['manualqc_rootdir', 'pc_rootdir'],
        'optional_output_dir_args': ['pc_rootdir'],
        'output_rootdir_from': 'rawqc_rootdir',
        'extra_args': ['cleanup_flag', 'compiled_flag', 'reproducible_flag', 'pc_flag',
                       'exclude_noncompliant_flag'],
        'mrt_flag': 1,
        # the merge is the only step whose inputs are other modules' output, so
        # it is the only one this notices. The rest waits for Airflow
        'stale_check': {
            'output': ('rawqc_rootdir',
                       '{pid}/{sid}/{pid}_{sid}_desc-pcqc_status.tsv'),
            'inputs': [
                ('pc_rootdir',
                 '{pid}/{sid}/{pid}_{sid}_desc-pc_status.tsv'),
                ('rawqc_rootdir',
                 '{pid}/{sid}/*/{pid}_{sid}_desc-rawqc_status.tsv'),
                ('manualqc_rootdir',
                 '{pid}/{sid}/{pid}_{sid}_desc-manualrawqc_status.tsv'),
            ],
        },
    },
    'cmig_analyze_bids_behav': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': [],
        'events_status': ('analysis_rootdir',
                          '{pid}/{sid}/{mod}/{pid}_{sid}_task-{task}_desc-behav_status.tsv'),
        'status_files': [],
        'dir_args': ['raw_rootdir', 'analysis_rootdir'],
        'output_rootdir_from': 'analysis_rootdir',
        'extra_args': ['task', 'ignore_delay_flag', 'exclude_noncompliant_flag'],
        'mrt_flag': 0,
        'no_delete': 1,
    },
    'cmig_pc_bids_abcd': {
        'scan_dir_arg': 'raw_rootdir',
        'required_subdirs': ['anat'],
        'status_files': [
            ('pc_rootdir', '{pid}/{sid}/{pid}_{sid}_desc-pc_status.tsv'),
        ],
        # its E-Prime session variables are read from the rawdata events
        'events_newer': ('pc_rootdir', '{pid}/{sid}/{pid}_{sid}_desc-pc_status.tsv'),
        'dir_args': ['raw_rootdir', 'pc_rootdir'],
        'output_rootdir_from': 'pc_rootdir',
        'extra_args': ['exclude_noncompliant_flag'],
        'mrt_flag': 0,
    },
}

# ---------------------------------------------------------------------------
# Stale steps
# ---------------------------------------------------------------------------
# The steps of cmig_pipeline_bids in run order: each one's status file, the
# steps whose outputs it reads (inputs), and the outputs a job deletes when it
# is stale: (dir_key, directory, file pattern or None for the whole directory).
# {p} is {pid}/{sid}, {s} is {pid}_{sid}
_STEP_TABLE = [
    ('proc_anat', ('proc_rootdir', '{p}/anat/{s}_desc-proc_status.tsv'), [],
     [('proc_rootdir', '{p}/anat', None)]),
    ('proc_func', ('proc_rootdir', '{p}/func/{s}_desc-proc_status.tsv'),
     ['proc_anat'],
     [('proc_rootdir', '{p}/func', None),
      ('proc_rootdir', '{p}/fmap', '{s}_acq-func_*'),
      ('proc_rootdir', '{p}/fmap/proc_fmap', None)]),
    ('proc_dwi', ('proc_rootdir', '{p}/dwi/{s}_desc-proc_status.tsv'),
     ['proc_anat'],
     [('proc_rootdir', '{p}/dwi', None),
      ('proc_rootdir', '{p}/fmap', '{s}_acq-dwi_*')]),
    ('fsurf', ('fsurf_rootdir', '{s}/{s}_desc-fsurf_status.tsv'),
     ['proc_anat'],
     [('fsurf_rootdir', '{s}', None)]),
    ('analyze_anat', ('analysis_rootdir', '{p}/anat/{s}_desc-analyze_status.tsv'),
     ['proc_anat', 'fsurf'],
     [('analysis_rootdir', '{p}/anat', None)]),
    # the deletion is set by _analyze_func_deletions()
    ('analyze_func', ('analysis_rootdir',
                      '{p}/func/{s}_task-{task}_desc-analyze_status.tsv'),
     ['proc_func', 'fsurf'], []),
    ('analyze_dwi', ('analysis_rootdir', '{p}/dwi/{s}_desc-analyze_status.tsv'),
     ['proc_dwi', 'fsurf'],
     [('analysis_rootdir', '{p}/dwi', None)]),
    ('postqc_anat', ('postqc_rootdir', '{p}/anat/{s}_desc-postqc_status.tsv'),
     ['proc_anat', 'fsurf'],
     [('postqc_rootdir', '{p}/anat', None)]),
    ('postqc_fsurf', ('postqc_rootdir', '{p}/fsurf/{s}_desc-postqc_status.tsv'),
     ['fsurf', 'proc_anat'],
     [('postqc_rootdir', '{p}/fsurf', None)]),
    ('postqc_func', ('postqc_rootdir', '{p}/func/{s}_desc-postqc_status.tsv'),
     ['proc_func', 'proc_anat', 'fsurf'],
     [('postqc_rootdir', '{p}/func', None)]),
    ('postqc_fmap_func', ('postqc_rootdir',
                          '{p}/fmap/{s}_acq-func_desc-postqc_status.tsv'),
     ['proc_func', 'proc_anat', 'fsurf'],
     [('postqc_rootdir', '{p}/fmap', '{s}_acq-func_*')]),
    ('postqc_dwi', ('postqc_rootdir', '{p}/dwi/{s}_desc-postqc_status.tsv'),
     ['proc_dwi', 'proc_anat', 'fsurf'],
     [('postqc_rootdir', '{p}/dwi', None)]),
    ('postqc_fmap_dwi', ('postqc_rootdir',
                         '{p}/fmap/{s}_acq-dwi_desc-postqc_status.tsv'),
     ['proc_dwi', 'proc_anat', 'fsurf'],
     [('postqc_rootdir', '{p}/fmap', '{s}_acq-dwi_*')]),
]
STALE_STEPS = dict((name, {'status': status, 'inputs': inputs, 'delete': delete})
                   for name, status, inputs, delete in _STEP_TABLE)
STEP_ORDER = [row[0] for row in _STEP_TABLE]
# the analysis func/ files of the behavioral step, kept when only one task's
#   analysis is deleted, since analyze_func reads its run selection
BEHAV_SUFFIXES = ('_behavioral.tsv', '_desc-select_runs.tsv', '_desc-behav_status.tsv')
STALE_COLUMNS = ['participant_id', 'session_id', 'step', 'task', 'reason',
                 'input', 'input_time', 'startdate', 'deleted']


# ---------------------------------------------------------------------------
# Session discovery and status checking
# ---------------------------------------------------------------------------

def _find_sessions(scan_dir, sublist=None):
    sessions = []
    if sublist is not None:
        # only the named participants are read, so a scan_dir holding a whole
        # release costs len(sublist) directory reads rather than all of them;
        # sub- is added where missing, and a repeat is listed once
        entries = sorted(set(p if p.startswith('sub-') else 'sub-' + p
                             for p in sublist))
    else:
        try:
            entries = sorted(os.listdir(scan_dir))
        except OSError as e:
            sys.exit(f'Error reading scan directory {scan_dir}: {e}')
    for pid in entries:
        if not pid.startswith('sub-'):
            continue
        pid_dir = os.path.join(scan_dir, pid)
        if not os.path.isdir(pid_dir):
            continue
        try:
            sids = sorted(os.listdir(pid_dir))
        except OSError:
            continue
        for sid in sids:
            if not sid.startswith('ses-'):
                continue
            if os.path.isdir(os.path.join(pid_dir, sid)):
                sessions.append((pid, sid))
    return sessions


def _session_subdir(pid, sid, subpath_template):
    """Return the session-level subdirectory from a status file path template.

    Templates with '{pid}/{sid}/...' return '{pid}/{sid}'.
    Templates with '{pid}_{sid}/...' return '{pid}_{sid}'.
    """
    if subpath_template.split('/')[0] == '{pid}':
        return f'{pid}/{sid}'
    return f'{pid}_{sid}'


def _module_subdir(pid, sid, subpath_template):
    """Return the failing module's subdirectory, which errflag 1 deletes
    where the failed step is not in the step table (_failed_step_deletions).

    '{pid}/{sid}/<mod>/<status file>' gives '{pid}/{sid}/<mod>' (func, anat,
    dwi, fmap, beh, fsurf), so one modality's failure leaves the session's
    other modalities alone; a session-level status ('{pid}/{sid}/<file>',
    '{pid}_{sid}/<file>') gives the session directory, as _session_subdir().
    """
    parts = subpath_template.split('/')
    if parts[0] == '{pid}' and len(parts) >= 4:
        return '/'.join([pid, sid] + [p.format(pid=pid, sid=sid) for p in parts[2:-1]])
    return _session_subdir(pid, sid, subpath_template)


def _func_included(config, tasks, opts):
    """True where a config's func_status_files apply: the session has task
    fMRI data (tasks) and func_flag is not 0."""
    if not config.get('func_status_files'):
        return False
    return bool(tasks) and (opts or {}).get('func_flag', 1) != 0


def _dwi_included(config, opts, pid=None, sid=None, dirs=None):
    """True where a config's dwi_status_files apply: the session has dwi data
    in rawdata and dwi_flag is not 0."""
    if not config.get('dwi_status_files'):
        return False
    if (opts or {}).get('dwi_flag', 1) == 0:
        return False
    return bool(dirs and dirs.get('raw_rootdir') and
                cmig_pipeline_bids.has_dwi(pid, sid, dirs['raw_rootdir']))


def _status_files(config, tasks, opts, pid=None, sid=None, dirs=None):
    """The config's status files for one session, func ones included where
    _func_included() and dwi ones where _dwi_included(); without a postqc fmap
    one where cmig_pipeline_bids skips that step (no processed and no usable
    field map of that acq type); first, the modality_status ones of the
    modalities in rawdata."""
    files = []
    if config.get('modality_status') and dirs and dirs.get('raw_rootdir'):
        dir_key, subpath_tmpl, mods = config['modality_status']
        files += [(dir_key, subpath_tmpl.replace('{mod}', mod)) for mod in mods
                  if os.path.isdir(os.path.join(dirs['raw_rootdir'], pid, sid, mod))]
    files += list(config['status_files'])
    have_dirs = bool(dirs and dirs.get('proc_rootdir') and dirs.get('raw_rootdir'))
    groups = []
    if _func_included(config, tasks, opts):
        groups.append(('func', config['func_status_files']))
    if _dwi_included(config, opts, pid, sid, dirs):
        groups.append(('dwi', config['dwi_status_files']))
    for acq, group in groups:
        if have_dirs and cmig_pipeline_bids.skip_postqc_fmap(
                pid, sid, dirs['proc_rootdir'], dirs['raw_rootdir'], acq_type=acq):
            group = [f for f in group if '/fmap/' not in f[1]]
        files += group
    return files


def _session_tasks(pid, sid, dirs, opts):
    """Return the task labels to expect status files for, for one session.

    An explicit --task wins. Otherwise the list is discovered from raw_rootdir
    exactly as cmig_pipeline_bids discovers it, so both sides agree without
    either passing the list to the other.
    """
    task = opts.get('task', None) if opts else None
    if isinstance(task, list):
        tasks = [t for t in task if t]
    elif task:
        tasks = [task]
    else:
        tasks = []
    if tasks:
        return tasks
    raw_rootdir = dirs.get('raw_rootdir')
    if not raw_rootdir:
        return []
    return cmig_get_bids_data.get_bids_tasks(pid, sid, raw_rootdir)


def _expand_status_files(status_files, tasks):
    """Expand task-bearing status file templates into groups of entries.

    Returns a list of groups, each a list of (dir_key, subpath_template, task).
    A template without a {task} placeholder yields a group of one. A
    task-bearing template yields one entry per task, and every entry must have
    succeeded for the group to be satisfied, so a session missing one task's
    output still gets a job.

    Single-task scripts such as cmig_analyze_bids_func are passed one task, so
    the group holds that one entry.
    """
    groups = []
    for dir_key, subpath_tmpl in status_files:
        if '{task}' in subpath_tmpl:
            groups.append([(dir_key, subpath_tmpl, tsk) for tsk in tasks])
        else:
            groups.append([(dir_key, subpath_tmpl, '')])
    return groups


def load_ignore_delay_flags(fname_csv=None, flags_dir=None):
    """{(participant_id, session_id): {task: 0/1}} from the REDCap tracking table,
    or, with flags_dir, from release-layout tables.

    eprime_mismatch_ok_{mid,sst,nback} hold both the automatic session-level
    exemption and staff overrides; a row set for a task in any of a session's
    rows sets it. Keys are BIDS, mapped as cmig_convert_abcd_eprime.py maps
    pGUID and EventName. flags_dir holds mr_y_qc__raw__tfmri__{task}.tsv with
    the released mr_y_qc__raw__tfmri__{task}__eprime__tdiff__ign_indicator,
    the release's own tables or those cmig_write_eprime_mismatch_tables.py
    writes from the tracking table.
    """
    import csv
    flags = {}
    if flags_dir is not None:
        if not os.path.isdir(flags_dir):
            return None
        for task in DELAY_TASKS:
            fname = os.path.join(flags_dir, 'mr_y_qc__raw__tfmri__%s.tsv' % task)
            col = 'mr_y_qc__raw__tfmri__%s__eprime__tdiff__ign_indicator' % task
            if not os.path.isfile(fname):
                return None
            with open(fname) as fid:
                for row in csv.DictReader(fid, delimiter='\t'):
                    if (row.get(col) or '').strip() == '1':
                        flags.setdefault((row['participant_id'], row['session_id']), {})[task] = 1
        return flags
    from cmig_convert_abcd_eprime import bids_sub, bids_ses
    if str(fname_csv).lower() == 'none':
        print('NOTE: --tracking_csv none -- ignore_delay_flag 0 for every session')
        return flags
    if not fname_csv or not os.path.isfile(fname_csv):
        # without it every flagged run would be excluded, with no error
        return None
    with open(fname_csv) as fid:
        for row in csv.DictReader(fid):
            try:
                key = (bids_sub(row['id_redcap']), bids_ses(row['redcap_event_name']))
            except (KeyError, ValueError):
                continue
            entry = flags.setdefault(key, {})
            for task in DELAY_TASKS:
                val = (row.get('eprime_mismatch_ok_' + task) or '').strip()
                if val == '1':
                    entry[task] = 1
    return flags


def get_ignore_delay_flag(flags, pid, sid, task):
    """The REDCap ignore_delay_flag for one BIDS participant, session and task.

    The one place the flag is looked up, so a precomputed eprime QC file in a
    BIDS derivatives directory can replace the tracking table later.
    """
    return flags.get((pid, sid), {}).get(task.lower(), 0)


def _events_tasks(pid, sid, raw_rootdir, tasks=None):
    """{task: [modality, ...]}: the rawdata func/ and beh/ holding its events,
    for the E-Prime tasks the behavioral step runs (events_tasks())."""
    found = cmig_select_bids_func_runs.events_tasks(raw_rootdir, pid, sid, tasks)
    return dict((t, found[t]) for t in sorted(found))


def _check_stale(pid, sid, dirs, config):
    """The input that ran after the output that read it (cmig_write_status.
    inputs_newer: input path, its time, the output's time), or None/False.

    Each input is another module's output, so its own status file records when
    it ran; cmig_write_status.inputs_newer makes the comparison, and
    cmig_check_event_pcqc makes the same one on itself before it skips.
    """
    check = config.get('stale_check')
    if not check:
        return False
    out_key, out_tmpl = check['output']
    out_root = dirs.get(out_key)
    if not out_root:
        return False
    status_path = os.path.join(out_root, out_tmpl.format(pid=pid, sid=sid))
    input_paths = []
    for dir_key, subpath_tmpl in check['inputs']:
        rootdir = dirs.get(dir_key) or out_root
        pattern = os.path.join(rootdir, subpath_tmpl.format(pid=pid, sid=sid))
        input_paths.extend(sorted(glob.glob(pattern)))
    return cmig_write_status.inputs_newer(status_path, input_paths)


def _step_path(name, pid, sid, dirs, task=''):
    """The status file of a STALE_STEPS step, or None without its root."""
    dir_key, tmpl = STALE_STEPS[name]['status']
    if not dirs.get(dir_key):
        return None
    return os.path.join(dirs[dir_key], tmpl.format(
        p=f'{pid}/{sid}', s=f'{pid}_{sid}', task=task))


def _run_steps(pid, sid, dirs, config, opts, tasks):
    """The (step, task) pairs of STALE_STEPS a session's job runs, in run
    order; [] for a script without stale_steps."""
    steps = config.get('stale_steps')
    opts = opts or {}
    if not steps:
        return []
    if steps == 'pipeline':
        names = ['proc_anat', 'fsurf', 'postqc_anat', 'postqc_fsurf']
        if opts.get('analyze_anat_flag', 1) != 0:
            names.append('analyze_anat')
        if _func_included(config, tasks, opts):
            names += ['proc_func', 'analyze_func', 'postqc_func', 'postqc_fmap_func']
        if _dwi_included(config, opts, pid, sid, dirs):
            names += ['proc_dwi', 'analyze_dwi', 'postqc_dwi', 'postqc_fmap_dwi']
    else:
        names = list(steps)
        acq_types = opts.get('acq_types')
        if acq_types:
            names = [n for n in names if not n.startswith('postqc_fmap_') or
                     n[len('postqc_fmap_'):] in acq_types]
    run = []
    for name in sorted(names, key=STEP_ORDER.index):
        if name == 'analyze_func':
            run += [(name, t) for t in tasks]
        else:
            run.append((name, ''))
    return run


def _scans_path(raw_rootdir, pid, sid):
    return os.path.join(raw_rootdir, pid, sid, f'{pid}_{sid}_scans.tsv')


def _read_scans(raw_rootdir, pid, sid):
    """The rows of a session's rawdata scans.tsv, or []."""
    import csv
    fname = _scans_path(raw_rootdir, pid, sid)
    if not os.path.isfile(fname):
        return []
    with open(fname) as fid:
        return list(csv.DictReader(fid, delimiter='\t'))


def _sidecar_uid(fname_json):
    import json
    try:
        with open(fname_json) as fid:
            return json.load(fid).get('SeriesInstanceUID', '')
    except (IOError, OSError, ValueError):
        return ''


def _raw_pair(raw_rootdir, pid, sid, filename):
    """(file stem, SeriesInstanceUID) of a scans.tsv filename."""
    stem = os.path.basename(filename).split('.nii')[0]
    return stem, _sidecar_uid(os.path.join(
        raw_rootdir, pid, sid, filename.split('.nii')[0] + '.json'))


def _proc_pairs(proc_rootdir, pid, sid, subdir, suffixes, prefix=''):
    """(file stem, SeriesInstanceUID) of the proc output sidecars in subdir
    whose names end with one of suffixes (and start with pid_sid_prefix).
    The stem drops desc-preproc, so it compares with the rawdata stem; a
    sidecar without it (written before desc-preproc) is read too."""
    pairs = set()
    for fname in glob.glob(os.path.join(proc_rootdir, pid, sid, subdir, '*.json')):
        name = os.path.basename(fname)
        if not name.startswith(f'{pid}_{sid}_{prefix}'):
            continue
        if not any(name.endswith(x) for x in suffixes):
            continue
        name = re.sub(r'_desc-preproc(_[^_]+\.json)$', r'\1', name)
        if '_desc-' in name or '_space-' in name:
            continue
        pairs.add((name[:-len('.json')], _sidecar_uid(fname)))
    return pairs


def _usable(rows, subdir, suffix, contains=''):
    return [r for r in rows if (r.get('usable') or '').strip() == '1'
            and r.get('filename', '').startswith(subdir + '/')
            and contains in r['filename'] and suffix + '.nii' in r['filename']]


def _freq_steps(info):
    """FrequencyEncodingSteps as cmig_proc_bids_func.m load_scan_info reads
    it: BaseResolution on Siemens, FrequencyEncodingSteps on Philips, else
    FrequencyEncodingSteps or AcquisitionMatrixPE."""
    manuf = (info.get('Manufacturer') or '').upper()
    if 'SIEMENS' in manuf:
        return info.get('BaseResolution')
    if 'PHILIPS' in manuf:
        return info.get('FrequencyEncodingSteps')
    val = info.get('FrequencyEncodingSteps')
    return val if val is not None else info.get('AcquisitionMatrixPE')


def _fmap_pair_ok(fwd, rev):
    """True where cmig_proc_bids_func.m check_fmap_pairs keeps a forward and
    reverse field map, by their json sidecars: phase-encode axis, sizes, echo
    spacing (relative tolerance 1e-4), orientation and shims. A scalar
    missing on either side passes; a vector missing on one side does not."""
    if (fwd.get('PhaseEncodingDirection') or 'U')[0] != \
            (rev.get('PhaseEncodingDirection') or 'U')[0]:
        return False
    for get in (lambda i: i.get('SliceThickness'), lambda i: i.get('ReconMatrixPE'),
                _freq_steps, lambda i: i.get('PhaseEncodingSteps')):
        a, b = get(fwd), get(rev)
        if a is not None and b is not None and a != b:
            return False
    a, b = fwd.get('EffectiveEchoSpacing'), rev.get('EffectiveEchoSpacing')
    if a is not None and b is not None and abs(a - b) > 1e-4 * max(abs(a), abs(b)):
        return False
    for key in ('ImageOrientationPatientDICOM', 'ShimSetting'):
        a, b = fwd.get(key), rev.get(key)
        if (a is None) != (b is None) or (a is not None and a != b):
            return False
    return True


def _paired_fmaps(raw_rootdir, pid, sid, rows):
    """The usable field map rows cmig_proc_bids_func.m pairs: one forward and
    one reverse (PhaseEncodingDirection with or without -) per
    B0FieldIdentifier, kept by the json checks of check_fmap_pairs
    (_fmap_pair_ok). Its NIfTI header checks (vox2ras, dimensions) are not
    made here, so a pair they reject reads as a change, costing one rerun."""
    import json
    groups = {}
    for r in rows:
        try:
            with open(os.path.join(raw_rootdir, pid, sid,
                                   r['filename'].split('.nii')[0] + '.json')) as fid:
                info = json.load(fid)
        except (IOError, OSError, ValueError):
            continue
        groups.setdefault(info.get('B0FieldIdentifier') or 'UNKNOWN', []).append((r, info))
    paired = []
    for members in groups.values():
        rev = [m for m in members if '-' in (m[1].get('PhaseEncodingDirection') or '')]
        fwd = [m for m in members if '-' not in (m[1].get('PhaseEncodingDirection') or '')]
        if len(rev) == 1 and len(fwd) == 1 and _fmap_pair_ok(fwd[0][1], rev[0][1]):
            paired += [fwd[0][0], rev[0][0]]
    return paired


def _selection_changed(name, pid, sid, dirs):
    """Why a proc step's scans from the current rawdata scans.tsv differ from
    those its output sidecars record, or '': (file stem, SeriesInstanceUID)
    pairs (a renamed run is a change). Compared only where scans.tsv changed
    after the step started."""
    raw, proc = dirs.get('raw_rootdir'), dirs.get('proc_rootdir')
    fstatus = _step_path(name, pid, sid, dirs)
    started = cmig_write_status.status_time(fstatus, 'startdate')
    fscans = _scans_path(raw or '', pid, sid)
    if not raw or not proc or started is None or not os.path.isfile(fscans):
        return ''
    if datetime.datetime.fromtimestamp(os.path.getmtime(fscans)) <= started:
        return ''
    rows = _read_scans(raw, pid, sid)
    if name == 'proc_anat':
        from cmig_select_bids_dwi_runs import pick
        now = set()
        for suffix in ('_T1w', '_T2w'):
            cands = _usable(rows, 'anat', suffix)
            if cands:
                r = pick(cands, lambda r: r.get('notes') or '',
                         lambda r: r.get('acq_time') or '')
                now.add(_raw_pair(raw, pid, sid, r['filename']))
        then = _proc_pairs(proc, pid, sid, 'anat', ('_T1w.json', '_T2w.json'))
    elif name == 'proc_func':
        now = set(_raw_pair(raw, pid, sid, r['filename'])
                  for r in _usable(rows, 'func', '_bold') +
                  _paired_fmaps(raw, pid, sid, _usable(rows, 'fmap', '_epi', '_acq-func_')))
        then = (_proc_pairs(proc, pid, sid, 'func', ('_bold.json',)) |
                _proc_pairs(proc, pid, sid, 'fmap', ('_epi.json',), 'acq-func_'))
    elif name == 'proc_dwi':
        from cmig_select_bids_dwi_runs import select_runs, NotProcessable
        try:
            sel = select_runs(os.path.join(raw, pid, sid))
            now = set((os.path.basename(sel[k]).split('.nii')[0],
                       _sidecar_uid(sel[k].split('.nii')[0] + '.json'))
                      for k in ('up_nifti', 'up_nifti2', 'down_nifti') if sel.get(k))
        except (NotProcessable, OSError, ValueError, KeyError):
            now = set()
        then = (_proc_pairs(proc, pid, sid, 'dwi', ('_dwi.json',)) |
                _proc_pairs(proc, pid, sid, 'fmap', ('_epi.json',), 'acq-dwi_dir-AP_'))
    else:
        return ''
    if now == then:
        return ''
    added = set(n for n, _ in now - then)
    removed = set(n for n, _ in then - now)
    # a name in both, with another series under it: a renamed run
    moved = sorted(added & removed)
    added, removed = sorted(added - removed), sorted(removed - added)
    parts = [label + ' ' + ' '.join(names) for label, names in
             (('now', added), ('was', removed), ('other series', moved)) if names]
    return 'scans.tsv selection changed: ' + '; '.join(parts)


def _events_copy_newer(pid, sid, dirs, task, fstatus):
    """(events copy, its mtime, startdate) where a task's events copy in proc
    func/ was made after analyze_func started, else None: analyze_func reads
    the copies, which proc_func remakes when the rawdata events change."""
    started = cmig_write_status.status_time(fstatus, 'startdate')
    if started is None or not dirs.get('proc_rootdir'):
        return None
    for fname in sorted(glob.glob(os.path.join(
            dirs['proc_rootdir'], pid, sid, 'func',
            f'{pid}_{sid}_task-{task}_*events.tsv'))):
        mtime = datetime.datetime.fromtimestamp(os.path.getmtime(fname))
        if mtime > started:
            return fname, mtime, started
    return None


def _deletion_cmd(rootdir, subdir, pattern, exclude=()):
    """A csh line deleting a directory, or the files matching pattern in it
    (find, as a csh glob with no match stops the line)."""
    path = os.path.join(rootdir, subdir)
    if pattern is None:
        return f'rm -rf {path}'
    nots = ''.join(f" ! -name '*{x}'" for x in exclude)
    return (f"if (-d {path}) find {path} -maxdepth 1 -name '{pattern}'{nots}"
            " -exec rm -rf {} +")


def _step_deletions(name, task, pid, sid, dirs, config, opts):
    """The csh deletion lines of a stale step."""
    fmt = dict(p=f'{pid}/{sid}', s=f'{pid}_{sid}', task=task)
    if name == 'analyze_func':
        root = dirs.get('analysis_rootdir')
        if not root:
            return []
        # a pipeline that runs every task reruns the behavioral step and all
        #   tasks, so the whole func/ goes; otherwise only this task's
        #   analysis, leaving the behavioral step's run selection it reads
        if config.get('stale_steps') == 'pipeline' and not (opts or {}).get('task'):
            return [_deletion_cmd(root, f'{pid}/{sid}/func', None)]
        return [_deletion_cmd(root, f'{pid}/{sid}/func',
                              f'{pid}_{sid}_task-{task}_*', BEHAV_SUFFIXES),
                _deletion_cmd(root, f'{pid}/{sid}/func/analyze_func_{task}', None)]
    cmds = []
    for dir_key, subdir, pattern in STALE_STEPS[name]['delete']:
        if dirs.get(dir_key):
            cmds.append(_deletion_cmd(dirs[dir_key], subdir.format(**fmt),
                                      pattern.format(**fmt) if pattern else None))
    return cmds


def _outputs_exist(name, task, pid, sid, dirs):
    """True where a step has a status file or any output its deletion names."""
    fstatus = _step_path(name, pid, sid, dirs, task)
    if fstatus and os.path.exists(fstatus):
        return True
    fmt = dict(p=f'{pid}/{sid}', s=f'{pid}_{sid}', task=task)
    if name == 'analyze_func':
        root = dirs.get('analysis_rootdir')
        return bool(root and glob.glob(os.path.join(
            root, pid, sid, 'func', f'{pid}_{sid}_task-{task}_*')))
    for dir_key, subdir, pattern in STALE_STEPS[name]['delete']:
        if not dirs.get(dir_key):
            continue
        path = os.path.join(dirs[dir_key], subdir.format(**fmt))
        if pattern is None and os.path.isdir(path):
            return True
        if pattern and glob.glob(os.path.join(path, pattern.format(**fmt))):
            return True
    return False


def _find_stale_steps(pid, sid, dirs, config, opts, tasks):
    """The stale steps of a session's job: a list of dicts with step, task,
    reason, input, input_time, startdate and the csh deletion lines (cmds).

    In run order, a step the job runs is stale where (i) an input step with
    errcode 0 ran after it started, (ii) for analyze_func, an events copy in
    proc func/ was made after it started, (iii) for a proc step, its scans
    from the current scans.tsv differ from those it processed, or (iv) an
    input step the job also runs is stale (the cascade, which also takes a
    step with outputs but no status file). Steps the job does not run are
    neither compared nor deleted."""
    run = _run_steps(pid, sid, dirs, config, opts, tasks)
    run_names = set(n for n, _ in run)
    stale_names = set()
    found = []
    for name, task in run:
        fstatus = _step_path(name, pid, sid, dirs, task)
        if fstatus is None:
            continue
        inputs = [i for i in STALE_STEPS[name]['inputs']]
        hit = None
        upstream = [i for i in inputs if i in stale_names and i in run_names]
        if upstream:
            if _outputs_exist(name, task, pid, sid, dirs):
                hit = ('after stale ' + upstream[0], upstream[0], '', '')
        elif os.path.isfile(fstatus):
            paths = [p for p in (_step_path(i, pid, sid, dirs) for i in inputs) if p]
            newer = cmig_write_status.inputs_newer(fstatus, paths, ok_only=True)
            if newer:
                fpath, ran, started = newer
                inp = next(i for i in inputs if _step_path(i, pid, sid, dirs) == fpath)
                hit = ('input ran after', inp, ran, started)
            elif name == 'analyze_func':
                newer = _events_copy_newer(pid, sid, dirs, task, fstatus)
                if newer:
                    hit = ('events copy made after', os.path.basename(newer[0]),
                           newer[1], newer[2])
            elif name in ('proc_anat', 'proc_func', 'proc_dwi'):
                why = _selection_changed(name, pid, sid, dirs)
                if why:
                    hit = (why, 'scans.tsv', '', cmig_write_status.status_time(
                        fstatus, 'startdate'))
        if hit:
            stale_names.add(name)
            reason, inp, ran, started = hit
            found.append({'step': name, 'task': task, 'reason': reason,
                          'input': inp, 'input_time': ran, 'startdate': started,
                          'cmds': _step_deletions(name, task, pid, sid, dirs,
                                                  config, opts)})
    return found


def _failed_step_deletions(pid, sid, dirs, config, opts, tasks):
    """errflag 1's deletion lines from the step table: each step the job runs
    whose status has an errcode other than 0, and each step the job runs
    downstream of one (with outputs), so one job reruns them all; None where
    none of the job's steps failed."""
    run = _run_steps(pid, sid, dirs, config, opts, tasks)
    run_names = set(n for n, _ in run)
    failed = False
    gone = set()
    cmds = []
    for name, task in run:
        fstatus = _step_path(name, pid, sid, dirs, task)
        if fstatus is None:
            continue
        status = cmig_write_status.read_status(fstatus)
        hit = status is not None and (status.get('errcode') or '').strip() != '0'
        failed = failed or hit
        if not hit and any(i in gone and i in run_names
                           for i in STALE_STEPS[name]['inputs']):
            hit = _outputs_exist(name, task, pid, sid, dirs)
        if hit:
            gone.add(name)
            for cmd in _step_deletions(name, task, pid, sid, dirs, config, opts):
                if cmd not in cmds:
                    cmds.append(cmd)
    return cmds if failed else None


def _stale_detail(found):
    """A one-line account of the stale steps, and the deletion lines without
    repeats."""
    names = []
    for f in found:
        label = f['step'] + (f'[{f["task"]}]' if f['task'] else '')
        if f['reason'] == 'input ran after':
            label += ' ({} ran {}, after {})'.format(
                f['input'], _fmt_time(f['input_time']), _fmt_time(f['startdate']))
        elif f['reason'] == 'events copy made after':
            label += ' ({} made {}, after {})'.format(
                f['input'], _fmt_time(f['input_time']), _fmt_time(f['startdate']))
        else:
            label += f' ({f["reason"]})'
        names.append(label)
    cmds = []
    for f in found:
        for c in f['cmds']:
            if c not in cmds:
                cmds.append(c)
    return 'stale: ' + '; '.join(names), cmds


def _fmt_time(t):
    if isinstance(t, datetime.datetime):
        return t.strftime('%Y-%m-%dT%H:%M:%S')
    return str(t or '')


def write_stale_rows(jobdir, pid, sid, found):
    """Append a session's stale steps to jobdir/stale.tsv."""
    fname = os.path.join(jobdir, 'stale.tsv')
    new = not os.path.isfile(fname)
    with open(fname, 'a') as fid:
        if new:
            fid.write('\t'.join(STALE_COLUMNS) + '\n')
        for f in found:
            row = [pid, sid, f['step'], f['task'], f['reason'], f['input'],
                   _fmt_time(f['input_time']), _fmt_time(f['startdate']),
                   ' '.join(f['cmds'])]
            fid.write('\t'.join(str(v).replace('\t', ' ') for v in row) + '\n')


def _check_session_status(pid, sid, dirs, config, opts=None, delay_flags=None):
    """Return (status, errmsg, failed_entry), one of four answers.

    stale   — an input ran after the output that read it
    failed  — every expected status file is there and one carries errcode != 0
    done    — every expected status file is there with errcode == 0
    pending — a status file is missing, so there is work left

    stale is tested first, so a session that is both stale and failed gets a
    job rather than sitting in the failure report: the failure is usually an
    input that was not there yet, which is what has since changed.

    failed_entry is (dir_key, subpath_template, delete_subdir) of the failing
    status file, else None; for a stale behavioral func/ status in the
    pipeline, the same, with its func/ directory to delete; for stale steps
    (stale_steps, with opts stale_flag not 0), a dict with found (the
    _find_stale_steps list) and cmds (the deletion lines, with the behavioral
    one where that is stale too). delete_subdir is the session subdirectory errflag 1
    deletes: None for the module's whole session directory; for the behavioral
    step's status files (events_status), only their own directory
    ({pid}/{sid}/func or beh), so a func/ failure reruns the func analysis and
    a beh/ one leaves func/ alone.
    """
    if _check_stale(pid, sid, dirs, config):
        return 'stale', None, None
    tasks = _session_tasks(pid, sid, dirs, opts)
    stale = None
    if config.get('stale_steps') and (opts or {}).get('stale_flag', 1) != 0:
        found = _find_stale_steps(pid, sid, dirs, config, opts, tasks)
        if found:
            why, cmds = _stale_detail(found)
            stale = {'found': found, 'cmds': cmds, 'why': why}
    if config.get('events_newer'):
        # an output computed from the rawdata events, older than them
        dir_key, subpath_tmpl = config['events_newer']
        if dirs.get(dir_key):
            fstatus = os.path.join(dirs[dir_key], subpath_tmpl.format(pid=pid, sid=sid))
            started = cmig_write_status.status_time(fstatus, 'startdate')
            newest = cmig_select_bids_func_runs.events_mtime(
                dirs['raw_rootdir'], pid, sid)
            if started is not None and newest is not None and newest > started:
                return 'stale', None, None

    groups = _expand_status_files(
        _status_files(config, tasks, opts, pid, sid, dirs), tasks)
    # a pipeline whose func steps are off has no behavioral step either; it
    #   runs, as the pipeline does, for events with or without bold data
    events_on = (not config.get('func_status_files') or
                 (opts or {}).get('func_flag', 1) != 0)
    if config.get('events_status') and events_on:
        # one more group: every rawdata directory holding a task's events must
        # have its status file, and all must have succeeded. A pipeline's
        # --task names bold tasks, so its behavioral tasks are all with events
        dir_key, subpath_tmpl = config['events_status']
        task = opts.get('task') if opts and not config['status_files'] else None
        etasks = task if isinstance(task, list) else ([task] if task else None)
        found = _events_tasks(pid, sid, dirs['raw_rootdir'], etasks)
        events_group = [(dir_key, subpath_tmpl.replace('{mod}', mod), tsk)
                        for tsk in found for mod in found[tsk]]
        # a behavioral status older than its events, or run with another
        #   REDCap flag, is stale; the behavioral step makes the same check
        for (gkey, gtmpl, tsk), mod in zip(
                events_group, [m for t2 in found for m in found[t2]]):
            if not dirs.get(gkey):
                continue
            fstatus = os.path.join(dirs[gkey], gtmpl.format(pid=pid, sid=sid, task=tsk))
            flag = (get_ignore_delay_flag(delay_flags, pid, sid, tsk)
                    if delay_flags is not None else None)
            why = cmig_select_bids_func_runs.behav_status_stale(
                fstatus, dirs['raw_rootdir'], pid, sid, tsk, mod, flag,
                int((opts or {}).get('exclude_noncompliant_flag') or 0))
            if why:
                # a func/ one has changed the run selection the func analysis
                #   read, whose outputs then stand unless removed: the
                #   pipeline (status_files beside it) deletes its func/ dir
                delete_subdir = ('{}/{}/{}'.format(pid, sid, mod)
                                 if mod == 'func' and config['status_files'] else None)
                if stale:
                    if delete_subdir:
                        cmd = 'rm -rf ' + os.path.join(dirs[gkey], delete_subdir)
                        if cmd not in stale['cmds']:
                            stale['cmds'].append(cmd)
                    return 'stale', stale['why'] + '; ' + why, stale
                return 'stale', why, (gkey, gtmpl, delete_subdir)
        if events_group:
            groups.append(events_group)
    else:
        events_group = None
    if stale:
        return 'stale', stale['why'], stale

    n_ok = 0
    for group in groups:
        n_expected = 0
        n_group_ok = 0
        n_present = 0
        failed_entry = None
        failed_errmsg = None
        for dir_key, subpath_tmpl, tsk in group:
            rootdir = dirs.get(dir_key)
            if not rootdir:
                continue
            n_expected += 1
            fpath = os.path.join(rootdir, subpath_tmpl.format(pid=pid, sid=sid, task=tsk))
            status = cmig_write_status.read_status(fpath)
            if status is None:
                continue
            n_present += 1
            errcode = int(status.get('errcode', -1))
            if errcode != 0:
                if failed_entry is None:
                    delete_subdir = (os.path.dirname(subpath_tmpl).format(pid=pid, sid=sid)
                                     if group is events_group else None)
                    failed_entry = (dir_key, subpath_tmpl, delete_subdir)
                    failed_errmsg = status.get('errmsg', '')
            else:
                n_group_ok += 1
        # every task in the group must have succeeded: a session whose sst
        # output is gone has work left, whatever its other tasks did. A group
        # with one absent is pending rather than failed, so the job is made
        # without the deletes errflag does for a failure
        if n_expected and n_group_ok == n_expected:
            n_ok += 1
        elif n_present == n_expected and failed_entry is not None:
            return 'failed', failed_errmsg, failed_entry
    if groups and n_ok == len(groups):
        return 'done', None, None
    return 'pending', None, None


# ---------------------------------------------------------------------------
# Job script writer
# ---------------------------------------------------------------------------

def _write_deletions(f, delete_dirs):
    if delete_dirs:
        # a directory, or a csh line of its own (stale steps)
        for d in delete_dirs:
            f.write(f'{d}\n' if ' ' in d else f'rm -rf {d}\n')
        f.write('\n')


def _bind_path(bind):
    """The host path of a bind that mounts it at the same path in the image
    (src or src:src[:opts]), else None."""
    parts = bind.split(':')
    if len(parts) == 1 or parts[1] == parts[0]:
        return os.path.normpath(parts[0])
    return None


def _image_binds(defaults, given):
    """The binds of an image job: the given ones, and each default not
    already under a given one."""
    paths = [p for p in (_bind_path(b) for b in given or []) if p]
    binds = [d for d in sorted(set(defaults))
             if not any(d == p or d.startswith(p.rstrip('/') + '/') for p in paths)]
    for b in given or []:
        if b not in binds:
            binds.append(b)
    return binds


def _write_image_job_script(fpath, script_name, pid, sid, dirs, config, opts,
                            delete_dirs=None):
    """A job running the script in the image (opts['image']) through run.py,
    in the BIDS app form; the deletions run on the host first. bids_dir is the
    script's scan root (rawdata, or the input root of a step that does not
    read rawdata); the roots are resolved (realpath), so a root given by a
    symbolic link is bound and passed at its real path."""
    keys = [k for k in config['dir_args'] + config.get('optional_dir_args', [])
            if dirs.get(k) and k != 'raw_rootdir']
    roots = dict((k, os.path.realpath(dirs[k])) for k in keys)
    bids_dir = os.path.realpath(dirs[config['scan_dir_arg']])
    outputs = [roots[k] for k in _output_roots(config) if k in roots]
    # output_dir holds no root of its own here: every root is given
    output_dir = os.path.commonpath(outputs) if outputs else bids_dir
    arg_lines = [f'  {bids_dir} {output_dir} participant',
                 f'  --module {script_name}',
                 f'  --participant_label {pid} --session_label {sid}']
    for key in keys:
        arg_lines.append(f'  --{IMAGE_ROOTS[key]} {roots[key]}')
    names = config.get('image_arg_names', {})
    for key in config['extra_args']:
        val = opts[key]
        if val is None or key in IMAGE_HOST_OPTS:
            continue
        arg_lines.append(f'  --{names.get(key, key)} {" ".join(str(v) for v in val) if isinstance(val, list) else val}')
    if opts.get('dataset_name'):
        arg_lines.append(f'  --dataset_name {shlex.quote(opts["dataset_name"])}')
    binds = _image_binds([bids_dir] + [roots[k] for k in keys if os.path.isdir(roots[k])],
                         opts.get('image_binds'))
    bind_args = ' '.join(f'-B {b}' for b in binds)
    if opts.get('fs_license'):
        bind_args += f' -B {opts["fs_license"]}:{IMAGE_FS_LICENSE}:ro'
        arg_lines.append(f'  --fs_license {IMAGE_FS_LICENSE}')
    with open(fpath, 'w') as f:
        f.write('#!/bin/csh\n\n')
        _write_deletions(f, delete_dirs)
        # --cleanenv keeps the host's environment out; the settings of
        #   reproducible_flag 1 are passed in, for the python steps too
        env = {}
        if opts.get('reproducible_flag') != 0:
            env.update(cmig_sing_utils.REPRODUCIBLE_ENV)
        # the image's sha256, for the status files (add_image_info)
        if opts.get('image_sha256'):
            env[cmig_sing_utils.IMAGE_SHA256_ENV] = opts['image_sha256']
        for key, val in env.items():
            f.write(f'setenv SINGULARITYENV_{key} {val}\n')
        if env:
            f.write('\n')
        f.write(f'singularity exec --cleanenv --no-home {bind_args} \\\n')
        f.write(f'  {opts["image"]} /run.py \\\n')
        for line in arg_lines[:-1]:
            f.write(f'{line} \\\n')
        f.write(f'{arg_lines[-1]}\n')
    os.chmod(fpath, 0o755)


def _write_job_script(fpath, script_name, pid, sid, dirs, config, opts, delete_dirs=None):
    if opts.get('image'):
        return _write_image_job_script(fpath, script_name, pid, sid, dirs, config,
                                       opts, delete_dirs=delete_dirs)
    if config.get('image_only'):
        raise ValueError(IMAGE_ONLY_MSG.format(script_name))
    # a MATLAB module runs through the generic wrapper, given its name
    if script_name in cmig_run_bids_module.list_names():
        script_ref = os.path.join('$pydir', f'cmig_run_bids_module.py {script_name}')
    else:
        script_ref = os.path.join('$pydir', f'{script_name}.py')

    output_key = config.get('output_rootdir_from')
    arg_lines = [
        f'  --participant_id {pid}',
        f'  --session_id {sid}',
    ]
    for key in config['dir_args'] + config.get('optional_dir_args', []):
        # an optional rootdir left unset is left off, so the script applies its
        # own default rather than being handed the string None
        if dirs.get(key) is None:
            continue
        if key == output_key:
            arg_lines.append(f'  --output_rootdir {dirs[key]}')
        else:
            arg_lines.append(f'  --{key} {dirs[key]}')
    for key in config['extra_args']:
        val = opts[key]
        if val is None:
            continue
        arg_lines.append(f'  --{key} {" ".join(str(v) for v in val) if isinstance(val, list) else val}')
    if config.get('mrt_flag', 1):
        for key in ('appver', 'appdir', 'pubsdir', 'mrtdir', 'mrtdir_target'):
            arg_lines.append(f'  --{key} {opts[key]}')

    with open(fpath, 'w') as f:
        f.write('#!/bin/csh\n\n')
        f.write(f'set pydir = {opts["pydir"]}\n')
        f.write('setenv PYTHONPATH "${pydir}:${CMIG_PIPELINE_DIR}/python"\n')
        if opts.get('matlab_addpath'):
            f.write(f'setenv CMIG_MATLAB_ADDPATH {opts["matlab_addpath"]}\n')
        f.write('\n')
        _write_deletions(f, delete_dirs)
        f.write(f'{opts.get("python_cmd") or DEFAULT_PYTHON_CMD} -u {script_ref} \\\n')
        for line in arg_lines[:-1]:
            f.write(f'{line} \\\n')
        f.write(f'{arg_lines[-1]}\n')

    os.chmod(fpath, 0o755)


def _output_roots(config):
    """Return the rootdir keys a script writes to, in dir_args order.

    Single-module scripts write one root, named by output_rootdir_from, and
    those of optional_output_dir_args. Pipeline scripts write every dir_arg
    except raw_rootdir, which is always an input.

    Note this is derived rather than read from status_files, which can
    under-report: a pipeline's func status files are listed apart from them.
    """
    output_key = config.get('output_rootdir_from')
    if output_key:
        return [output_key] + config.get('optional_output_dir_args', [])
    return [key for key in config.get('dir_args', []) if key != 'raw_rootdir']


def _input_roots(config):
    """Return the rootdir keys a script reads from, in dir_args order."""
    outputs = _output_roots(config)
    return [key for key in config.get('dir_args', []) if key not in outputs]


# ---------------------------------------------------------------------------
# Per-session job and root setup, used by the batch loop below and by
# cmig_run_bids_rounds.py
# ---------------------------------------------------------------------------

def check_dirs(script_name, dirs):
    """Check the root dirs a script needs: create its output roots, require its
    input roots to exist. Raises ValueError."""
    config = SCRIPT_CONFIG[script_name]
    output_keys = _output_roots(config)
    for key in config['dir_args']:
        if not dirs.get(key):
            raise ValueError(f'{script_name} requires --{key}')
        if key in output_keys:
            try:
                os.makedirs(dirs[key], exist_ok=True)
            except OSError as e:
                raise ValueError(f'cannot create directory: {dirs[key]} (--{key}): {e}')
        elif not os.path.isdir(dirs[key]):
            raise ValueError(f'directory not found: {dirs[key]} (--{key})')
    for key in config.get('optional_output_dir_args', []):
        if dirs.get(key):
            try:
                os.makedirs(dirs[key], exist_ok=True)
            except OSError as e:
                raise ValueError(f'cannot create directory: {dirs[key]} (--{key}): {e}')


def add_image_info(opts, log=print):
    """With opts['image'], set opts['image_sha256'] (the image file's
    sha256, written into each job) and opts['image_versions'] (its
    versions.json, for dataset_description.json), once per batch; a failure
    is a warning, leaving '' and {}."""
    image = opts.get('image')
    if not image or 'image_sha256' in opts:
        return
    log(f'computing the sha256 of {image}')
    try:
        opts['image_sha256'] = cmig_sing_utils.image_sha256(image)
    except OSError as e:
        opts['image_sha256'] = ''
        log(f'WARNING: sha256 of {image} not computed: {e}')
    opts['image_versions'], msg = cmig_sing_utils.image_versions(image)
    if msg:
        log('WARNING: ' + msg)


def write_root_files(script_name, dirs, dataset_name='',
                     pipeline_name=cmig_write_bids_dataset_description.DEFAULT_PIPELINE_NAME,
                     software_name=cmig_write_bids_dataset_description.DEFAULT_SOFTWARE_NAME,
                     software_version=cmig_write_bids_dataset_description.DEFAULT_SOFTWARE_VERSION,
                     FS_version=DEFAULT_FS_VERSION, relpath_flag=0, atlas_flag=1,
                     versions=None):
    """Write dataset_description.json into the script's output roots (when
    dataset_name is given) and the atlas into proc_rootdir, once per set of
    roots rather than per job. versions, the package versions for
    GeneratedBy (an image's, from add_image_info; None: those of
    $CMIG_VERSIONS_FILE). Returns the messages to report: a root whose
    directory name BIDS does not accept (check_root_name), a file not
    written, or no dataset_name."""
    dd = cmig_write_bids_dataset_description
    output_keys = [k for k in _output_roots(SCRIPT_CONFIG[script_name]) if dirs.get(k)]
    msgs = []
    if dataset_name:
        kinds = dict((row[0], row[2]) for row in dd.ROOTS)
        for key in output_keys:
            if key in kinds:
                msg = dd.check_root_name(dirs[key], kinds[key], pipeline_name)
                if msg:
                    msgs.append('WARNING: ' + msg)
        # realpaths, so SourceDatasets and DatasetLinks resolve in the image
        #   (bound by realpath, --no-home), not through host links
        roots = dict((k, os.path.realpath(v)) for k, v in dirs.items()
                     if v and (k in kinds or k == 'raw_rootdir'))
        errcode, errmsg = dd.write_dataset_descriptions(
            dataset_name, write_roots=output_keys,
            pipeline_name=pipeline_name, software_name=software_name,
            software_version=software_version, FS_version=FS_version,
            relpath_flag=relpath_flag, versions=versions, **roots)
        if errcode:
            msgs.append('WARNING: ' + str(errmsg))
    else:
        msgs.append('no --dataset_name supplied; skipping dataset_description.json')
    # the atlas is the geometric reference named by cmig_proc_bids_anat, so it is
    # needed only where that output lands
    if atlas_flag and 'proc_rootdir' in output_keys:
        errcode, errmsg = cmig_write_bids_atlas.write_bids_atlas(dirs['proc_rootdir'])
        if errcode:
            msgs.append('WARNING: atlas not written: ' + str(errmsg))
    return msgs


def make_opts(gradunwarp_flag=2, cleanup_flag=1, task=None, acq_types=None,
              stim_times_model=None, stim_times_model_mid=None,
              stim_times_model_sst=None, stim_times_model_nback=None,
              stim_times_flag=None, baseline_flag=None, skipTRs=None,
              despike_flag=None, exclude_noncompliant_flag=None, T2_flag=0,
              FS_version=DEFAULT_FS_VERSION, compiled_flag=1,
              reproducible_flag=None, pc_flag=0, func_flag=None,
              dwi_flag=None, dwi_n_cpus=None, analyze_anat_flag=None,
              std_orient=None, appver=DEFAULT_APPVER, appdir=DEFAULT_APPDIR,
              pubsdir=DEFAULT_PUBSDIR, mrtdir=DEFAULT_MRTDIR,
              mrtdir_target=DEFAULT_MRTDIR_TARGET, pydir=DEFAULT_PYDIR,
              python_cmd=DEFAULT_PYTHON_CMD, stale_flag=1, image=None,
              image_binds=None, fs_license=None, dataset_name=None,
              matlab_addpath=None):
    """The job options make_session_job writes into a job script, and
    stale_flag, which it reads itself. With image (a singularity image file
    of the CMIG BIDS app), the job runs the script in it through run.py;
    image_binds, host directories bound in besides the job's roots;
    fs_license, the host's FreeSurfer license file, bound in; dataset_name,
    passed to run.py; matlab_addpath, dirs set as CMIG_MATLAB_ADDPATH in host
    jobs (compiled_flag 0: run working copies). image and fs_license are made
    absolute and must exist. Raises ValueError."""
    if image:
        image = os.path.abspath(image)
        if not os.path.exists(image):
            raise ValueError(f'image not found: {image}')
    if fs_license:
        fs_license = os.path.abspath(fs_license)
        if not os.path.isfile(fs_license):
            raise ValueError(f'FreeSurfer license not found: {fs_license}')
    return dict(
        gradunwarp_flag=gradunwarp_flag,
        cleanup_flag=cleanup_flag, task=task, acq_types=acq_types,
        stim_times_model=stim_times_model,
        stim_times_model_mid=stim_times_model_mid,
        stim_times_model_sst=stim_times_model_sst,
        stim_times_model_nback=stim_times_model_nback,
        stim_times_flag=stim_times_flag, baseline_flag=baseline_flag,
        skipTRs=skipTRs,
        despike_flag=despike_flag,
        exclude_noncompliant_flag=exclude_noncompliant_flag,
        T2_flag=T2_flag, FS_version=FS_version, compiled_flag=compiled_flag,
        reproducible_flag=reproducible_flag,
        pc_flag=pc_flag, func_flag=func_flag, dwi_flag=dwi_flag,
        dwi_n_cpus=dwi_n_cpus,
        analyze_anat_flag=analyze_anat_flag, std_orient=std_orient,
        ignore_delay_flag=None, ignore_delay_tasks=None,
        appver=appver, appdir=appdir, pubsdir=pubsdir,
        mrtdir=mrtdir, mrtdir_target=mrtdir_target, pydir=pydir,
        python_cmd=python_cmd, stale_flag=stale_flag, image=image,
        image_binds=image_binds, fs_license=fs_license,
        dataset_name=dataset_name, matlab_addpath=matlab_addpath,
    )


def load_delay_flags(script_name, tracking_csv=None, delay_flags_dir=None):
    """The REDCap ignore_delay flags for a script that takes them, else None.
    Raises ValueError when the script takes them and none can be read."""
    config = SCRIPT_CONFIG[script_name]
    if not set(['ignore_delay_flag', 'ignore_delay_tasks']) & set(config['extra_args']):
        return None
    delay_flags = load_ignore_delay_flags(tracking_csv, delay_flags_dir)
    if delay_flags is None:
        if delay_flags_dir is not None:
            raise ValueError(f'no mr_y_qc__raw__tfmri__{{mid,sst,nback}}.tsv in '
                             f'{delay_flags_dir}')
        raise ValueError(f'tracking csv not found: {tracking_csv} '
                         '(--tracking_csv none to run with every flag 0)')
    return delay_flags


def _status_step(subpath):
    """The step name of a status file: its datatype directory and its stem
    without the {pid}_{sid}_ prefix and _status.tsv suffix, e.g. anat/desc-proc,
    func/task-mid_desc-analyze; fsurf/desc-fsurf for the recon status
    ({pid}_{sid}/...), session/<stem> for a session-level one."""
    parts = subpath.split('/')
    stem = parts[-1]
    if stem.endswith('_status.tsv'):
        stem = stem[:-len('_status.tsv')]
    # the stem's sub- and ses- entities, whatever the ids
    stem = '_'.join(e for e in stem.split('_')
                    if not e.startswith('sub-') and not e.startswith('ses-'))
    if len(parts) >= 4:
        datatype = parts[-2]
    elif stem == 'desc-fsurf':
        datatype = 'fsurf'
    else:
        datatype = 'session'
    return '{0}/{1}'.format(datatype, stem)


def status_rows(script_name, pid, sid, dirs, opts=None):
    """Every status file a session is expected to have for a script, read,
    writing nothing: a list of dicts with step (as _status_step), dir_key,
    task, path, and status (the row as a dict, or None where the file is
    absent or unreadable). Unlike check_session(), which stops at the first
    failing file, this gives them all; the error catalog reads it."""
    config = SCRIPT_CONFIG[script_name]
    tasks = _session_tasks(pid, sid, dirs, opts)
    entries = [e for group in _expand_status_files(
        _status_files(config, tasks, opts, pid, sid, dirs), tasks) for e in group]
    # the behavioral step's files, as in _check_session_status()
    events_on = (not config.get('func_status_files') or
                 (opts or {}).get('func_flag', 1) != 0)
    if config.get('events_status') and events_on and dirs.get('raw_rootdir'):
        dir_key, subpath_tmpl = config['events_status']
        task = opts.get('task') if opts and not config['status_files'] else None
        etasks = task if isinstance(task, list) else ([task] if task else None)
        found = _events_tasks(pid, sid, dirs['raw_rootdir'], etasks)
        entries += [(dir_key, subpath_tmpl.replace('{mod}', mod), tsk)
                    for tsk in found for mod in found[tsk]]
    rows = []
    for dir_key, subpath_tmpl, tsk in entries:
        rootdir = dirs.get(dir_key)
        if not rootdir:
            continue
        subpath = subpath_tmpl.format(pid=pid, sid=sid, task=tsk)
        fpath = os.path.join(rootdir, subpath)
        rows.append({'step': _status_step(subpath), 'dir_key': dir_key,
                     'task': tsk, 'path': fpath,
                     'status': cmig_write_status.read_status(fpath)})
    return rows


def check_session(script_name, pid, sid, dirs, opts, delay_flags=None):
    """One session's status from its status files, writing nothing: (status,
    errmsg), status one of done, failed, stale, pending."""
    status, errmsg, _ = _check_session_status(
        pid, sid, dirs, SCRIPT_CONFIG[script_name], opts, delay_flags=delay_flags)
    return status, errmsg


def make_session_job(script_name, pid, sid, dirs, opts, jobdir, job_id,
                     errflag=0, newflag=True, delay_flags=None, verbose=True):
    """Check one session's status and, where it needs a job, write its job
    script jobdir/<job_id>.csh. Returns (status, errmsg, fpath):
      status  missing  a required input directory, or (behavioral step) any
                       task events, is absent; no job
              failed   a module failed; a job only with errflag >= 1, with
                       errflag's deletes
              done     every module succeeded; a job only with newflag False
              stale    an input ran after an output; a job
              pending  work left; a job
      fpath   the job script written, or None
    """
    config = SCRIPT_CONFIG[script_name]
    scan_dir = dirs[config['scan_dir_arg']]
    task = opts.get('task')

    missing = [s for s in config.get('required_subdirs', [])
               if not os.path.isdir(os.path.join(scan_dir, pid, sid, s))]
    if missing:
        if verbose:
            print(f'skipping {pid} {sid}: missing input directory: {missing[0]}')
        return 'missing', f'missing input directory: {missing[0]}', None

    job_opts = opts
    if 'ignore_delay_tasks' in config['extra_args']:
        # the tasks whose REDCap flag is set, or nothing
        job_opts = dict(opts, ignore_delay_tasks=[
            t for t in DELAY_TASKS
            if get_ignore_delay_flag(delay_flags or {}, pid, sid, t)] or None)
    elif 'ignore_delay_flag' in config['extra_args'] and \
            not config.get('events_status'):
        # a single-task job: its task's REDCap flag
        one = task[0] if isinstance(task, list) else task
        job_opts = dict(opts, ignore_delay_flag=(
            get_ignore_delay_flag(delay_flags or {}, pid, sid, one)
            if one else None))
    elif config.get('events_status'):
        tasks = task if isinstance(task, list) else ([task] if task else None)
        found = _events_tasks(pid, sid, scan_dir, tasks)
        if not found:
            if verbose:
                print(f'skipping {pid} {sid}: no task events in func/ or beh/')
            return 'missing', 'no task events in func/ or beh/', None
        job_tasks = list(found)
        job_opts = dict(opts, task=job_tasks, ignore_delay_flag=[
            get_ignore_delay_flag(delay_flags or {}, pid, sid, t)
            for t in job_tasks])

    session_status, errmsg, failed_entry = _check_session_status(
        pid, sid, dirs, config, opts, delay_flags=delay_flags)

    if session_status == 'failed' and errflag == 0:
        if verbose:
            print(f'skipping {pid} {sid}: prior error: {errmsg}')
        return 'failed', errmsg, None

    if session_status == 'done' and newflag:
        return 'done', None, None

    delete_dirs = []
    if session_status == 'stale':
        if verbose:
            print(f'{pid} {sid}: stale, {errmsg or "an input ran after the output"}')
        if isinstance(failed_entry, dict):
            delete_dirs += failed_entry['cmds']
            write_stale_rows(jobdir, pid, sid, failed_entry['found'])
        elif failed_entry and failed_entry[2] and not config.get('no_delete'):
            delete_dirs.append(os.path.join(dirs[failed_entry[0]], failed_entry[2]))
    if session_status == 'failed' and errflag >= 1 and not config.get('no_delete'):
        if errflag == 1:
            dir_key, subpath_tmpl, delete_subdir = failed_entry
            # a failed step of the step table: its outputs and those of the
            #   steps after it that the job runs; else (a behavioral step, a
            #   script without the step table) the failing module's directory
            cmds = None
            if not delete_subdir and config.get('stale_steps'):
                cmds = _failed_step_deletions(
                    pid, sid, dirs, config, opts,
                    _session_tasks(pid, sid, dirs, opts))
            if cmds:
                delete_dirs += cmds
            else:
                subdir = delete_subdir or _module_subdir(pid, sid, subpath_tmpl)
                delete_dirs.append(os.path.join(dirs[dir_key], subdir))
        else:
            seen = set()
            for dir_key, subpath_tmpl in _status_files(
                    config, _session_tasks(pid, sid, dirs, opts), opts,
                    pid, sid, dirs):
                if dir_key in seen or dir_key not in dirs:
                    continue
                seen.add(dir_key)
                subdir = _session_subdir(pid, sid, subpath_tmpl)
                delete_dirs.append(os.path.join(dirs[dir_key], subdir))

    fpath = os.path.join(jobdir, f'{job_id}.csh')
    _write_job_script(fpath, script_name, pid, sid, dirs, config, job_opts,
                      delete_dirs=delete_dirs or None)
    if verbose:
        print(f'created job: {job_id}')
    return session_status, errmsg, fpath


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def cmig_mkjobs_bids_pipeline(
    script_name,
    batchdir=DEFAULT_BATCHDIR,
    raw_rootdir=None,
    proc_rootdir=None,
    fsurf_rootdir=None,
    analysis_rootdir=None,
    rawqc_rootdir=None,
    postqc_rootdir=None,
    pc_rootdir=None,
    manualqc_rootdir=None,
    proc_anat_rootdir=None,
    gradunwarp_flag=2,
    cleanup_flag=1,
    task=None,
    acq_types=None,
    stim_times_model=None,
    stim_times_model_mid=None,
    stim_times_model_sst=None,
    stim_times_model_nback=None,
    stim_times_flag=None,
    baseline_flag=None,
    skipTRs=None,
    despike_flag=None,
    exclude_noncompliant_flag=None,
    T2_flag=0,
    FS_version=DEFAULT_FS_VERSION,
    compiled_flag=1,
    reproducible_flag=None,
    pc_flag=0,
    func_flag=None,
    dwi_flag=None,
    dwi_n_cpus=None,
    analyze_anat_flag=None,
    std_orient=None,
    appver=DEFAULT_APPVER,
    appdir=DEFAULT_APPDIR,
    pubsdir=DEFAULT_PUBSDIR,
    mrtdir=DEFAULT_MRTDIR,
    mrtdir_target=DEFAULT_MRTDIR_TARGET,
    pydir=DEFAULT_PYDIR,
    python_cmd=DEFAULT_PYTHON_CMD,
    matlab_addpath=None,
    newflag=True,
    errflag=0,
    max_jobs=None,
    seslist=None,
    sublist=None,
    study='',
    dataset_name='',
    pipeline_name=cmig_write_bids_dataset_description.DEFAULT_PIPELINE_NAME,
    software_name=cmig_write_bids_dataset_description.DEFAULT_SOFTWARE_NAME,
    software_version=cmig_write_bids_dataset_description.DEFAULT_SOFTWARE_VERSION,
    relpath_flag=0,
    atlas_flag=1,
    tracking_csv=None,
    delay_flags_dir=None,
    stale_flag=1,
    image=None,
    image_binds=None,
    fs_license=None,
):
    """Generate csh batch job scripts for a CMIG BIDS pipeline or module script.

    Parameters
    ----------
    script_name : str
        name of the Python script (without .py), e.g. 'cmig_pipeline_bids', or
        of a MATLAB module run through cmig_run_bids_module.py
    batchdir : str, optional
        parent directory; a job subdirectory is created inside it
        default: ~/batchdirs
    raw_rootdir : str, optional
        root directory of raw BIDS data
    proc_rootdir : str, optional
        root directory for processed output
    fsurf_rootdir : str, optional
        root directory for FreeSurfer output
    analysis_rootdir : str, optional
        root directory for analysis output
    rawqc_rootdir : str, optional
        root directory for raw QC output
        default: none
    postqc_rootdir : str, optional
        root directory for post-processing QC output
        default: none
    pc_rootdir : str, optional
        root directory for protocol compliance output
        default: none
    manualqc_rootdir : str, optional
        root directory of manual raw QC output
        default: none (the script's own default)
    proc_anat_rootdir : str, optional
        root directory of the processed T1w, for cmig_proc_bids_dwi
        default: none (the script's own default, proc_rootdir)
    cleanup_flag : int, optional
        passed to scripts that accept it
        default: 1
    task : str or list, optional
        task label(s) for functional scripts (str for single-task scripts, list for cmig_pipeline_bids)
        default: None
    acq_types : list of str, optional
        acquisition types to process for cmig_postqc_bids_fmap (e.g. ['func'] or ['func','dwi'])
        default: None (uses script default: both func and dwi)
    T2_flag : int, optional
        T2w handling for FreeSurfer: 0=ignore, 1=use if available, 2=require
        default: 0
    FS_version : int, optional
        FreeSurfer version number
        default: DEFAULT_FS_VERSION
    compiled_flag : int, optional
        run the MATLAB modules as MATLAB code (0) or compiled app (1)
        default: 1
    reproducible_flag : int, optional
        run the MATLAB modules with OpenBLAS and MKL settings that give the
        same results on every CPU type (1) or the inherited environment (0)
        default: None (not passed; the scripts' own default, 1)
    appver : str, optional
        container image version tag
        default: mmps_dev
    appdir : str, optional
        directory containing the compiled app
        default: /usr/pubsw/packages/MMPS/apps
    pubsdir : str, optional
        public software directory
        default: /usr/pubsw
    mrtdir : str, optional
        MATLAB runtime directory on host
        default: /usr/pubsw/packages/matlab/R2023a
    mrtdir_target : str, optional
        MATLAB runtime mount point inside container
        default: /opt/matlabruntime/R2023a
    pydir : str, optional
        directory containing the Python scripts
        default: directory of this file
    python_cmd : str, optional
        the python the job scripts run the script with, e.g. the
        cmig_python3 wrapper (csh/ of cmig_pipeline) of the python
        singularity image
        default: python3
    newflag : bool, optional
        if True, skip sessions already successfully completed
        default: True
    errflag : int, optional
        0 = skip sessions with prior errors (default)
        1 = create jobs for sessions with prior errors; delete the failed module's
            modality dir ({pid}/{sid}/<mod>), or session dir for a session-level status
        2 = create jobs for sessions with prior errors; delete all session output dirs
    max_jobs : int, optional
        maximum number of job scripts to create
        default: no limit
    seslist : list of str, optional
        list of session IDs to include (e.g. ['ses-00A', 'ses-02A'])
        default: all sessions
    sublist : list of str, optional
        list of participant IDs to include (e.g. ['sub-003RTV85'])
        default: all participants
    study : str, optional
        prefix for the job subdirectory name
        default: none
    dataset_name : str, optional
        dataset name for dataset_description.json, written at each output root
        and suffixed per root (e.g. 'ABCD Study' -> 'ABCD Study - cmig_proc');
        if not supplied, no dataset_description.json is written; with image,
        also passed to run.py
        default: none
    pipeline_name : str, optional
        name of the pipeline that produced the output, recorded as the first
        GeneratedBy.Name in dataset_description.json, which BIDS requires to
        be in each root's directory name
        default: cmig
    software_name : str, optional
        not written; accepted for compatibility
        default: cmig_pipeline
    software_version : str, optional
        version of the pipeline, recorded as GeneratedBy.Version
        default: $CMIG_PIPELINE_VERSION
    relpath_flag : int, optional
        write SourceDatasets paths relative to each root rather than absolute,
        assuming one master root dir
        default: 0
    atlas_flag : int, optional
        install the ABCD T1 atlas into proc_rootdir, for scripts that write there;
        also the way to opt out of the FreeSurfer dependency mri_convert brings
        default: 1
    tracking_csv : str, optional
        REDCap tracking table whose eprime_mismatch_ok_{task} flags become each
        job's ignore_delay_flag, for scripts that take one; 'none' sets every
        flag 0; this or delay_flags_dir required by those scripts
        default: None
    delay_flags_dir : str, optional
        directory of release-layout mr_y_qc__raw__tfmri__{task}.tsv tables whose
        __eprime__tdiff__ign_indicator columns give the same flags (the
        release's phenotype directory, or the output of
        cmig_write_eprime_mismatch_tables.py); used in place of tracking_csv
        default: None
    stale_flag : int, optional
        1 = a step whose input steps ran after it started (or, for a proc step,
        whose scans from scans.tsv changed) is stale: its job deletes its
        outputs and those of the steps after it that the job runs, listed in
        jobdir/stale.tsv; 0 = not checked
        default: 1
    image : str, optional
        singularity image file of the CMIG BIDS app: each job runs the script
        in it through run.py (bids_dir output_dir participant --module
        script_name, the session's labels and the roots as options; bids_dir
        is the script's scan root: rawdata, or the input root of a script
        that does not read it), rather than with the host's python and
        compiled app; the host settings (compiled_flag, pydir, appver,
        appdir, pubsdir, mrtdir, mrtdir_target, python_cmd) are not used in
        the job; the roots are resolved (realpath); its sha256, computed
        once, is passed to each job for the status files, and its
        versions.json is read for dataset_description.json. Required by
        cmig_proc_bids_dwi
        default: None
    image_binds : list of str, optional
        host directories bound into the image at the same paths, besides
        the defaults (bids_dir and the job's roots that exist, resolved); a
        default under one of them is dropped
        default: None
    fs_license : str, optional
        FreeSurfer license file, bound into the image and passed to run.py
        default: None (none; the image has no FS_LICENSE of its own)
    """
    if script_name not in SCRIPT_CONFIG:
        raise ValueError(f'unknown script: {script_name}. '
                         f'Known scripts: {sorted(SCRIPT_CONFIG)}')

    config = SCRIPT_CONFIG[script_name]

    if config.get('image_only') and not image:
        print('ERROR: ' + IMAGE_ONLY_MSG.format(script_name), file=sys.stderr)
        sys.exit(1)

    # cmig_analyze_bids_func takes one --task (required, single valued):
    #   a batch with none or several would make jobs that fail at argparse
    if script_name == 'cmig_analyze_bids_func':
        tlist = task if isinstance(task, list) else ([task] if task else [])
        if len(tlist) != 1:
            print('ERROR: cmig_analyze_bids_func takes exactly one --task '
                  '(given: {})'.format(' '.join(tlist) or 'none'), file=sys.stderr)
            sys.exit(1)

    dirs = {
        'raw_rootdir': raw_rootdir,
        'proc_rootdir': proc_rootdir,
        'fsurf_rootdir': fsurf_rootdir,
        'analysis_rootdir': analysis_rootdir,
        'rawqc_rootdir': rawqc_rootdir,
        'postqc_rootdir': postqc_rootdir,
        'pc_rootdir': pc_rootdir,
        'manualqc_rootdir': manualqc_rootdir,
        'proc_anat_rootdir': proc_anat_rootdir,
    }

    try:
        opts = make_opts(
            gradunwarp_flag=gradunwarp_flag, cleanup_flag=cleanup_flag, task=task,
            acq_types=acq_types, stim_times_model=stim_times_model,
            stim_times_model_mid=stim_times_model_mid,
            stim_times_model_sst=stim_times_model_sst,
            stim_times_model_nback=stim_times_model_nback,
            stim_times_flag=stim_times_flag, baseline_flag=baseline_flag,
            skipTRs=skipTRs, despike_flag=despike_flag,
            exclude_noncompliant_flag=exclude_noncompliant_flag, T2_flag=T2_flag,
            FS_version=FS_version, compiled_flag=compiled_flag,
            reproducible_flag=reproducible_flag, pc_flag=pc_flag,
            func_flag=func_flag, dwi_flag=dwi_flag, dwi_n_cpus=dwi_n_cpus,
            analyze_anat_flag=analyze_anat_flag,
            std_orient=std_orient, appver=appver,
            appdir=appdir, pubsdir=pubsdir, mrtdir=mrtdir,
            mrtdir_target=mrtdir_target, pydir=pydir, python_cmd=python_cmd,
            stale_flag=stale_flag, image=image, image_binds=image_binds,
            fs_license=fs_license, dataset_name=dataset_name,
            matlab_addpath=matlab_addpath)
    except ValueError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(1)

    check_dirs(script_name, dirs)
    add_image_info(opts)

    # describe the roots just created, before the session scan: a run that finds
    # no sessions still leaves the output roots behind, and they should not be
    # left undescribed
    for msg in write_root_files(
            script_name, dirs, dataset_name=dataset_name,
            pipeline_name=pipeline_name, software_name=software_name,
            software_version=software_version, FS_version=FS_version,
            relpath_flag=relpath_flag, atlas_flag=atlas_flag,
            versions=opts.get('image_versions')):
        print(msg)

    scan_dir = dirs[config['scan_dir_arg']]
    sessions = _find_sessions(scan_dir, sublist)
    if seslist is not None:
        sessions = [(pid, sid) for pid, sid in sessions if sid in seslist]
    if not sessions:
        print(f'No sub-*/ses-* pairs found in {scan_dir}')
        return

    if study:
        jobdir_name = f'{study}_{script_name}'
    else:
        jobdir_name = script_name
    jobdir = os.path.join(batchdir, jobdir_name)
    os.makedirs(jobdir, exist_ok=True)

    scriptlist_path = os.path.join(jobdir, 'scriptlist.txt')
    width = max(3, len(str(len(sessions))))

    # this run's stale steps only
    stale_path = os.path.join(jobdir, 'stale.tsv')
    if os.path.isfile(stale_path):
        os.remove(stale_path)

    try:
        delay_flags = load_delay_flags(script_name, tracking_csv, delay_flags_dir)
    except ValueError as e:
        print(f'ERROR: {e}', file=sys.stderr)
        sys.exit(1)

    n_created = 0
    n_done = 0
    n_failed = 0
    n_skipped = 0
    n_missing = 0
    n_stale = 0

    with open(scriptlist_path, 'w') as scriptlist:
        for pid, sid in sessions:
            if max_jobs is not None and n_created >= max_jobs:
                break
            job_num = str(n_created + 1).zfill(width)
            job_id = f'job_{job_num}_{pid}_{sid}'
            status, _, fpath = make_session_job(
                script_name, pid, sid, dirs, opts, jobdir, job_id,
                errflag=errflag, newflag=newflag, delay_flags=delay_flags)
            if status == 'missing':
                n_missing += 1
            elif status == 'failed' and fpath is None:
                n_failed += 1
            elif status == 'done' and fpath is None:
                n_done += 1
                n_skipped += 1
            if status == 'stale':
                n_stale += 1
            if fpath is not None:
                scriptlist.write(f'{job_id}\n')
                n_created += 1

    print(f'\n{os.path.basename(__file__)}: {n_created} jobs created in {jobdir}')
    if n_skipped:
        print(f'  skipped {n_done} already-completed sessions')
    if n_failed:
        print(f'  skipped {n_failed} sessions with prior module errors')
    if n_missing:
        print(f'  skipped {n_missing} sessions with missing input directories')
    if n_stale:
        print(f'  {n_stale} of the jobs are for stale sessions')
    if os.path.isfile(stale_path):
        import csv
        with open(stale_path) as fid:
            steps = [r['step'] for r in csv.DictReader(fid, delimiter='\t')]
        print('  stale steps ({}): {}'.format(stale_path, ', '.join(
            f'{n} {steps.count(n)}' for n in STEP_ORDER if n in steps)))


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Generate csh batch job scripts for a CMIG BIDS pipeline or module script.'
    )
    parser.add_argument('script_name',
                        help='Script name (e.g. cmig_pipeline_bids)')
    parser.add_argument('--batchdir', default=DEFAULT_BATCHDIR,
                        help='Parent directory for job subdirectory (default: %(default)s)')
    parser.add_argument('--raw_rootdir',
                        help='Root directory of raw BIDS data')
    parser.add_argument('--proc_rootdir',
                        help='Root directory for processed output')
    parser.add_argument('--fsurf_rootdir',
                        help='Root directory for FreeSurfer output')
    parser.add_argument('--analysis_rootdir',
                        help='Root directory for analysis output')
    parser.add_argument('--rawqc_rootdir',
                        help='Root directory for raw QC output')
    parser.add_argument('--postqc_rootdir',
                        help='Root directory for post-processing QC output')
    parser.add_argument('--pc_rootdir',
                        help='Root directory for protocol compliance output')
    parser.add_argument('--manualqc_rootdir',
                        help='Root directory for manual raw QC output '
                             '(default: the QC output rootdir)')
    parser.add_argument('--proc_anat_rootdir',
                        help='cmig_proc_bids_dwi: root directory of the processed '
                             'T1w (default: the proc root)')
    parser.add_argument('--gradunwarp_flag', type=int, default=2, choices=[0, 1, 2],
                        help='Gradient unwarp: 0=skip, 1=apply if gradwarp type known, '
                             '2=fail on unknown gradwarp type (default: 2)')
    parser.add_argument('--cleanup_flag', type=int, default=1, choices=[0, 1],
                        help='Cleanup intermediate files (default: 1)')
    parser.add_argument('--task', nargs='+', default=None,
                        help='Task label(s) for functional scripts (default: none)')
    parser.add_argument('--acq_types', nargs='+', default=None,
                        help='Acquisition types for cmig_postqc_bids_fmap (default: func dwi)')
    parser.add_argument('--stim_times_model', default=None,
                        help='Stimulus HRF model for 3dDeconvolve, e.g. SPMG1 or SPMG2 (default: none)')
    parser.add_argument('--stim_times_model_mid', default=None,
                        help='stim_times_model override for the MID task (default: none)')
    parser.add_argument('--stim_times_model_sst', default=None,
                        help='stim_times_model override for the SST task (default: none)')
    parser.add_argument('--stim_times_model_nback', default=None,
                        help='stim_times_model override for the nback task (default: none)')
    parser.add_argument('--stim_times_flag', type=int, default=None, choices=[0, 1, 2, 3],
                        help='Stimulus time file format for 3dDeconvolve (default: none)')
    parser.add_argument('--baseline_flag', type=int, default=None, choices=[0, 1],
                        help='Extract model-based polort-0 baseline for post-normalization (default: none)')
    parser.add_argument('--skipTRs', type=int, default=None,
                        help='TRs to skip at the start of each run (default: none, set from manufacturer)')
    parser.add_argument('--despike_flag', type=int, default=None, choices=[0, 1],
                        help='Run 3dDespike on each run before extraction (default: none)')
    parser.add_argument('--exclude_noncompliant_flag', type=int, default=None, choices=[0, 1],
                        help='Run selection and pc valid count treat runs with compliant 0 in the '
                             'rawdata _scans.tsv as invalid (default: none)')
    parser.add_argument('--T2_flag', type=int, default=0, choices=[0, 1, 2],
                        help='T2w handling for FreeSurfer: 0=ignore, 1=use if available, 2=require (default: 0)')
    parser.add_argument('--FS_version', type=int, default=DEFAULT_FS_VERSION,
                        help='FreeSurfer version number '
                             '(default: %d)' % DEFAULT_FS_VERSION)
    parser.add_argument('--compiled_flag', type=int, default=1, choices=[0, 1],
                        help='0 = MATLAB code, 1 = compiled app (default: 1)')
    parser.add_argument('--reproducible_flag', type=int, default=None, choices=[0, 1],
                        help='1 = OpenBLAS and MKL settings that give the same '
                             'results on every CPU type, 0 = the inherited '
                             'environment (default: none, the scripts\' own '
                             'default 1)')
    parser.add_argument('--std_orient', default=None,
                        help='cmig_proc_bids_func output orientation, e.g. LAS or LPS '
                             '(default: none, the module\'s own, LAS)')
    parser.add_argument('--func_flag', type=int, default=None, choices=[0, 1],
                        help='cmig_pipeline_bids: run the func steps (default: '
                             'none, the pipeline\'s own default 1)')
    parser.add_argument('--dwi_flag', type=int, default=None, choices=[0, 1],
                        help='cmig_pipeline_bids: run the dwi steps (default: '
                             'none, the pipeline\'s own default 1)')
    parser.add_argument('--dwi_n_cpus', type=int, default=None,
                        help='cmig_pipeline_bids: threads for the dwi processing tools '
                             '(default: none, the pipeline\'s own default 8)')
    parser.add_argument('--analyze_anat_flag', type=int, default=None, choices=[0, 1],
                        help='cmig_pipeline_bids: run cmig_analyze_bids_anat (default: '
                             'none, the pipeline\'s own default 1)')
    parser.add_argument('--pc_flag', type=int, default=0, choices=[0, 1],
                        help='Run protocol compliance check after rawqc (default: 0)')
    flags = parser.add_mutually_exclusive_group(required=True)
    flags.add_argument('--tracking_csv',
                       help='REDCap tracking table of the eprime_mismatch_ok '
                            'flags (DAL_ABCD_QC_daic_tracking.csv), or none '
                            'to set every flag 0')
    flags.add_argument('--delay_flags_dir',
                       help='directory of release-layout '
                            'mr_y_qc__raw__tfmri__{mid,sst,nback}.tsv holding '
                            'the same flags as __eprime__tdiff__ign_indicator '
                            '(a release phenotype directory, or the output of '
                            'cmig_write_eprime_mismatch_tables.py)')
    parser.add_argument('--appver', default=DEFAULT_APPVER,
                        help='Container image version tag (default: %(default)s)')
    parser.add_argument('--appdir', default=DEFAULT_APPDIR,
                        help='Directory containing the compiled app (default: %(default)s)')
    parser.add_argument('--pubsdir', default=DEFAULT_PUBSDIR,
                        help='Public software directory (default: %(default)s)')
    parser.add_argument('--mrtdir', default=DEFAULT_MRTDIR,
                        help='MATLAB runtime directory on host (default: %(default)s)')
    parser.add_argument('--mrtdir_target', default=DEFAULT_MRTDIR_TARGET,
                        help='MATLAB runtime mount point inside container (default: %(default)s)')
    parser.add_argument('--pydir', default=DEFAULT_PYDIR,
                        help='Directory containing Python scripts (default: %(default)s)')
    parser.add_argument('--matlab_addpath', default=None,
                        help='Directories (separated by :) the host jobs put first '
                             'on the MATLAB path with --compiled_flag 0, e.g. '
                             'working copies (sets CMIG_MATLAB_ADDPATH)')
    parser.add_argument('--python_cmd', default=DEFAULT_PYTHON_CMD,
                        help='Python the job scripts run, e.g. the cmig_python3 '
                             'wrapper of the python singularity image (default: '
                             '%(default)s)')
    parser.add_argument('--image', default=None,
                        help='Singularity image file of the CMIG BIDS app: the jobs '
                             'run the script in it through run.py, in the BIDS app '
                             'form, not with the host python and compiled app '
                             '(default: none)')
    parser.add_argument('--image_binds', nargs='+', default=None,
                        help='With --image, host directories bound in at the same '
                             'paths besides bids_dir and the job\'s roots that '
                             'exist, which are always bound (default: none)')
    parser.add_argument('--fs_license', default=None,
                        help='With --image, the FreeSurfer license file, bound in '
                             'and passed to run.py (default: none)')
    parser.add_argument('--stale_flag', type=int, default=1, choices=[0, 1],
                        help='Make jobs for steps whose input steps ran after them, '
                             'or whose scans in scans.tsv changed, deleting their '
                             'outputs and those of the steps after them in the job '
                             '(listed in stale.tsv) (default: 1)')
    parser.add_argument('--newflag', type=int, default=1, choices=[0, 1],
                        help='Skip already-completed sessions (default: 1)')
    parser.add_argument('--errflag', type=int, default=0, choices=[0, 1, 2],
                        help='0=skip prior-error sessions, 1=create job and delete the failed '
                             'module\'s modality dir, 2=create job and delete all session output '
                             '(default: 0)')
    parser.add_argument('--max_jobs', type=int, default=None,
                        help='Maximum number of jobs to create (default: no limit)')
    parser.add_argument('--seslist', nargs='+', default=None,
                        help='Session IDs to include (default: all sessions)')
    parser.add_argument('--sublist', nargs='+', default=None,
                        help='Participant IDs to include (default: all participants)')
    parser.add_argument('--study', default='',
                        help='Prefix for the job subdirectory name (default: none)')
    parser.add_argument('--dataset_name', default='',
                        help='Dataset name for dataset_description.json, suffixed per root '
                             "(e.g. 'ABCD Study' -> 'ABCD Study - cmig_proc'); "
                             'if not supplied, no dataset_description.json is written; '
                             'with --image, also passed to run.py')
    parser.add_argument('--pipeline_name',
                        default=cmig_write_bids_dataset_description.DEFAULT_PIPELINE_NAME,
                        help='Name of the pipeline that produced the output, recorded in '
                             'dataset_description.json (default: %(default)s)')
    parser.add_argument('--software_name',
                        default=cmig_write_bids_dataset_description.DEFAULT_SOFTWARE_NAME,
                        help='Name of the versioned software the pipeline ships inside '
                             '(default: %(default)s)')
    parser.add_argument('--software_version',
                        default=cmig_write_bids_dataset_description.DEFAULT_SOFTWARE_VERSION,
                        help='Version of that software (default: $CMIG_PIPELINE_VERSION)')
    parser.add_argument('--relpath_flag', type=int, default=0, choices=[0, 1],
                        help='Write SourceDatasets paths relative to each root (default: 0)')
    parser.add_argument('--atlas_flag', type=int, default=1, choices=[0, 1],
                        help='Install the ABCD T1 atlas into proc_rootdir (default: 1)')
    return parser.parse_args()


def main():
    args = parse_args()
    cmig_mkjobs_bids_pipeline(
        script_name=args.script_name,
        batchdir=args.batchdir,
        raw_rootdir=args.raw_rootdir,
        proc_rootdir=args.proc_rootdir,
        fsurf_rootdir=args.fsurf_rootdir,
        analysis_rootdir=args.analysis_rootdir,
        rawqc_rootdir=args.rawqc_rootdir,
        postqc_rootdir=args.postqc_rootdir,
        pc_rootdir=args.pc_rootdir,
        manualqc_rootdir=args.manualqc_rootdir,
        proc_anat_rootdir=args.proc_anat_rootdir,
        gradunwarp_flag=args.gradunwarp_flag,
        cleanup_flag=args.cleanup_flag,
        task=args.task,
        acq_types=args.acq_types,
        stim_times_model=args.stim_times_model,
        stim_times_model_mid=args.stim_times_model_mid,
        stim_times_model_sst=args.stim_times_model_sst,
        stim_times_model_nback=args.stim_times_model_nback,
        stim_times_flag=args.stim_times_flag,
        baseline_flag=args.baseline_flag,
        skipTRs=args.skipTRs,
        despike_flag=args.despike_flag,
        exclude_noncompliant_flag=args.exclude_noncompliant_flag,
        T2_flag=args.T2_flag,
        FS_version=args.FS_version,
        compiled_flag=args.compiled_flag,
        reproducible_flag=args.reproducible_flag,
        pc_flag=args.pc_flag,
        func_flag=args.func_flag,
        dwi_flag=args.dwi_flag,
        dwi_n_cpus=args.dwi_n_cpus,
        analyze_anat_flag=args.analyze_anat_flag,
        std_orient=args.std_orient,
        tracking_csv=args.tracking_csv,
        delay_flags_dir=args.delay_flags_dir,
        stale_flag=args.stale_flag,
        image=args.image,
        image_binds=args.image_binds,
        fs_license=args.fs_license,
        appver=args.appver,
        appdir=args.appdir,
        pubsdir=args.pubsdir,
        mrtdir=args.mrtdir,
        mrtdir_target=args.mrtdir_target,
        pydir=args.pydir,
        python_cmd=args.python_cmd,
        matlab_addpath=args.matlab_addpath,
        newflag=bool(args.newflag),
        errflag=args.errflag,
        max_jobs=args.max_jobs,
        seslist=args.seslist,
        sublist=args.sublist,
        study=args.study,
        dataset_name=args.dataset_name,
        pipeline_name=args.pipeline_name,
        software_name=args.software_name,
        software_version=args.software_version,
        relpath_flag=args.relpath_flag,
        atlas_flag=args.atlas_flag,
    )


if __name__ == '__main__':
    main()
