#!/usr/bin/env python3
"""
Shared utilities for CMIG singularity wrapper scripts.

Created:  05/01/2026 by Don Hagler
Prev mod: 05/06/2026 by Don Hagler
Prev mod: 06/17/2026 by Don Hagler
Prev mod: 07/24/2026 by Don Hagler
Prev mod: 09/16/2026 by Don Hagler
Prev mod: 10/06/2026 by Don Hagler
Prev mod: 10/07/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/03/2026 by Don Hagler
Last mod: 11/12/2026 by Don Hagler
"""

import json
import os
import re
import shutil
import subprocess
import sys


def run_singularity(sif, pubsdir, mrtdir, mrtdir_target, indirs, outdirs, container_args):
    """Build and run a singularity command. indirs are bound read-only; outdirs are read-write."""
    cmd = [
        'singularity', 'run', '-e',
        '-B', f'{pubsdir}:/usr/pubsw:ro',
        '-B', f'{mrtdir}:{mrtdir_target}:ro',
    ]
    for d in indirs:
        cmd += ['-B', f'{d}:{d}:ro']
    for d in outdirs:
        cmd += ['-B', f'{d}:{d}']
    cmd += [sif] + container_args
    subprocess.run(cmd)


def check_singularity(sif):
    """Verify singularity is in PATH and the SIF image exists. Exit on failure."""
    if shutil.which('singularity') is None:
        sys.exit('Error: singularity not found in PATH')
    if not os.path.isfile(sif):
        sys.exit(f'Error: container image not found: {sif}')


def _compiled_subdir(appdir, appname, appver=''):
    """Return the compiled app subdirectory path, mirroring cmig_write_runtime_script.m logic.

    The mmps_ prefix is stripped anchored and the remainder lowercased, so that
    appver 'MMPS_DJH' and 'mmps_djh' both resolve to compiled_djh; the compiled
    dirs are lowercase because cmig_build_docker/cmig_compile_mmps derive them
    with lower(mmps_version).
    """
    subdir = os.path.join(appdir, appname)
    if not os.path.isdir(subdir):
        if appver:
            stripped = re.sub(r'^mmps_', '', appver.lower())
            compdir = os.path.join(appdir, f'compiled_{stripped}')
        else:
            compdir = os.path.join(appdir, 'compiled')
        subdir = os.path.join(compdir, appname)
    return subdir


def check_compiled(appdir, appname, appver=''):
    """Verify the compiled app run script exists. Exit on failure."""
    subdir = _compiled_subdir(appdir, appname, appver)
    if not os.path.isdir(subdir):
        sys.exit(f'Error: compiled app directory not found for {appname} in {appdir}')
    cmdstr = os.path.join(subdir, f'run_{appname}.sh')
    if not os.path.isfile(cmdstr):
        sys.exit(f'Error: run script not found: {cmdstr}')


def run_compiled(appdir, appname, mrtdir, container_args, appver='', env=None):
    """Run a compiled MATLAB app via its run_{appname}.sh script, in env
    (default: the inherited environment)."""
    subdir = _compiled_subdir(appdir, appname, appver)
    cmdstr = os.path.join(subdir, f'run_{appname}.sh')
    cmd = [cmdstr, mrtdir] + container_args
    subprocess.run(cmd, env=env)


# OpenBLAS and MKL settings that give the same results on every CPU type:
#   OpenBLAS 0.3.3 (FSL 6.0.4's) does not recognize newer CPUs (e.g. AMD Zen 3)
#   and falls back to other kernels, changing topup's results; HASWELL kernels
#   and MKL_CBWR=AVX2 need a CPU with AVX2; MATLAB's own results depend on its
#   thread count, so cmig_run_bids_module.m fixes it to CMIG_MATLAB_NTHREADS
REPRODUCIBLE_ENV = {
    'OPENBLAS_CORETYPE': 'HASWELL',
    'OPENBLAS_NUM_THREADS': '1',
    'MKL_CBWR': 'AVX2',
    'CMIG_MATLAB_NTHREADS': '1',
}


def reproducible_env():
    """Return a copy of the environment with REPRODUCIBLE_ENV added; a
    variable already set is left as it is."""
    env = dict(os.environ)
    for key, val in REPRODUCIBLE_ENV.items():
        env.setdefault(key, val)
    return env


def cmig_getenv(name):
    """Return an environment variable the pipeline needs; stop with an error
    naming it if it is not set."""
    val = os.environ.get(name)
    if not val:
        raise RuntimeError('environment variable {0} not set'.format(name))
    return val


DEFAULT_PIPELINE_NAME = 'cmig'
DEFAULT_SOFTWARE_NAME = 'cmig_pipeline'
# None: $CMIG_PIPELINE_VERSION, read when used
DEFAULT_SOFTWARE_VERSION = None

# the image's package versions, written at build time (Dockerfile)
VERSIONS_FILE_ENV = 'CMIG_VERSIONS_FILE'
# the image's tag (Dockerfile ENV) and the .sif's sha256 (set in the job by
#   mkjobs and the rounds)
IMAGE_TAG_ENV = 'CMIG_IMAGE_TAG'
IMAGE_SHA256_ENV = 'CMIG_IMAGE_SHA256'

# the packages of versions.json: key, status column, GeneratedBy Name; the
#   cmig_pipeline keys are software_version and the cmig GeneratedBy entry
PACKAGES = [
    ('freesurfer', 'freesurfer_version', 'FreeSurfer'),
    ('fsl', 'fsl_version', 'FSL'),
    ('afni', 'afni_version', 'AFNI'),
    ('ants', 'ants_version', 'ANTs'),
    ('c3d', 'c3d_version', 'Convert3D'),
    ('mcr', 'mcr_version', 'MATLAB Runtime'),
    ('tortoise', 'tortoise_version', 'TORTOISE'),
    ('synb0_disco', 'synb0_version', 'Synb0-DISCO'),
    ('bids_validator', 'bids_validator_version', 'bids-validator'),
    ('deno', 'deno_version', 'Deno'),
    ('atlases', 'atlases_version', 'cmig_pipeline_atlases'),
    ('gradient_distortion_fields', 'gradient_distortion_fields_version',
     'gradient_distortion_fields'),
    ('python', 'python_version', 'Python'),
]
PACKAGE_COLUMNS = [p[1] for p in PACKAGES]

# the status files' empty value
NA = 'n/a'

# columns appended to the status row; the MATLAB writes the row, python amends it,
# since cmig_write_status has no version argument and the modules are not passed one
VERSION_COLUMNS = (['app_name', 'app_version', 'app_sha256', 'pipeline_name',
                    'software_name', 'software_version'] + PACKAGE_COLUMNS)


def read_versions(fname=None):
    """The package versions ({key: version}) of versions.json (default
    $CMIG_VERSIONS_FILE); {} where there is none (the host) or it cannot be
    read."""
    fname = fname or os.environ.get(VERSIONS_FILE_ENV)
    if not fname or not os.path.isfile(fname):
        return {}
    try:
        with open(fname) as f:
            versions = json.load(f)
    except (OSError, ValueError):
        return {}
    if not isinstance(versions, dict):
        return {}
    return dict((str(k), str(v)) for k, v in versions.items() if v not in (None, ''))


def image_sha256(image):
    """The sha256 of an image file (hex)."""
    import hashlib
    h = hashlib.sha256()
    with open(image, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def image_versions(image):
    """The versions.json of a singularity image, read with singularity exec;
    (versions, None), or ({}, a message) where it cannot be read."""
    try:
        proc = subprocess.run(
            ['singularity', 'exec', '--cleanenv', '--no-home', image,
             'cat', '/opt/cmig/versions.json'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=600)
    except (OSError, subprocess.SubprocessError) as e:
        return {}, 'versions.json not read from {0}: {1}'.format(image, e)
    if proc.returncode:
        return {}, 'versions.json not read from {0}: {1}'.format(
            image, ' '.join(proc.stderr.split()) or 'exit {0}'.format(proc.returncode))
    try:
        versions = json.loads(proc.stdout)
    except ValueError as e:
        return {}, 'versions.json of {0} not valid: {1}'.format(image, e)
    if not isinstance(versions, dict):
        return {}, 'versions.json of {0} not an object'.format(image)
    return dict((str(k), str(v)) for k, v in versions.items() if v not in (None, '')), None


def add_status_versions(status_file, app_name=None, app_version=None,
                        pipeline_name=DEFAULT_PIPELINE_NAME,
                        software_name=DEFAULT_SOFTWARE_NAME,
                        software_version=DEFAULT_SOFTWARE_VERSION):
    """Append version columns to a status file written by MATLAB.

    Values that do not apply are written as n/a rather than omitted, so the
    header is the same shape everywhere. app_version does not apply where no
    compiled app or .sif is invoked; software_version applies everywhere, since
    the code ships in cmig_pipeline however it is invoked; None means
    $CMIG_PIPELINE_VERSION.

    In the image ($CMIG_IMAGE_TAG set), app_version is the image tag and
    app_sha256 the .sif's sha256 ($CMIG_IMAGE_SHA256, n/a where not set);
    the package columns come from $CMIG_VERSIONS_FILE, n/a without it.

    Readers zip header to values and look up by name, so the added columns are
    ignored by anything not expecting them.
    """
    if not os.path.isfile(status_file):
        return
    if software_version is None:
        software_version = cmig_getenv('CMIG_PIPELINE_VERSION')
    try:
        with open(status_file) as f:
            lines = f.read().splitlines()
        if len(lines) < 2:
            return
        header = lines[0].split('\t')
        values = lines[1].split('\t')
        if any(col in header for col in VERSION_COLUMNS):
            return
        image_tag = os.environ.get(IMAGE_TAG_ENV)
        if image_tag:
            app_version = image_tag
        versions = read_versions()
        added = [app_name, str(app_version) if app_version else '',
                 os.environ.get(IMAGE_SHA256_ENV, ''),
                 pipeline_name, software_name, software_version]
        added += [versions.get(key, '') for key, _, _ in PACKAGES]
        added = [v or NA for v in added]
        with open(status_file, 'w') as f:
            f.write('\t'.join(header + VERSION_COLUMNS) + '\n')
            f.write('\t'.join(values + added) + '\n')
    except OSError:
        return


def check_status(status_file, appname, app_name=None, app_version=None):
    """Read status TSV written by the container, report result, and return errcode.

    Also records which code produced the output, since the status row carries only
    procdate. app_name defaults to appname; pass it explicitly where the thing that
    did the work is not the wrapper, as with recon-all.
    """
    if not os.path.isfile(status_file):
        # a missing or truncated status file is this module's failure, not the
        #   caller's: a pipeline running several tasks carries on with the rest
        errmsg = 'status file not found: ' + status_file
        print(appname + ': ERROR: ' + errmsg, file=sys.stderr)
        return 1, errmsg
    with open(status_file) as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        errmsg = 'status file is incomplete: ' + status_file
        print(appname + ': ERROR: ' + errmsg, file=sys.stderr)
        return 1, errmsg
    header = lines[0].split('\t')
    values = lines[1].split('\t')
    status = dict(zip(header, values))
    add_status_versions(status_file, app_name=app_name or appname,
                        app_version=app_version)
    errcode = int(status.get('errcode', -1))
    # n/a is the writers' empty value
    errmsg = status.get('errmsg')
    if not errmsg or errmsg == NA:
        errmsg = None
    if errcode == 0:
        print(appname + ': completed successfully')
        if errmsg:
            print(appname + ': warning: ' + errmsg)
    else:
        print(appname + ': errcode=' + str(errcode), file=sys.stderr)
        if errmsg:
            print(appname + ': errmsg=' + errmsg, file=sys.stderr)
    return errcode, errmsg
