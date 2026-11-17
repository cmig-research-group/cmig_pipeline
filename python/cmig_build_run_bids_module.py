#!/usr/bin/env python3
"""
Build the compiled app cmig_run_bids_module, its option help, and the
CMIG BIDS app image (Docker, then Singularity).

Steps, from the same MATLAB tree ({pipeline_dir}/matlab, with
--matlab_funcdir, if given, taking precedence, as in the compile):
  1. read the option help of every module in cmig_run_bids_module.MODULES
     from its MATLAB help header, and stop if the help and the table
     disagree (an option in one but not the other)
  2. with --compile_flag 1, compile cmig_run_bids_module with
     cmig_compile_matlab (only --toolbox_dir, --gradient_dir and
     {pipeline_dir}/matlab on the path, no MMPS), through a MATLAB script
     run with run_matlab_script, into
     {appdir}/compiled_{ver}/cmig_run_bids_module
  3. write the help there as cmig_run_bids_module_help.json, for
     cmig_run_bids_module.py --help, which runs where the MATLAB code is not
     installed
  4. with --docker_flag 1, assemble the build context in --context_dir (the
     Dockerfile, files/ and documentation/ of --docker_srcdir, the parts of
     --pipeline_dir the image needs, with --python_dir, if given, taking
     precedence for python/, and in local/ the pieces the Dockerfile does
     not download by tag: the compiled app, the four atlas folders of
     --atlas_dir, the maps of --gradient_dir, and the TORTOISE programs and
     settings of --tortoise_dir named in docker/files/keep_tortoise.txt),
     then docker build it as
     {image}:{appver}, with CMIG_PIPELINE_VERSION from git describe of
     --pipeline_dir and CMIG_IMAGE_TAG {image}:{appver}; --prune_flag 1 (default) prunes FSL, AFNI, the MATLAB
     Runtime and FreeSurfer to their keep-lists (docker/files/keep_*.txt).
     The Dockerfile's pins (CMIG_PIPELINE_TAG/APP_SHA256, ATLAS_TAG/
     ATLAS_SHA256, GDF_TAG/GDF_SHA256, TORTOISE_TAG/TORTOISE_ASSET/
     TORTOISE_SHA256) apply unless --release_tag, --atlas_tag, --gdf_tag or
     --tortoise_tag (each with its --*_sha256, and --tortoise_asset)
     override them; an empty tag ('' or none) takes that piece from local/
  5. with --sing_flag 1, save that Docker image as
     {appdir}/images/{image}-{appver}.tar and build
     {appdir}/images/{image}-{appver}.sif from the tar
  6. with --python_sing_flag 1, docker build the Dockerfile's base stage
     (the OS, system packages, Python 3.12 and its packages) as
     {python_image}:{appver}, save it as
     {appdir}/images/{python_image}-{appver}.tar and build
     {appdir}/images/{python_image}-{appver}.sif from the tar: the python
     image the hosts run the pipeline python in (csh/cmig_python3)
With --assets_dir DIR and --release_tag TAG, write the release assets
cmig_run_bids_module_TAG.tar.gz (the compiled app) and
cmig_pipeline_atlases_TAG.tar.gz (the atlas folders) into DIR and print
their sha256, for the Dockerfile's pins at the release step.

In a help header, an option is "%   'name': description", continued by more
deeply indented comment lines, with "{default = ...}" giving its default; a
{task}-style placeholder documents a family of options.

Usage:
  cmig_build_run_bids_module.py [--appdir DIR] [--appver VER]
    [--pipeline_dir DIR] [--toolbox_dir DIR] [--gradient_dir DIR]
    [--toolbox_data_dir DIR]
    [--matlab_funcdir DIR] [--compile_flag 0|1]
    [--docker_flag 0|1 --docker_srcdir DIR] [--prune_flag 0|1] [--atlas_dir DIR]
    [--python_dir DIR]
    [--sing_flag 0|1] [--python_sing_flag 0|1] [--python_image NAME]
    [--release_tag TAG --app_sha256 SHA]
    [--atlas_tag TAG --atlas_sha256 SHA] [--gdf_tag TAG --gdf_sha256 SHA]
    [--tortoise_tag TAG --tortoise_asset NAME --tortoise_sha256 SHA]
    [--tortoise_dir DIR]
    [--assets_dir DIR --release_tag TAG]

Created:  10/13/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/27/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/02/2026 by Don Hagler
Prev mod: 11/03/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Prev mod: 11/15/2026 by Don Hagler
Last mod: 11/16/2026 by Don Hagler
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

import cmig_run_bids_module
import cmig_sing_utils
import cmig_write_scripts

APPNAME = 'cmig_run_bids_module'
# the tree this script is in (its parent, holding setup/versions.txt), else
#   the clone
TREE = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
# the folder holding the repo clones: CMIG_REPOS_DIR, else the one holding
#   this tree
REPOS = os.environ.get('CMIG_REPOS_DIR') or os.path.dirname(TREE)
DEFAULT_PIPELINE_DIR = (TREE if os.path.isfile(os.path.join(TREE, 'setup', 'versions.txt'))
                        else os.path.join(REPOS, 'cmig_pipeline_dev'))
DEFAULT_TOOLBOX_DIR = os.path.join(REPOS, 'surface_morphometry_toolbox')
DEFAULT_GRADIENT_DIR = os.path.join(REPOS, 'gradient_distortion_fields')
# the toolbox data, not in the toolbox repo, in its folder layout; the MMPS
#   copy until the toolbox data release asset exists
DEFAULT_TOOLBOX_DATA_DIR = '/usr/pubsw/packages/MMPS/external/external.2025.03.20/matlab/Projects'
# the atlases the modules read (the cmig_pipeline_atlases install: the files
#   the code reads, and NIfTI versions of T1T2), for builds that take them
#   from local/ and for the release asset
DEFAULT_ATLAS_DIR = '/usr/pubsw/packages/MMPS/cmig_pipeline_atlases/cmig_pipeline_atlases.0.1.0'
ATLAS_SUBDIRS = ['T1_Atlas', 'T1T2', 'DTI_Atlas', 'fuzzy_clusters']
# the license of the atlas files, in docker/assets, at the top of the atlas
#   release asset and of the local-build atlases
ATLAS_LICENSE = 'LICENSE_atlases.txt'
# the TORTOISE programs and settings/ for builds that take them from local/
#   (--tortoise_tag ''): the dwi step's bundle, whose bin/ has no DRTAMAS/
#   folder, so it no longer matches keep_tortoise.txt (release 4.3.0's layout)
#   (CMIG_TORTOISE_DIR, else none: --tortoise_dir is then required)
DEFAULT_TORTOISE_DIR = os.environ.get('CMIG_TORTOISE_DIR')
DEFAULT_MATLAB_VERSION = 'R2023a'
DEFAULT_IMAGE = 'cmig_bids'
DEFAULT_PYTHON_IMAGE = 'cmig_python'
# the parts of the cmig_pipeline tree the image reads or runs (the MATLAB
#   code is in the compiled app): the python, the dwi step's bash, the
#   parms, the binaries the modules run (reg, applyUnwarpB0, epic, with
#   RUNPATH $ORIGIN/../lib) and their libraries and license texts
PIPELINE_SUBDIRS = ['python', 'bash', 'parms', 'bin', 'lib']
# the Dockerfile's pins of the downloaded pieces: (tag ARG, sha256 ARG)
PINS = {'app': ('CMIG_PIPELINE_TAG', 'APP_SHA256'),
        'atlases': ('ATLAS_TAG', 'ATLAS_SHA256'),
        'gradient_distortion_fields': ('GDF_TAG', 'GDF_SHA256'),
        'tortoise': ('TORTOISE_TAG', 'TORTOISE_SHA256')}
# the Dockerfile's other pin ARGs, by piece
PIN_EXTRAS = {'tortoise': 'TORTOISE_ASSET'}


def parse_args():
    parser = argparse.ArgumentParser(
        description='Compile cmig_run_bids_module and write its option help.')
    parser.add_argument('--appdir', default=cmig_run_bids_module.DEFAULT_APPDIR,
                        help='apps directory; output in compiled_{ver}/%s '
                             '(default: %%(default)s)' % APPNAME)
    parser.add_argument('--appver', default='mmps_djh',
                        help='app version, e.g. mmps_djh for compiled_djh '
                             '(default: %(default)s)')
    parser.add_argument('--pipeline_dir', default=DEFAULT_PIPELINE_DIR,
                        help='cmig_pipeline tree (default: %(default)s)')
    parser.add_argument('--toolbox_dir', default=DEFAULT_TOOLBOX_DIR,
                        help='surface_morphometry_toolbox tree (default: %(default)s)')
    parser.add_argument('--gradient_dir', default=DEFAULT_GRADIENT_DIR,
                        help='gradient_distortion_fields tree (default: %(default)s)')
    parser.add_argument('--toolbox_data_dir', default=DEFAULT_TOOLBOX_DATA_DIR,
                        help='toolbox data in the toolbox folder layout, bundled '
                             'into the compiled app (default: %(default)s)')
    parser.add_argument('--atlas_dir', default=DEFAULT_ATLAS_DIR,
                        help='directory with the atlas folders the image needs '
                             '(default: %(default)s)')
    parser.add_argument('--matlab_funcdir', default=None,
                        help='directory of MATLAB functions taking precedence '
                             'over the pipeline tree, e.g. working copies (default: none)')
    parser.add_argument('--matlab_version', default=DEFAULT_MATLAB_VERSION,
                        help='MATLAB version for run_matlab_script (default: %(default)s)')
    parser.add_argument('--compile_flag', type=int, default=None, choices=[0, 1],
                        help='compile (1), or only check and write the help (0) '
                             '(default: 1, or 0 for a docker build whose app is '
                             'downloaded by tag or a build of the python image '
                             'alone)')
    parser.add_argument('--docker_flag', type=int, default=0, choices=[0, 1],
                        help='assemble the build context and docker build the image '
                             '(default: %(default)s)')
    parser.add_argument('--prune_flag', type=int, default=1, choices=[0, 1],
                        help='prune FSL, AFNI, the MATLAB Runtime and FreeSurfer '
                             'in the image to their keep-lists (1), or keep them '
                             'whole (0) (default: %(default)s)')
    parser.add_argument('--sing_flag', type=int, default=0, choices=[0, 1],
                        help='build a Singularity image from the Docker image '
                             '(default: %(default)s)')
    parser.add_argument('--python_sing_flag', type=int, default=0, choices=[0, 1],
                        help='docker build the Dockerfile\'s base stage as '
                             '{python_image}:{appver} and a Singularity image '
                             'from it (default: %(default)s)')
    parser.add_argument('--python_image', default=DEFAULT_PYTHON_IMAGE,
                        help='name of the python image; tagged with appver '
                             '(default: %(default)s)')
    parser.add_argument('--docker_srcdir', default=None,
                        help='directory with the Dockerfile and files/ '
                             '(default: docker/ of --pipeline_dir)')
    parser.add_argument('--context_dir', default=None,
                        help='build context directory, updated in place '
                             '(default: {appdir}/docker_context_{ver})')
    parser.add_argument('--python_dir', default=None,
                        help='python files taking precedence over the pipeline '
                             "tree's python/, e.g. working copies (default: none)")
    for piece, opt in [('app', 'release'), ('atlases', 'atlas'),
                       ('gradient_distortion_fields', 'gdf'),
                       ('tortoise', 'tortoise')]:
        parser.add_argument(f'--{opt}_tag', default=None,
                            help=f'override the Dockerfile\'s tag for the {piece} '
                                 '(\'\' or none: take it from local/) (default: '
                                 'the Dockerfile\'s pin)')
    parser.add_argument('--app_sha256', default=None,
                        help='sha256 of the app asset of --release_tag')
    parser.add_argument('--atlas_sha256', default=None,
                        help='sha256 of the atlas asset of --atlas_tag')
    parser.add_argument('--gdf_sha256', default=None,
                        help='sha256 of the --gdf_tag source archive')
    parser.add_argument('--tortoise_asset', default=None,
                        help='name of the TORTOISE release asset of --tortoise_tag')
    parser.add_argument('--tortoise_sha256', default=None,
                        help='sha256 of the --tortoise_asset')
    parser.add_argument('--tortoise_dir', default=DEFAULT_TORTOISE_DIR,
                        help='TORTOISE package directory (bin/, settings/) for '
                             'a build with no TORTOISE tag (default: %(default)s)')
    parser.add_argument('--assets_dir', default=None,
                        help='with --release_tag, write the release assets here')
    parser.add_argument('--image', default=DEFAULT_IMAGE,
                        help='image name; tagged with appver (default: %(default)s)')
    return parser.parse_args()


def find_mfile(module, matlab_dirs):
    """Return the path of module.m, searching each directory tree in order."""
    for d in matlab_dirs:
        for root, _, files in os.walk(d):
            if module + '.m' in files:
                return os.path.join(root, module + '.m')
    return None


def parse_help(fname):
    """Return {option: {'description': str, 'default': str}} from a help header."""
    with open(fname, errors='replace') as f:
        text = f.read().replace('\r', '')
    header = text.split('%%%%')[0]
    opts = {}
    name = None
    indent = 0
    for line in header.splitlines():
        m = re.match(r"^%(\s+)'([A-Za-z_0-9{}]+)':\s*(.*)$", line)
        if m:
            name = m.group(2)
            indent = len(m.group(1))
            opts[name] = {'description': m.group(3).strip(), 'default': ''}
            continue
        if name is None:
            continue
        m = re.match(r'^%(\s+)(.*)$', line)
        if not m or len(m.group(1)) <= indent or not m.group(2).strip():
            name = None
            continue
        cont = m.group(2).strip()
        d = re.match(r'^\{\s*default\s*=\s*(.*?)\s*\}$', cont)
        if d:
            opts[name]['default'] = d.group(1)
        else:
            opts[name]['description'] += ' ' + cont
    return opts


def module_help(matlab_dirs):
    """Return the option help of every module and the table-help mismatches."""
    out = {}
    problems = []
    for module, spec in sorted(cmig_run_bids_module.MODULES.items()):
        fname = find_mfile(module, matlab_dirs)
        if fname is None:
            problems.append(f'{module}.m not found')
            continue
        doc = parse_help(fname)
        entries = {}
        for key in spec['options']:
            entry = doc.get(key)
            if entry is None:
                # a {task}-style placeholder documents a family of options
                for pat, e in doc.items():
                    if '{' in pat and re.fullmatch(
                            re.sub(r'\\\{\w+\\\}', r'\\w+', re.escape(pat)), key):
                        entry = e
                        break
            if entry is None:
                problems.append(f'{module}: option {key} is not in its help')
                continue
            entries[key] = entry
        for pat in doc:
            if '{' in pat or pat in spec['options'] or pat in spec.get('excluded', []):
                continue
            problems.append(f'{module}: documented option {pat} is not in the table')
        out[module] = {'options': entries}
    return out, problems


def compile_app(args, outdir):
    """Compile cmig_run_bids_module into outdir; return True if it was written."""
    parms = {'funcname': APPNAME, 'pipeline_dir': args.pipeline_dir,
             'toolbox_dir': args.toolbox_dir, 'gradient_dir': args.gradient_dir,
             'toolbox_data_dir': args.toolbox_data_dir, 'outdir': outdir}
    if args.matlab_funcdir:
        parms['matlab_funcdir'] = args.matlab_funcdir
    fd, fname_m = tempfile.mkstemp(prefix='cmig_build_run_bids_module_', suffix='.m')
    os.close(fd)
    t0 = time.time()
    try:
        # cmig_compile_matlab is not on the MMPS path MATLAB starts with
        dirs = [os.path.join(args.pipeline_dir, 'matlab', 'cmig_utils')]
        if args.matlab_funcdir:
            dirs.append(args.matlab_funcdir)
        header = '\n'.join("addpath('{0}');".format(d) for d in dirs)
        cmig_write_scripts.cmig_write_matlab_script(
            fname_m, 'cmig_compile_matlab', [], tags=list(parms), parms=parms,
            header=header)
        subprocess.run(['run_matlab_script', fname_m, args.matlab_version])
    finally:
        os.remove(fname_m)
    fname_run = os.path.join(outdir, f'run_{APPNAME}.sh')
    return os.path.isfile(fname_run) and os.path.getmtime(fname_run) >= t0


def sync(src, dst):
    """Copy the tree src to dst, following links, deleting what src lacks;
    python caches, which git ignores, are left out."""
    os.makedirs(dst, exist_ok=True)
    # u+w: copies of read-only installs (e.g. the atlases) can be removed
    run(['rsync', '-aL', '--delete', '--chmod=u+w', '--exclude', '__pycache__',
         '--exclude', '*.pyc', src.rstrip('/') + '/', dst.rstrip('/') + '/'])


def rmtree(path):
    """Remove the tree path, read-only entries included."""
    def onerror(func, p, _):
        os.chmod(os.path.dirname(p), 0o755)
        if os.path.exists(p) and not os.path.isdir(p):
            os.chmod(p, 0o644)
        func(p)
    shutil.rmtree(path, onerror=onerror)


def run(cmd):
    print('  ' + ' '.join(cmd))
    if subprocess.run(cmd).returncode:
        sys.exit('Error: command failed: ' + ' '.join(cmd))


def dockerfile_pins(fname):
    """{ARG name: default} of the Dockerfile's ARG lines with a default."""
    pins = {}
    with open(fname) as f:
        for line in f:
            m = re.match(r'ARG\s+(\w+)=(\S*)\s*$', line)
            if m:
                pins[m.group(1)] = m.group(2)
    return pins


def build_pins(args):
    """{piece: (tag, sha256)}: the Dockerfile's pins with the options'
    overrides; an empty tag means the piece comes from local/. Also
    {ARG: value} of the pieces' other pin ARGs (PIN_EXTRAS)."""
    pins = dockerfile_pins(os.path.join(args.docker_srcdir, 'Dockerfile'))
    overrides = {'app': (args.release_tag, args.app_sha256),
                 'atlases': (args.atlas_tag, args.atlas_sha256),
                 'gradient_distortion_fields': (args.gdf_tag, args.gdf_sha256),
                 'tortoise': (args.tortoise_tag, args.tortoise_sha256)}
    extra_overrides = {'tortoise': args.tortoise_asset}
    out = {}
    extras = {}
    for piece, (tag_arg, sha_arg) in PINS.items():
        tag, sha = pins.get(tag_arg, ''), pins.get(sha_arg, '')
        extra = pins.get(PIN_EXTRAS[piece], '') if piece in PIN_EXTRAS else None
        if overrides[piece][0] is not None:
            tag = '' if overrides[piece][0] in ('', 'none') else overrides[piece][0]
            sha = overrides[piece][1] or ''
            if extra is not None:
                extra = extra_overrides[piece] or ''
        out[piece] = (tag, sha)
        if extra is not None:
            extras[PIN_EXTRAS[piece]] = extra if tag else ''
    return out, extras


def read_keep_list(fname):
    """[(type, path)] of a keep-list (files/keep_<pkg>.txt): d <dir> (its
    files, not its subdirectories), f <file>."""
    entries = []
    with open(fname) as f:
        for line in f:
            line = line.split('#')[0].strip()
            if line:
                typ, path = line.split(None, 1)
                entries.append((typ, path))
    return entries


def copy_keep_list(srcdir, dstdir, fname_keep):
    """Copy the keep-list entries of srcdir into dstdir (links followed),
    stopping on an entry srcdir lacks."""
    for typ, path in read_keep_list(fname_keep):
        src = os.path.join(srcdir, path)
        dst = os.path.join(dstdir, path)
        if typ == 'd':
            if not os.path.isdir(src):
                sys.exit(f'Error: {src} not found ({fname_keep})')
            os.makedirs(dst, exist_ok=True)
            for name in sorted(os.listdir(src)):
                if os.path.isfile(os.path.join(src, name)):
                    shutil.copy2(os.path.join(src, name), os.path.join(dst, name))
        elif typ == 'f':
            if not os.path.isfile(src):
                sys.exit(f'Error: {src} not found ({fname_keep})')
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
        else:
            sys.exit(f'Error: bad line in {fname_keep}: {typ} {path}')


def assemble_context(args, appdir_app, context_dir, pins, pieces=True):
    """Assemble the Docker build context the Dockerfile expects, from empty
    (an earlier context is removed, so nothing stale reaches the build): in
    local/, only the pieces with no tag; with pieces False (a build of the
    base stage alone), only the Dockerfile, files/ and documentation/."""
    if os.path.isdir(context_dir) and os.listdir(context_dir):
        if not os.path.isfile(os.path.join(context_dir, 'Dockerfile')):
            sys.exit(f'Error: {context_dir} is not empty and holds no Dockerfile; '
                     'not removing it')
        rmtree(context_dir)
    os.makedirs(context_dir)
    for name in ['Dockerfile', 'files', 'documentation']:
        src = os.path.join(args.docker_srcdir, name)
        if os.path.isdir(src):
            sync(src, os.path.join(context_dir, name))
        else:
            shutil.copy2(src, os.path.join(context_dir, name))
    local = os.path.join(context_dir, 'local')
    for piece in PINS:
        os.makedirs(os.path.join(local, piece))
    if not pieces:
        return
    if not pins['app'][0]:
        sync(appdir_app, os.path.join(local, 'app', APPNAME))
    if not pins['atlases'][0]:
        for sub in ATLAS_SUBDIRS:
            sync(os.path.join(args.atlas_dir, sub), os.path.join(local, 'atlases', sub))
        shutil.copy2(os.path.join(args.docker_srcdir, 'assets', ATLAS_LICENSE),
                     os.path.join(local, 'atlases', ATLAS_LICENSE))
    if not pins['gradient_distortion_fields'][0]:
        # the maps only, not the code beside them
        run(['rsync', '-aL', '--delete', '--chmod=u+w', '--include', '*.mgh', '--exclude', '*',
             args.gradient_dir.rstrip('/') + '/',
             os.path.join(local, 'gradient_distortion_fields') + '/'])
    if not pins['tortoise'][0]:
        if not args.tortoise_dir:
            sys.exit('Error: no TORTOISE tag: --tortoise_dir (or CMIG_TORTOISE_DIR) is required')
        # the programs and settings the image keeps, not the whole package
        copy_keep_list(args.tortoise_dir, os.path.join(local, 'tortoise'),
                       os.path.join(args.docker_srcdir, 'files', 'keep_tortoise.txt'))
    for sub in PIPELINE_SUBDIRS:
        sync(os.path.join(args.pipeline_dir, sub),
             os.path.join(context_dir, 'cmig_pipeline', sub))
    if args.python_dir:
        # overlaid without --delete: it holds only the files that differ
        run(['rsync', '-aL', '--include', '*.py', '--exclude', '*',
             args.python_dir.rstrip('/') + '/',
             os.path.join(context_dir, 'cmig_pipeline', 'python') + '/'])


def sha256sum(fname):
    h = hashlib.sha256()
    with open(fname, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_assets(args, appdir_app):
    """The release assets of --release_tag in --assets_dir, with their
    sha256: the compiled app (cmig_run_bids_module/ at the top) and the
    atlas folders with their license (at the top)."""
    tag = args.release_tag
    os.makedirs(args.assets_dir, exist_ok=True)
    assets_src = os.path.join(args.docker_srcdir, 'assets')
    assets = [(f'{APPNAME}_{tag}.tar.gz', [], os.path.dirname(appdir_app), [APPNAME]),
              (f'cmig_pipeline_atlases_{tag}.tar.gz', ['-C', assets_src, ATLAS_LICENSE],
               args.atlas_dir, ATLAS_SUBDIRS)]
    for name, first, topdir, members in assets:
        fname = os.path.join(args.assets_dir, name)
        # --dereference: the atlas folders may be links
        run(['tar', '-czhf', fname] + first + ['-C', topdir] + members)
        print(f'{name}  sha256 {sha256sum(fname)}')


def main():
    args = parse_args()
    if args.docker_srcdir is None:
        args.docker_srcdir = os.path.join(args.pipeline_dir, 'docker')
    # a docker build whose app is downloaded by tag needs no compiled app here
    app_pinned = bool(args.docker_flag and build_pins(args)[0]['app'][0])
    # a build of the python image alone needs no compiled app either
    python_only = bool(args.python_sing_flag and not args.docker_flag)
    if args.compile_flag is None:
        args.compile_flag = 0 if app_pinned or python_only else 1
    matlab_dirs = ([args.matlab_funcdir] if args.matlab_funcdir else []) + [
        os.path.join(args.pipeline_dir, 'matlab')]
    outdir = cmig_sing_utils._compiled_subdir(args.appdir, APPNAME, args.appver)

    doc, problems = module_help(matlab_dirs)
    if problems:
        for p in problems:
            print(f'ERROR: {p}', file=sys.stderr)
        sys.exit(f'{len(problems)} mismatches between the table in '
                 'cmig_run_bids_module.py and the MATLAB help; nothing built')

    if args.compile_flag:
        print(f'compiling {APPNAME} from {" then ".join(matlab_dirs)} into {outdir}')
        if not compile_app(args, outdir):
            sys.exit(f'Error: compile did not write run_{APPNAME}.sh in {outdir}')
    if args.compile_flag or not (app_pinned or python_only):
        os.makedirs(outdir, exist_ok=True)
        fname_help = os.path.join(outdir, cmig_run_bids_module.HELP_NAME)
        with open(fname_help, 'w') as f:
            json.dump(doc, f, indent=1, sort_keys=True)
            f.write('\n')
        print(f'wrote {fname_help}: {len(doc)} modules')

    if args.assets_dir:
        if not args.release_tag:
            sys.exit('Error: --assets_dir needs --release_tag')
        write_assets(args, outdir)

    tag = f'{args.image}:{args.appver}'
    python_tag = f'{args.python_image}:{args.appver}'
    if args.docker_flag or args.python_sing_flag:
        pins, extras = build_pins(args)
        if args.docker_flag and not pins['app'][0] and \
                not os.path.isfile(os.path.join(outdir, f'run_{APPNAME}.sh')):
            sys.exit(f'Error: no compiled app in {outdir}')
        context_dir = args.context_dir or os.path.join(
            args.appdir, 'docker_context_' + re.sub(r'^mmps_', '', args.appver.lower()))
        print(f'assembling the build context in {context_dir}')
        assemble_context(args, outdir, context_dir, pins, pieces=bool(args.docker_flag))
        version = subprocess.run(
            ['git', '-C', args.pipeline_dir, 'describe', '--tags', '--always', '--dirty'],
            stdout=subprocess.PIPE, universal_newlines=True).stdout.strip() or 'unknown'
        pin_args = []
        for piece, (tag_arg, sha_arg) in PINS.items():
            pin_args += ['--build-arg', f'{tag_arg}={pins[piece][0]}',
                         '--build-arg', f'{sha_arg}={pins[piece][1]}']
        for arg, val in extras.items():
            pin_args += ['--build-arg', f'{arg}={val}']
        # the image tag, which the image's status files record as app_version
        build_args = ['--build-arg', f'CMIG_PIPELINE_VERSION={version}',
                      '--build-arg', f'CMIG_IMAGE_TAG={tag}',
                      '--build-arg', f'PRUNE={args.prune_flag}'] + pin_args
    if args.docker_flag:
        print(f'building {tag} (CMIG_PIPELINE_VERSION {version}, '
              f'PRUNE {args.prune_flag}; ' + ', '.join(
                  f'{p} {t or "local"}' for p, (t, _) in pins.items()) + ')')
        run(['docker', 'build'] + build_args + ['-t', tag, context_dir])
    if args.python_sing_flag:
        print(f'building {python_tag} (the base stage)')
        run(['docker', 'build', '--target', 'base'] + build_args +
            ['-t', python_tag, context_dir])
    images = []
    if args.sing_flag:
        images.append((tag, args.image))
    if args.python_sing_flag:
        images.append((python_tag, args.python_image))
    for image_tag, name in images:
        fstem = os.path.join(args.appdir, 'images', f'{name}-{args.appver}')
        os.makedirs(os.path.dirname(fstem), exist_ok=True)
        # from a saved tar: singularity's docker-daemon client is older than
        #   the API the Docker daemon accepts
        print(f'saving {image_tag} as {fstem}.tar')
        run(['docker', 'save', image_tag, '-o', fstem + '.tar'])
        print(f'building {fstem}.sif from {fstem}.tar')
        run(['singularity', 'build', '--force', fstem + '.sif',
             f'docker-archive://{fstem}.tar'])


if __name__ == '__main__':
    main()
