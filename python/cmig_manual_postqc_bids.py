#!/usr/bin/env python3
"""
Manual post-processing QC review of BIDS derivatives: the T2w (t2), the fMRI
(fmri), the dMRI (dmri) and the FreeSurfer recon (fsurf).

--modality takes any of the four, reviewed in turn for each session (default:
all four, those whose input roots are given). For each, opens tkmedit (or
freeview) on the volume to review with the T1w as its aux volume, running
the tcl script ABCD used for that modality (MMIL_postQC_Exam), which plays a
scan through the volume and flicks to the T1w on loading; for dmri, the
AtlasTrack fiber segmentation is overlaid. With --images 1, the images
cmig_postqc_bids_* write are shown beside it, in a matplotlib window, or with
--image_viewer (eog) where matplotlib has no window toolkit. Takes a
pass/fail, a 0-3 score per MMIL_postQC_Exam artifact type and notes (a
score of 3 fails the review without asking for the pass/fail, and a
reviewer's pass with a score of 3 counts as a fail in the consensus), and
writes them, as cmig_manual_rawqc_bids does for raw QC, into per-session
files under output_rootdir/{pid}/{sid}:

  {pid}_{sid}_desc-<modality>_manualpostqc.tsv  one row per reviewer, in
                                             the mr_y_qc__post__man__
                                             <modality>__review names
  {pid}_{sid}_desc-session_manualpostqc.tsv  the consensus, in the
                                             mr_y_qc__post__man__ names,
                                             recomputed on each review as
                                             MMIL_Summarize_postQC does
  {pid}_{sid}_desc-manualpostqc_status.tsv   procdate, moved by each review

The reviews in files of the older names ({pid}_{sid}_desc-manualpostqc
<modality>_qc.tsv, and a single {pid}_{sid}_desc-manualpostqc_qc.tsv with one
row per modality and reviewer and short column names) are read with the
others; the next save of the session writes them to the current files and
then removes the old ones, with the old desc-manualpostqc_session.tsv.

Sessions come from the input roots or from --selected_file; a modality
column there limits each session to the modalities listed for it, and a
pipeline_rootdir column gives the cmig_pipeline root holding the session's
inputs (cmig_proc, cmig_postqc, cmig_analysis, freesurfer), used in place of
the root options where it has a value. Other columns (the selection's
reasons and priority) are not read or shown. With --random_order 1
(default) the sessions are shuffled and the seed printed at the start; with
0 they are in the file's order (or sorted, from the roots).

A review takes a lock on the session and modality
({pid}_{sid}_desc-manualpostqc_{modality}_review.lock) when it starts, so a
second reviewer skips it meanwhile; it expires after --lock_hours and
--override_lock 1 takes it regardless. Each save holds
{pid}_{sid}_desc-manualpostqc_save.lock for its read-replace-write. Lock files
are created with O_EXCL, not flock, since reviewers share NFS across hosts.

tkmedit runs only under FreeSurfer 5.3 (and, on some servers, with MATLAB's
opengl libraries on LD_LIBRARY_PATH), which the caller sets up.

Created:  10/08/2026 by Don Hagler
Prev mod: 10/08/2026 by Don Hagler
Prev mod: 10/09/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/13/2026 by Don Hagler
"""

import argparse
import copy
import csv
import datetime
import glob
import os
import random
import re
import socket
import subprocess
import sys
import termios
import time
import tty
from collections import OrderedDict

import cmig_sing_utils

APPNAME = 'cmig_manual_postqc_bids'
NA = 'n/a'
MODALITIES = ['t2', 'fmri', 'dmri', 'fsurf']

# the artifact types rated per modality, from MMIL_postQC_Exam, and the
#   released column each is averaged into
FTYPES = {
    't2':    ['regT1', 'motion', 'inhomogeneity', 'Artifact'],
    'fmri':  ['regT1', 'B0warp', 'imgQual', 'cutoff'],
    'dmri':  ['regT1', 'B0warp', 'imgQual', 'cutoff', 'fseg'],
    'fsurf': ['motion', 'PialOver', 'WMUnder', 'inhomogeneity', 'Artifact'],
}
FTYPE_LABEL = {
    'regT1': 'registration to T1w', 'motion': 'motion',
    'inhomogeneity': 'inhomogeneity', 'Artifact': 'artifact',
    'B0warp': 'residual B0 distortion', 'imgQual': 'image quality',
    'cutoff': 'cutoff', 'fseg': 'fiber segmentation',
    'PialOver': 'pial overestimation', 'WMUnder': 'white matter underestimation',
}
FTYPE_COLUMN = {
    'regT1': 'regist', 'motion': 'mot', 'inhomogeneity': 'inhomog',
    'Artifact': 'artifact', 'B0warp': 'distort', 'imgQual': 'quality',
    'cutoff': 'cutoff', 'fseg': 'segmnt', 'PialOver': 'pialover',
    'WMUnder': 'wmunder',
}
# an artifact score at or above this fails the review
SEVERE_SCORE = 3
ALL_FTYPES = ['regT1', 'motion', 'inhomogeneity', 'Artifact', 'B0warp',
              'imgQual', 'cutoff', 'fseg', 'PialOver', 'WMUnder']
PREFIX = 'mr_y_qc__post__man__'

# the tcl scripts ABCD used (abcd_manual_postqc.m), in $CMIG_PIPELINE_DIR/tcl
TCL_SCRIPTS = {
    't2': 'qc_sMRI_T2_quick.tcl',
    'fmri': 'qc_fMRI_IM_quick.tcl',
    'dmri': 'qc_dMRI_fseg_quick.tcl',
    'fsurf': 'qc_recon.tcl',
}

# the FreeSurfer recon files the fsurf review loads
FSURF_FILES = ['mri/nu.mgz', 'mri/brainmask.mgz', 'mri/aparc+aseg.mgz',
               'surf/lh.white', 'surf/rh.white', 'surf/lh.pial',
               'surf/rh.pial']

# the key columns of the per-modality per-reviewer files
REVIEW_KEYS = ['participant_id', 'session_id', 'reviewer_id']
STATUS_COLUMNS = ['participant_id', 'session_id', 'procdate', 'errcode',
                  'errmsg']

LOCK_WAIT = 60     # seconds a save waits for the save lock
LOCK_STALE = 300   # seconds after which a save lock is broken

#############################################################################
# sessions and inputs
#############################################################################

# the input roots under a cmig_pipeline root (pipeline_rootdir)
PIPELINE_SUBDIRS = {'proc_rootdir': 'cmig_proc',
                    'postqc_rootdir': 'cmig_postqc',
                    'analysis_rootdir': 'cmig_analysis',
                    'fsurf_rootdir': 'freesurfer'}


def read_selected_file(fname):
    """(sessions, modalities, pipeline_roots, errmsg): the (pid, sid) of
    --selected_file in the file's order, once each; if it has a modality
    column, a dict of the modalities listed for each session, else None; if
    it has a pipeline_rootdir column, a dict of each session's non-empty
    value, else None. Only participant_id, session_id, modality and
    pipeline_rootdir are read."""
    sessions = []
    selected = {}
    pipeline_roots = {}
    with open(fname, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        if 'participant_id' not in fields or 'session_id' not in fields:
            return [], None, None, \
                '{0} lacks participant_id or session_id'.format(fname)
        has_mod = 'modality' in fields
        has_root = 'pipeline_rootdir' in fields
        for row in reader:
            pid = (row.get('participant_id') or '').strip()
            sid = (row.get('session_id') or '').strip()
            if not (pid and sid):
                continue
            if has_mod:
                mod = (row.get('modality') or '').strip()
                if mod not in MODALITIES:
                    print('{0}: WARNING: {1}_{2}: modality "{3}" is not one '
                          'of {4}; row ignored'.format(
                              APPNAME, pid, sid, mod, ', '.join(MODALITIES)))
                    continue
                selected.setdefault((pid, sid), set()).add(mod)
            root = (row.get('pipeline_rootdir') or '').strip()
            if has_root and root and root != 'n/a':
                pipeline_roots[(pid, sid)] = root
            if (pid, sid) not in sessions:
                sessions.append((pid, sid))
    return (sessions, (selected if has_mod else None),
            (pipeline_roots if has_root else None), '')


def get_sessions(args):
    """(pid, sid) of the sessions with the modality's input, narrowed by
    --participant_id and --session_id."""
    pid_pat = args.participant_id or 'sub-*'
    sid_pat = args.session_id or 'ses-*'
    if args.modality == 'fsurf':
        dirs = glob.glob(os.path.join(args.fsurf_rootdir,
                                      '{0}_{1}'.format(pid_pat, sid_pat)))
        sessions = []
        for d in dirs:
            name = os.path.basename(d)
            idx = name.find('_ses-')
            if os.path.isdir(d) and idx > 0:
                sessions.append((name[:idx], name[idx + 1:]))
        return sorted(set(sessions))
    rootdir, subdir = {'t2': (args.proc_rootdir, 'anat'),
                       'fmri': (args.postqc_rootdir, 'func'),
                       'dmri': (args.analysis_rootdir, 'dwi')}[args.modality]
    dirs = glob.glob(os.path.join(rootdir, pid_pat, sid_pat, subdir))
    return sorted(set((os.path.basename(os.path.dirname(os.path.dirname(d))),
                       os.path.basename(os.path.dirname(d)))
                      for d in dirs if os.path.isdir(d)))


def first_match(pattern):
    fnames = sorted(glob.glob(pattern))
    return fnames[0] if fnames else None


def get_proc_anat(args, pid, sid, suffix):
    """The desc-preproc T1w or T2w, else one of the older name without it."""
    anatdir = os.path.join(args.proc_rootdir, pid, sid, 'anat')
    stem = '{0}_{1}_'.format(pid, sid)
    return (first_match(os.path.join(anatdir, stem + '*_desc-preproc_' + suffix + '.nii.gz')) or
            first_match(os.path.join(anatdir, stem + '*_' + suffix + '.nii.gz')))


def get_t1w(args, pid, sid):
    return get_proc_anat(args, pid, sid, 'T1w')


def get_inputs(args, pid, sid):
    """(inputs, errmsg): the volume to review, the T1w and, for dmri, the
    fiber segmentation, as paths, or None and why not."""
    if args.modality == 'fsurf':
        fsurf_dir = os.path.join(args.fsurf_rootdir, '{0}_{1}'.format(pid, sid))
        missing = [f for f in FSURF_FILES
                   if not os.path.isfile(os.path.join(fsurf_dir, f))]
        if missing:
            return None, 'missing recon files: ' + ', '.join(missing)
        return {'fsurf_dir': fsurf_dir}, ''
    t1w = get_t1w(args, pid, sid)
    if not t1w:
        return None, 'no T1w in {0}'.format(args.proc_rootdir)
    stem = '{0}_{1}'.format(pid, sid)
    if args.modality == 't2':
        vol = get_proc_anat(args, pid, sid, 'T2w')
        if not vol:
            return None, 'no T2w'
        return {'volume': vol, 't1w': t1w}, ''
    if args.modality == 'fmri':
        vol = first_match(os.path.join(args.postqc_rootdir, pid, sid, 'func',
                                       stem + '_space-T1w_desc-meanNorm_boldref.nii.gz'))
        if not vol:
            return None, 'no _space-T1w_desc-meanNorm_boldref'
        return {'volume': vol, 't1w': t1w}, ''
    # dmri: the RSI RND and AtlasTrack's fiber segmentation, both written on
    #   the T1w grid by cmig_analyze_bids_dwi
    dwidir = os.path.join(args.analysis_rootdir, pid, sid, 'dwi')
    vol = first_match(os.path.join(dwidir, stem + '_space-T1w_rsiRND.nii.gz'))
    if not vol:
        return None, 'no _space-T1w_rsiRND in {0}'.format(dwidir)
    out = {'volume': vol, 't1w': t1w}
    dseg = first_match(os.path.join(dwidir,
                                    stem + '_space-T1w_atlas-at_dseg.nii.gz'))
    if dseg:
        out['fseg'] = dseg
    else:
        print('{0}: WARNING: no AtlasTrack dseg in {1}; no fiber overlay'.format(
            APPNAME, dwidir))
    return out, ''


def get_images(args, pid, sid):
    """The jpgs cmig_postqc_bids_* wrote for the modality, if any."""
    if not args.postqc_rootdir:
        return []
    sesdir = os.path.join(args.postqc_rootdir, pid, sid)
    pats = {'t2': ['anat/*T2w*.jpg'],
            'fmri': ['func/*.jpg', 'fmap/*acq-func*.jpg'],
            'dmri': ['dwi/*.jpg', 'fmap/*acq-dwi*.jpg'],
            'fsurf': ['fsurf/*.jpg']}[args.modality]
    fnames = []
    for pat in pats:
        fnames += sorted(glob.glob(os.path.join(sesdir, pat)))
    return fnames

#############################################################################
# viewers
#############################################################################

def find_tcl(args):
    tcl_dir = args.tcl_dir or os.path.join(cmig_sing_utils.cmig_getenv('CMIG_PIPELINE_DIR'), 'tcl')
    fname = os.path.join(tcl_dir, TCL_SCRIPTS[args.modality])
    if not os.path.isfile(fname):
        print('{0}: WARNING: tcl script not found: {1}'.format(APPNAME, fname))
        return None
    return fname


def find_lut(modality):
    """The color lookup table of the segmentation overlay: for dmri,
    CMIG_ColorLUT.txt, whose 21000 + fiber number codes the dseg carries."""
    if modality == 'dmri':
        return os.path.join(cmig_sing_utils.cmig_getenv('CMIG_PIPELINE_DIR'), 'parms',
                            'CMIG_ColorLUT.txt')
    return os.path.join(os.environ.get('FREESURFER_HOME', ''),
                        'FreeSurferColorLUT.txt')


def launch_viewer(args, pid, sid, inputs):
    """The viewer process on the review volume, with the T1w as aux."""
    if args.viewer == 'none':
        return None
    devnull = subprocess.DEVNULL if args.suppress else None
    env = dict(os.environ)
    if args.modality == 'fsurf':
        subj = os.path.basename(inputs['fsurf_dir'])
        if args.viewer == 'tkmedit':
            cmd = ['tkmedit', subj, 'nu.mgz', 'lh.white', '-aux', 'brainmask.mgz',
                   '-aux-surface', 'rh.white', '-segmentation', 'aparc+aseg.mgz']
            lut = find_lut('fsurf')
            if os.path.isfile(lut):
                cmd.append(lut)
            cmd += ['-mm-main', '0', '255', '-mm-aux', '0', '255']
            tcl = find_tcl(args)
            if tcl:
                cmd += ['-tcl', tcl]
            env['SUBJECTS_DIR'] = os.path.dirname(inputs['fsurf_dir'])
        else:
            mridir = os.path.join(inputs['fsurf_dir'], 'mri')
            surfdir = os.path.join(inputs['fsurf_dir'], 'surf')
            cmd = ['freeview',
                   '-v', os.path.join(mridir, 'nu.mgz'),
                   '-v', os.path.join(mridir, 'brainmask.mgz') +
                   ':colormap=heat:opacity=0.3',
                   '-v', os.path.join(mridir, 'aparc+aseg.mgz') +
                   ':colormap=lut:opacity=0.3',
                   '-f', os.path.join(surfdir, 'lh.white'),
                   os.path.join(surfdir, 'rh.white'),
                   os.path.join(surfdir, 'lh.pial') + ':edgecolor=red',
                   os.path.join(surfdir, 'rh.pial') + ':edgecolor=red',
                   '-zoom', '1.5']
        return subprocess.Popen(cmd, env=env, stdout=devnull, stderr=devnull)
    if args.viewer == 'tkmedit':
        cmd = ['tkmedit', '-f', inputs['volume'], '-aux', inputs['t1w']]
        if inputs.get('fseg'):
            cmd += ['-segmentation', inputs['fseg'], find_lut('dmri')]
        tcl = find_tcl(args)
        if tcl:
            cmd += ['-tcl', tcl]
    else:
        cmd = ['freeview', '-v', inputs['t1w'], inputs['volume']]
        if inputs.get('fseg'):
            cmd += ['-v', inputs['fseg'] + ':colormap=lut:lut=' +
                    find_lut('dmri') + ':opacity=0.4']
    return subprocess.Popen(cmd, env=env, stdout=devnull, stderr=devnull)


# matplotlib backends that cannot show a window
NONINTERACTIVE = ('agg', 'pdf', 'ps', 'svg', 'cairo', 'template')


def interactive_backend():
    """Whether matplotlib can open a window here: it falls back to agg where
    no GUI toolkit loads (on cmig42, tkinter lacks libtk8.5)."""
    try:
        import matplotlib
        import matplotlib.pyplot  # noqa: F401, selects the backend
        return matplotlib.get_backend().lower() not in NONINTERACTIVE
    except Exception:
        return False


def launch_image_viewer(args, fnames):
    """The external image viewer process on the images, or None."""
    try:
        return subprocess.Popen([args.image_viewer] + fnames,
                                stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
    except OSError as e:
        print('  WARNING: cannot run {0}: {1}'.format(args.image_viewer, e))
        return None


class ImageViewer(object):
    """Show a list of images one at a time, with prev/next and zoom, as
    cmig_manual_rawqc_bids does."""

    def __init__(self, fnames, title=''):
        self.fnames = fnames
        self.index = 0
        self.fig = None
        self.title = title

    def open(self):
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Button
        self.plt = plt
        self.fig = plt.figure(figsize=(11, 8))
        self.ax = self.fig.add_axes([0.02, 0.10, 0.96, 0.84])
        ax_prev = self.fig.add_axes([0.35, 0.02, 0.12, 0.05])
        ax_next = self.fig.add_axes([0.53, 0.02, 0.12, 0.05])
        self.b_prev = Button(ax_prev, 'Previous')
        self.b_next = Button(ax_next, 'Next')
        self.b_prev.on_clicked(self.prev)
        self.b_next.on_clicked(self.next)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.show()
        plt.show(block=False)

    def show(self):
        import matplotlib.image as mpimg
        self.ax.clear()
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        fname = self.fnames[self.index]
        try:
            self.ax.imshow(mpimg.imread(fname))
        except Exception as e:
            self.ax.set_title('failed to read image: {0}'.format(e))
        else:
            self.ax.set_title('{0}   [{1}/{2}]  {3}'.format(
                os.path.basename(fname), self.index + 1, len(self.fnames),
                self.title), fontsize=9)
        self.fig.canvas.draw_idle()

    def prev(self, event=None):
        self.index = (self.index - 1) % len(self.fnames)
        self.show()

    def next(self, event=None):
        self.index = (self.index + 1) % len(self.fnames)
        self.show()

    def on_key(self, event):
        if event.key in ('left', 'backspace'):
            self.prev()
        elif event.key in ('right', ' '):
            self.next()

    def refresh(self):
        # keep the window responsive while the terminal waits for input
        if self.fig is not None:
            try:
                self.fig.canvas.flush_events()
            except Exception:
                pass

    def close(self):
        if self.fig is not None:
            try:
                self.plt.close(self.fig)
            except Exception:
                pass
            self.fig = None

#############################################################################
# prompts
#############################################################################

def prompt_key(prompt, keys):
    """One keypress from keys ('\n' for Enter), without waiting for Enter;
    a line of input where stdin is not a terminal."""
    if not sys.stdin.isatty():
        val = prompt_valid(prompt, [k.strip() for k in keys])
        return '\n' if val == '' else val
    sys.stdout.write(prompt)
    sys.stdout.flush()
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    try:
        tty.setcbreak(fd)
        while True:
            key = sys.stdin.read(1)
            if key == '\r':
                key = '\n'
            if key in ('', '\x04'):
                raise KeyboardInterrupt
            if key in keys:
                break
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old)
    print('' if key == '\n' else key)
    return key


class SkipReview(Exception):
    """'s' at a scoring prompt: discard the scores and move on."""


def prompt_valid(prompt, valid, images=None, skip=False):
    while True:
        if images is not None:
            images.refresh()
        try:
            val = input(prompt).strip()
        except EOFError:
            raise KeyboardInterrupt
        if val in valid:
            return val
        if skip and val == 's':
            raise SkipReview
        print('  invalid input, please enter one of: {0}'.format(
            ', '.join(valid)))


def prompt_notes(required, images=None):
    while True:
        if images is not None:
            images.refresh()
        label = 'required' if required else 'optional'
        try:
            val = input('  notes ({0}): '.format(label)).strip()
        except EOFError:
            raise KeyboardInterrupt
        if val == 's':
            raise SkipReview
        val = re.sub(r'[\t\n]', ' ', val)
        if val:
            try:
                float(val)
            except ValueError:
                return val
            print('  numeric notes are invalid')
            continue
        if required:
            print('  notes are required when QC=0 with no artifact scored '
                  '{0}'.format(SEVERE_SCORE))
            continue
        return val

#############################################################################
# locks
#############################################################################

class LockError(Exception):
    """The save lock could not be taken."""


def lock_content(reviewer_id):
    """What a lock file holds: who, where and when, and a token that makes
    each lock unique, so a holder can tell its own lock from a newer one."""
    now = time.time()
    return ('reviewer_id\t{0}\nhost\t{1}\npid\t{2}\nstarted\t{3}\n'
            'epoch\t{4:.3f}\ntoken\t{5}\n').format(
        reviewer_id, socket.gethostname(), os.getpid(),
        fmt_time(datetime.datetime.fromtimestamp(now)), now,
        os.urandom(8).hex())


def read_lock(path):
    """The lock file's content, or None if there is none."""
    try:
        with open(path) as f:
            return f.read()
    except (IOError, OSError):
        return None


def lock_info(content):
    info = {}
    for line in (content or '').splitlines():
        key, _, val = line.partition('\t')
        info[key] = val
    return info


def describe_lock(content):
    info = lock_info(content)
    return '{0} on {1} (pid {2}) since {3}'.format(
        info.get('reviewer_id', '?'), info.get('host', '?'),
        info.get('pid', '?'), info.get('started', '?'))


def lock_age(path, content):
    """Seconds since the lock was taken, from its epoch, else its mtime."""
    try:
        return time.time() - float(lock_info(content)['epoch'])
    except (KeyError, ValueError):
        try:
            return time.time() - os.path.getmtime(path)
        except OSError:
            return 0


def create_lock(path, content):
    """Create the lock file atomically; False if it exists."""
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o664)
    except FileExistsError:
        return False
    with os.fdopen(fd, 'w') as f:
        f.write(content)
    return True


def break_lock(path, seen):
    """Remove the lock file if it still holds seen: renamed aside first, so
    only one of two reviewers breaking it at once removes it, and a newer
    lock renamed by mistake is linked back."""
    aside = '{0}.{1}.{2}'.format(path, socket.gethostname(), os.getpid())
    try:
        os.rename(path, aside)
    except OSError:
        return
    if read_lock(aside) != seen:
        try:
            os.link(aside, path)
        except OSError:
            pass
    try:
        os.remove(aside)
    except OSError:
        pass


def release_lock(path, content):
    """Remove the lock file only if it is still ours."""
    if content and read_lock(path) == content:
        try:
            os.remove(path)
        except OSError:
            pass


def take_lock(path, content, max_age, override=False):
    """(True, '') once the lock file holds content, else (False, the
    holder's content); a lock older than max_age seconds, or any with
    override, is broken and replaced."""
    for _ in range(5):
        if create_lock(path, content):
            return True, ''
        seen = read_lock(path)
        if seen is None:
            continue
        if not override and lock_age(path, seen) < max_age:
            return False, seen
        break_lock(path, seen)
    return False, read_lock(path) or ''


def acquire_save_lock(path, reviewer_id, wait=LOCK_WAIT):
    """The save lock's content once taken, retrying for wait seconds."""
    content = lock_content(reviewer_id)
    deadline = time.time() + wait
    while True:
        ok, holder = take_lock(path, content, LOCK_STALE)
        if ok:
            return content
        if time.time() > deadline:
            raise LockError('could not take {0} within {1} s; held by {2}'.format(
                path, wait, describe_lock(holder)))
        time.sleep(0.05 + 0.2 * random.random())


def lock_paths(args, pid, sid):
    """The save lock of the session and the review lock of the modality."""
    outdir = os.path.join(args.output_rootdir, pid, sid)
    stem = os.path.join(outdir, '{0}_{1}_desc-manualpostqc'.format(pid, sid))
    return stem + '_save.lock', '{0}_{1}_review.lock'.format(stem, args.modality)


def take_review_lock(args, pid, sid):
    """(content, '') once this reviewer holds the review lock, else
    (None, the holder's content)."""
    _, review_lock = lock_paths(args, pid, sid)
    os.makedirs(os.path.dirname(review_lock), exist_ok=True)
    content = lock_content(args.reviewer_id)
    ok, holder = take_lock(review_lock, content, args.lock_hours * 3600,
                           args.override_lock)
    return (content, '') if ok else (None, holder)

#############################################################################
# output
#############################################################################

def output_stem(output_rootdir, pid, sid):
    """The session's file stem, {pid}_{sid} in its output directory."""
    return os.path.join(output_rootdir, pid, sid, '{0}_{1}'.format(pid, sid))


def output_paths(args, pid, sid):
    """outdir, the file stem, and the session and status files."""
    stem = output_stem(args.output_rootdir, pid, sid)
    return (os.path.dirname(stem), stem,
            stem + '_desc-session_manualpostqc.tsv',
            stem + '_desc-manualpostqc_status.tsv')


def detail_path(stem, modality):
    """The per-reviewer file of a modality."""
    return '{0}_desc-{1}_manualpostqc.tsv'.format(stem, modality)


def old_modality_path(stem, modality):
    """The per-reviewer file of a modality under its older name."""
    return '{0}_desc-manualpostqc{1}_qc.tsv'.format(stem, modality)


def old_detail_path(stem):
    """The single per-reviewer file written before the per-modality ones."""
    return stem + '_desc-manualpostqc_qc.tsv'


def old_paths(stem):
    """The files of the older names that exist, removed once converted."""
    fnames = [old_modality_path(stem, mod) for mod in MODALITIES] + \
        [old_detail_path(stem), stem + '_desc-manualpostqc_session.tsv']
    return [f for f in fnames if os.path.isfile(f)]


def review_names(modality):
    """Short key -> release name in the modality's per-reviewer file,
    table mr_y_qc__post__man__<modality>__review."""
    t = PREFIX + modality + '__review'
    names = OrderedDict([('QC', t + '_score')])
    for ftype in FTYPES[modality]:
        names['QU_' + ftype] = '{0}__{1}_score'.format(t, FTYPE_COLUMN[ftype])
    names['notes'] = t + '_notes'
    names['started'] = t + '__start_dtt'
    names['completed'] = t + '__end_dtt'
    names['duration'] = t + '_dur'
    return names


def read_tsv_rows(fname):
    """The rows of a table, with n/a read as empty (files written before
    n/a have empty cells)."""
    if not os.path.isfile(fname):
        return []
    with open(fname, newline='') as f:
        return [dict((k, '' if isinstance(v, str) and v.strip() == NA else v)
                     for k, v in r.items())
                for r in csv.DictReader(f, delimiter='\t')]


def read_review_rows(stem):
    """Every review of the session, in the short keys with modality: the
    per-modality files' (of the older name where the current one is
    missing), then those of an old single file that the per-modality files
    do not have."""
    rows = []
    for mod in MODALITIES:
        keys = dict((v, k) for k, v in review_names(mod).items())
        fname = detail_path(stem, mod)
        if not os.path.isfile(fname):
            fname = old_modality_path(stem, mod)
        for r in read_tsv_rows(fname):
            row = dict((keys.get(k, k), v) for k, v in r.items())
            row['modality'] = mod
            rows.append(row)
    have = set((r['modality'], r.get('reviewer_id')) for r in rows)
    for r in read_tsv_rows(old_detail_path(stem)):
        if r.get('modality') in MODALITIES and \
                (r['modality'], r.get('reviewer_id')) not in have:
            rows.append(r)
    return rows


def write_review_rows(stem, modality, rows):
    """The modality's rows to its per-reviewer file."""
    names = review_names(modality)
    out = []
    for r in rows:
        if r.get('modality') != modality:
            continue
        o = dict((c, r.get(c)) for c in REVIEW_KEYS)
        for k, v in names.items():
            o[v] = r.get(k)
        out.append(o)
    write_tsv(detail_path(stem, modality),
              REVIEW_KEYS + list(names.values()), out)


def write_tsv(fname, columns, rows, na=NA):
    """rows to a table, empty values as na."""
    with open(fname, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(
                (c, na if row.get(c) is None or str(row.get(c)).strip() == ''
                 else row.get(c)) for c in columns))


def to_num(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def round_half_up(x):
    return int(x + 0.5) if x >= 0 else -int(-x + 0.5)


def consensus_qc(row, mod):
    """A reviewer's QC for the consensus: a pass with an artifact scored
    SEVERE_SCORE or more counts as a fail."""
    qc = to_num(row.get('QC'))
    if qc == 1 and any((to_num(row.get('QU_' + f)) or 0) >= SEVERE_SCORE
                       for f in FTYPES[mod]):
        return 0
    return qc


def compute_consensus(pid, sid, rows):
    """The session row in the released names, the reviewers of each
    modality combined as MMIL_Summarize_postQC does: rev_count the reviews
    with a QC; rev__disp_indicator 1 where any two QC differ; _score the
    mean QC less 0.01, rounded, so a split is a fail; each artifact score
    the mean, a missing one as 0; notes 'R1: ... R2: ...'. A reviewer's QC
    1 with an artifact scored 3 counts as 0 (consensus_qc), in _score and
    rev__disp_indicator."""
    out = {'participant_id': pid, 'session_id': sid}
    columns = ['participant_id', 'session_id']
    for mod in MODALITIES:
        mrows = [r for r in rows if r.get('modality') == mod and
                 to_num(r.get('QC')) is not None]
        if not mrows:
            continue
        p = PREFIX + mod
        qcs = [consensus_qc(r, mod) for r in mrows]
        nrev = len(qcs)
        out[p + '_score'] = round_half_up(sum(qcs) / nrev - 0.01)
        out[p + '__rev_count'] = nrev
        out[p + '__rev__disp_indicator'] = int(len(set(qcs)) > 1)
        columns += [p + '_score', p + '__rev_count', p + '__rev__disp_indicator']
        for ftype in FTYPES[mod]:
            col = '{0}__{1}_score'.format(p, FTYPE_COLUMN[ftype])
            vals = [to_num(r.get('QU_' + ftype)) for r in mrows]
            out[col] = round(sum(v or 0 for v in vals) / nrev, 4)
            columns.append(col)
        notes = ['R{0}: {1}'.format(k + 1, r['notes'])
                 for k, r in enumerate(mrows) if r.get('notes')]
        out[p + '_notes'] = ' '.join(notes)
        columns.append(p + '_notes')
    return columns, out


def write_status(status_path, pid, sid, errcode, errmsg):
    now = datetime.datetime.now()
    procdate = now.strftime('%Y-%m-%dT%H:%M:%S.') + \
        '{0:03d}'.format(now.microsecond // 1000)
    write_tsv(status_path, STATUS_COLUMNS,
              [{'participant_id': pid, 'session_id': sid, 'procdate': procdate,
                'errcode': errcode, 'errmsg': errmsg or ''}], na='')
    cmig_sing_utils.add_status_versions(status_path, app_name=APPNAME)


def write_review(args, pid, sid, review):
    """Replace this reviewer's row for the modality, recompute the
    consensus and move the status procdate."""
    outdir, stem, session_path, status_path = output_paths(args, pid, sid)
    os.makedirs(outdir, exist_ok=True)
    save_lock, _ = lock_paths(args, pid, sid)
    lock = acquire_save_lock(save_lock, args.reviewer_id)
    try:
        rows = [r for r in read_review_rows(stem)
                if not (r.get('modality') == args.modality and
                        r.get('reviewer_id') == args.reviewer_id)]
        # as it reads back: a value of n/a is empty
        rows.append(dict((k, '' if isinstance(v, str) and v.strip() == NA
                          else v) for k, v in review.items()))
        # the reviews in files of the older names go to the current files;
        #   the old files are removed once those are written
        old_fnames = old_paths(stem)
        for mod in MODALITIES:
            if mod == args.modality or \
                    (old_fnames and any(r.get('modality') == mod for r in rows)):
                write_review_rows(stem, mod, rows)
        columns, consensus = compute_consensus(pid, sid, rows)
        write_tsv(session_path, columns, [consensus])
        for fname in old_fnames:
            os.remove(fname)
        write_status(status_path, pid, sid, 0, '')
    finally:
        release_lock(save_lock, lock)
    return detail_path(stem, args.modality)

#############################################################################
# review
#############################################################################

def fmt_time(dt):
    return dt.strftime('%Y-%m-%dT%H:%M:%S')


def already_reviewed(args, pid, sid):
    """Why the session is skipped, or '' to review it."""
    if args.force:
        return ''
    _, stem, _, _ = output_paths(args, pid, sid)
    revs = set(r.get('reviewer_id') for r in read_review_rows(stem)
               if r.get('modality') == args.modality)
    if args.reviewer_id in revs:
        return 'already reviewed by {0}'.format(args.reviewer_id)
    if args.max_nrev > 0 and len(revs) >= args.max_nrev:
        return 'max reviewers ({0}) reached'.format(args.max_nrev)
    return ''


def recheck_reviewed(args, pid, sid):
    """already_reviewed under the save lock, so no save is half written."""
    save_lock, _ = lock_paths(args, pid, sid)
    lock = acquire_save_lock(save_lock, args.reviewer_id)
    try:
        return already_reviewed(args, pid, sid)
    finally:
        release_lock(save_lock, lock)


def review_session(args, pid, sid):
    subj = '{0}_{1}'.format(pid, sid)
    why = already_reviewed(args, pid, sid)
    if why:
        if args.verbose:
            print('{0}: {1}, skipping'.format(subj, why))
        return
    inputs, errmsg = get_inputs(args, pid, sid)
    if inputs is None:
        print('{0}: WARNING: skipping {1}: {2}'.format(APPNAME, subj, errmsg))
        return
    lock, holder = take_review_lock(args, pid, sid)
    if lock is None:
        print('{0} ({1}): under review by {2}, skipping'.format(
            subj, args.modality, describe_lock(holder)))
        return
    _, review_lock = lock_paths(args, pid, sid)
    try:
        # again under the locks, as another reviewer may have saved meanwhile
        try:
            why = recheck_reviewed(args, pid, sid)
        except LockError as e:
            why = str(e)
        if why:
            if args.verbose:
                print('{0}: {1}, skipping'.format(subj, why))
            return
        run_review(args, pid, sid, subj, inputs, review_lock, lock)
    finally:
        release_lock(review_lock, lock)


def run_review(args, pid, sid, subj, inputs, review_lock, lock):
    """The prompts and viewers for one session and modality, then the save,
    with the review lock held."""
    print('\n{0} ({1})'.format(subj, args.modality))
    if prompt_key('press Enter to begin, or s to skip... ',
                  ['\n', 's']) == 's':
        print('skipping...')
        return
    print('  s at any prompt skips, discarding the scores')
    images = get_images(args, pid, sid) if args.images else []
    while True:
        started = datetime.datetime.now()
        proc = launch_viewer(args, pid, sid, inputs)
        viewer = None
        image_proc = None
        if images and interactive_backend():
            viewer = ImageViewer(images, title=subj)
            try:
                viewer.open()
            except Exception as e:
                print('  WARNING: cannot show the postqc images: {0}'.format(e))
                viewer = None
        if images and viewer is None:
            image_proc = launch_image_viewer(args, images)
        scores = {}
        try:
            for ftype in FTYPES[args.modality]:
                scores[ftype] = int(prompt_valid(
                    '  rate {0}? [0|1|2|3]: '.format(FTYPE_LABEL[ftype]),
                    ['0', '1', '2', '3'], viewer, skip=True))
            severe = [f for f in FTYPES[args.modality]
                      if scores[f] >= SEVERE_SCORE]
            if severe:
                # a severe artifact fails the review; the verdict is not
                #   asked, and the score stands for the reason
                qc = 0
                print('  -> failed automatically: {0} scored {1}'.format(
                    ', '.join(FTYPE_LABEL[f] for f in severe), SEVERE_SCORE))
            else:
                qc = int(prompt_valid('  accept? [1|0]: ', ['1', '0'],
                                      viewer, skip=True))
            notes = prompt_notes(qc == 0 and not severe, viewer)
            print('\n--- summary ---')
            for ftype in FTYPES[args.modality]:
                print('  {0} = {1}'.format(FTYPE_LABEL[ftype], scores[ftype]))
            print('  qc = {0}'.format(qc))
            if notes:
                print('  notes = {0}'.format(notes))
            choice = prompt_valid('confirm, repeat, or skip? [c|r|s]: ',
                                  ['c', 'r', 's'], viewer)
        except SkipReview:
            choice = 's'
        if viewer is not None:
            viewer.close()
        for p in (proc, image_proc):
            if p is None:
                continue
            try:
                p.terminate()
                p.wait(timeout=5)
            except Exception:
                pass
        completed = datetime.datetime.now()
        if choice == 'r':
            print('repeating...')
            continue
        if choice == 's':
            print('skipping...')
            return
        review = {'participant_id': pid, 'session_id': sid,
                  'modality': args.modality,
                  'reviewer_id': args.reviewer_id, 'QC': qc,
                  'notes': notes, 'started': fmt_time(started),
                  'completed': fmt_time(completed),
                  'duration': round((completed - started).total_seconds(), 1)}
        for ftype, val in scores.items():
            review['QU_' + ftype] = val
        if read_lock(review_lock) != lock:
            print('  WARNING: the review lock was taken over by {0}'.format(
                describe_lock(read_lock(review_lock))))
        while True:
            try:
                fpath = write_review(args, pid, sid, review)
                break
            except LockError as e:
                print('{0}: ERROR: not saved: {1}'.format(APPNAME, e))
                if prompt_valid('retry the save? [y|n]: ',
                                ['y', 'n']) == 'n':
                    return
        if args.verbose:
            print('{0}: saved {1}'.format(subj, fpath))
        return

#############################################################################
# CLI
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Manual post-processing QC review of BIDS derivatives')
    parser.add_argument('--modality', nargs='+', default=MODALITIES,
                        choices=MODALITIES,
                        help='what to review, in turn for each session: any '
                             'of t2, fmri, dmri and fsurf (default: all four, '
                             'those whose input roots are given)')
    parser.add_argument('--output_rootdir', required=True,
                        help='root directory of the review files')
    parser.add_argument('--reviewer_id', required=True,
                        help='reviewer initials (e.g. DH)')
    parser.add_argument('--proc_rootdir',
                        help='cmig_proc root: the T1w and the T2w (all but '
                             'fsurf)')
    parser.add_argument('--postqc_rootdir',
                        help='cmig_postqc root: the fMRI mean image, and the '
                             'images shown beside the viewer')
    parser.add_argument('--analysis_rootdir',
                        help='cmig_analysis root: the T1w-space RSI RND and '
                             'AtlasTrack dseg (dmri)')
    parser.add_argument('--fsurf_rootdir',
                        help='FreeSurfer root, {pid}_{sid} recons (fsurf)')
    parser.add_argument('--selected_file',
                        help='TSV with participant_id and session_id '
                             'columns and, optionally, modality: each session '
                             'is then reviewed only for the modalities listed '
                             'for it; and, optionally, pipeline_rootdir: the '
                             'cmig_pipeline root of the session\'s inputs, in '
                             'place of the root options where given; other '
                             'columns are not read')
    parser.add_argument('--participant_id',
                        help='review this participant; alone, all its sessions')
    parser.add_argument('--session_id',
                        help='review this session; alone, every participant '
                             'with it')
    parser.add_argument('--random_order', type=int, default=1, choices=[0, 1],
                        help='review the sessions in random order (default: '
                             '1); 0 for the order of --selected_file, or '
                             'sorted')
    parser.add_argument('--random_seed', type=int,
                        help='seed of the random order (default: drawn, '
                             'printed and written in the qc tsv)')
    parser.add_argument('--viewer', default='tkmedit',
                        choices=['tkmedit', 'freeview', 'none'],
                        help='volume viewer (default: tkmedit); none for no '
                             'viewer')
    parser.add_argument('--images', type=int, default=0, choices=[0, 1],
                        help='show the cmig_postqc_bids_* images (default: 0)')
    parser.add_argument('--image_viewer', default='eog',
                        help='program that shows the images where matplotlib '
                             'cannot open a window (default: eog)')
    parser.add_argument('--tcl_dir',
                        help='directory of the tcl scripts (default: '
                             '$CMIG_PIPELINE_DIR/tcl)')
    parser.add_argument('--max_nrev', type=int, default=2,
                        help='maximum number of reviewers per session '
                             '(default: 2); 0 = no limit')
    parser.add_argument('--force', type=int, default=0, choices=[0, 1],
                        help='review sessions already reviewed (default: 0)')
    parser.add_argument('--lock_hours', type=float, default=4,
                        help='hours after which another reviewer\'s review '
                             'lock expires and is taken over (default: 4)')
    parser.add_argument('--override_lock', type=int, default=0, choices=[0, 1],
                        help='take over review locks held by other reviewers '
                             '(default: 0)')
    parser.add_argument('--suppress', type=int, default=1, choices=[0, 1],
                        help='suppress viewer output (default: 1)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='print status messages (default: 0)')
    return parser.parse_args()


# the input roots each modality needs
NEEDS = {'t2': ['proc_rootdir'], 'fmri': ['proc_rootdir', 'postqc_rootdir'],
         'dmri': ['proc_rootdir', 'analysis_rootdir'],
         'fsurf': ['fsurf_rootdir']}


def for_modality(args, modality, pipeline_rootdir=None):
    """args with modality set to one of the requested ones, and, with
    pipeline_rootdir, the input roots under it."""
    margs = copy.copy(args)
    margs.modality = modality
    if pipeline_rootdir:
        for name, subdir in PIPELINE_SUBDIRS.items():
            setattr(margs, name, os.path.join(pipeline_rootdir, subdir))
    return margs


def check_args(args, roots_in_file=False):
    """(modalities, errmsg): the requested modalities whose input roots are
    given; a missing root is an error for a modality named explicitly and
    drops it, with a warning, from the default list. With roots_in_file
    (a pipeline_rootdir column) every requested modality is kept, and a
    session lacking both is skipped at review."""
    if args.selected_file and (args.participant_id or args.session_id):
        return [], '--selected_file excludes --participant_id and --session_id'
    if roots_in_file:
        return [m for m in MODALITIES if m in args.modality], ''
    explicit = '--modality' in sys.argv
    modalities = []
    for mod in MODALITIES:
        if mod not in args.modality:
            continue
        missing = ['--' + n for n in NEEDS[mod] if not getattr(args, n)]
        if not missing:
            modalities.append(mod)
        elif explicit:
            return [], '--modality {0} needs {1}'.format(
                mod, ' and '.join(missing))
        elif args.verbose:
            print('{0}: no {1}, so not reviewing {2}'.format(
                APPNAME, ' or '.join(missing), mod))
    if not modalities:
        return [], 'no modality has its input roots'
    return modalities, ''


def session_order(args, sessions):
    """The sessions shuffled with args.random_seed (drawn if not given), or
    as they are; sets args.random_seed ('' when not random)."""
    sessions = list(sessions)
    if not args.random_order:
        args.random_seed = ''
        return sessions
    if args.random_seed is None:
        args.random_seed = random.SystemRandom().randint(1, 999999)
    random.Random(args.random_seed).shuffle(sessions)
    print('{0}: random order, seed {1}'.format(APPNAME, args.random_seed))
    return sessions


def run(args):
    # each session with any requested modality, reviewed for each in turn
    selected = None
    pipeline_roots = None
    if args.selected_file:
        sessions, selected, pipeline_roots, errmsg = read_selected_file(
            args.selected_file)
        if errmsg:
            return 1, errmsg
    modalities, errmsg = check_args(args, pipeline_roots is not None)
    if errmsg:
        return 1, errmsg
    if not args.selected_file:
        sessions = set()
        for mod in modalities:
            sessions.update(get_sessions(for_modality(args, mod)))
        sessions = sorted(sessions)
    todo = []
    for pid, sid in sessions:
        mods = [m for m in modalities
                if selected is None or m in selected.get((pid, sid), ())]
        if mods:
            todo.append(((pid, sid), mods))
    if not todo:
        return 1, 'no sessions found'
    if args.verbose:
        print('{0}: {1} sessions to review for {2}'.format(
            APPNAME, len(todo), ', '.join(modalities)))
    todo = session_order(args, todo)
    try:
        for (pid, sid), mods in todo:
            root = (pipeline_roots or {}).get((pid, sid))
            for mod in mods:
                margs = for_modality(args, mod, root)
                missing = ['--' + n for n in NEEDS[mod]
                           if not getattr(margs, n)]
                if missing:
                    print('{0}: WARNING: skipping {1}_{2} ({3}): no '
                          'pipeline_rootdir and no {4}'.format(
                              APPNAME, pid, sid, mod, ' or '.join(missing)))
                    continue
                review_session(margs, pid, sid)
    except KeyboardInterrupt:
        print('\n{0}: stopped'.format(APPNAME))
    return 0, ''


def main():
    args = parse_args()
    errcode, errmsg = run(args)
    if errcode:
        print('{0}: ERROR: {1}'.format(APPNAME, errmsg))
    sys.exit(errcode)


if __name__ == '__main__':
    main()
