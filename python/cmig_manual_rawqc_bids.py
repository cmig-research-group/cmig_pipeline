#!/usr/bin/env python3
"""
Manual raw QC review of BIDS data, using the images from cmig_rawqc_bids.

Presents the pre-rendered images for each scan, takes a pass/fail and a
severity per artifact, and writes the review into the same per-session files
cmig_copy_manual_rawqc_bids_abcd.py writes, marked source=reviewer.

A review takes a lock on each modality of the session it covers
({pid}_{sid}_desc-manualrawqc_{modality}_review.lock), so a second reviewer
skips those meanwhile; it expires after --lock_hours and --override_lock 1
takes it regardless. Each save holds {pid}_{sid}_desc-manualrawqc_save.lock
for its read-replace-write. Lock files are created with O_EXCL, not flock,
since reviewers share NFS across hosts.

Created:  09/21/2026 by Don Hagler
Prev mod: 09/21/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Last mod: 11/12/2026 by Don Hagler
"""

import argparse
import csv
import datetime
import glob
import os
import random
import re
import socket
import sys
import termios
import time
import tty

import cmig_copy_manual_rawqc_bids_abcd as copyqc
import cmig_sing_utils

APPNAME = 'cmig_manual_rawqc_bids'

# artifact -> the column it is written to; the shared schema's map, so the
# two scripts cannot disagree about where a score goes
QU_COLUMN_MAP = copyqc.QU_COLUMN_MAP

# what the reviewer is asked to rate, per modality: tiers 1 and 2, to be
# narrowed after the testing round
FTYPE_LISTS = {
    'sMRI': ['mo', 'sus'],
    'dMRI': ['sus', 'fov', 'line', 'fa', 'gh'],
    'fMRI': ['sus', 'fov', 'line', 'gh', 'moire'],
}

# how each artifact is put to the reviewer
FTYPE_LABEL = {
    'mo': 'motion', 'sus': 'susceptibility', 'fov': 'cutoff',
    'line': 'line artifact', 'fa': 'FA artifact', 'gh': 'ghosting',
    'moire': 'moire',
}

# the severity at which the scan fails regardless of the verdict given
SEVERE_SCORE = copyqc.SEVERE_SCORE

MODALITIES = ['anat', 'dwi', 'fmap', 'func']

LOCK_WAIT = 60     # seconds a save waits for the save lock
LOCK_STALE = 300   # seconds after which a save lock is broken

# the images shown first for each modality, in order; anything else the scan
# has is shown after these
INFIX_ORDER = ['desc-3view', 'desc-sagMS', 'desc-corMS', 'desc-horMS',
               'desc-sagMF', 'desc-corMF', 'desc-horMF', 'desc-pepolar',
               'fovcutoff', 'lines', 'motion']

#############################################################################
# series types
#############################################################################

def get_series_type(filename):
    # derive the csv SeriesType from the BIDS filename, matching the values
    # cmig_copy_manual_rawqc_bids_abcd.py writes
    base = os.path.basename(filename)
    stem = re.sub(r'\.nii(\.gz)?$', '', base)
    parts = stem.split('_')
    suffix = parts[-1]
    ents = {}
    for part in parts[:-1]:
        if '-' in part:
            key, val = part.split('-', 1)
            ents[key] = val
    norm = (ents.get('rec') or '').lower() == 'norm'
    if suffix == 'T1w':
        return 'T1_NORM' if norm else 'T1'
    if suffix == 'T2w':
        return 'T2_NORM' if norm else 'T2'
    if suffix == 'dwi':
        return 'dMRI'
    if suffix == 'bold':
        task = (ents.get('task') or '').lower()
        return {'rest':  'rsfMRI',
                'mid':   'fMRI_MID_task',
                'sst':   'fMRI_SST_task',
                'nback': 'fMRI_nBack_task'}.get(task)
    if suffix == 'epi':
        acq = (ents.get('acq') or '').lower()
        direc = (ents.get('dir') or '').upper()
        if acq == 'func' and direc in ('AP', 'PA'):
            return 'fMRI_FM_' + direc
        if acq == 'dwi' and direc in ('AP', 'PA'):
            return 'dMRI_FM_' + direc
    return None

#############################################################################
# sessions and scans
#############################################################################

def get_sessions(rawqc_rootdir, participant_id=None, session_id=None,
                 selected_file=None):
    # sessions to review, from a selection file or the rawqc tree
    if selected_file:
        sessions = []
        with open(selected_file, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                pid = (row.get('participant_id') or '').strip()
                sid = (row.get('session_id') or '').strip()
                if pid and sid:
                    sessions.append((pid, sid))
        return sessions
    pid_pats = participant_id if participant_id else ['sub-*']
    sid_pat = session_id if session_id else 'ses-*'
    dirs = []
    for pid_pat in pid_pats:
        dirs += glob.glob(os.path.join(rawqc_rootdir, pid_pat, sid_pat))
    sessions = []
    for d in sorted(set(dirs)):
        if os.path.isdir(d):
            sessions.append((os.path.basename(os.path.dirname(d)),
                             os.path.basename(d)))
    return sessions


def get_scans(rawqc_rootdir, pid, sid, modalities, verbose=0):
    # one entry per scan, read from the per-scan rawqc tsvs
    scans = []
    for modality in modalities:
        moddir = os.path.join(rawqc_rootdir, pid, sid, modality)
        if not os.path.isdir(moddir):
            continue
        for tsv in sorted(glob.glob(os.path.join(moddir, '*_rawqc.tsv'))):
            if re.search(r'desc-(session|rawqc_status)', os.path.basename(tsv)):
                continue
            try:
                with open(tsv, newline='') as f:
                    rows = [r for r in csv.DictReader(f, delimiter='\t')]
            except (IOError, OSError):
                print('{0}: WARNING: failed to read {1}'.format(APPNAME, tsv))
                continue
            for row in rows:
                fname = strip_quotes(row.get('filename'))
                if not fname:
                    continue
                stype = get_series_type(fname)
                if stype is None:
                    print('{0}: WARNING: unrecognized series for {1}'.format(
                        APPNAME, fname))
                    continue
                scans.append({
                    'filename':    fname,
                    'series_type': stype,
                    'modality':    modality,
                    'moddir':      moddir,
                    'acq_time':    strip_quotes(row.get('acq_time')),
                    'series_uid':  get_series_uid(row),
                })
    scans.sort(key=lambda s: (s['acq_time'] or '', s['filename']))
    return scans


def strip_quotes(val):
    if val is None:
        return None
    return val.strip().strip('"').strip()


def get_series_uid(row):
    # the rawqc tsvs carry the uid under a metric-style column name
    for col, val in row.items():
        if col and col.endswith('series_uid'):
            uid = strip_quotes(val)
            if uid and uid != copyqc.NA:
                return uid
    return None


def get_images(scan):
    # the rendered images for this scan: same entity prefix, and for anat the
    # suffix too, since T1w and T2w share a prefix
    base = os.path.basename(scan['filename'])
    stem = re.sub(r'\.nii(\.gz)?$', '', base)
    parts = stem.split('_')
    suffix = parts[-1]
    prefix = '_'.join(parts[:-1])
    flist = glob.glob(os.path.join(scan['moddir'], prefix + '_*.jpg'))
    if scan['modality'] == 'anat':
        flist = [f for f in flist if suffix in os.path.basename(f)]
    return sorted(flist, key=image_sort_key)


def image_sort_key(fpath):
    name = os.path.basename(fpath)
    for i, infix in enumerate(INFIX_ORDER):
        if infix in name:
            return (i, name)
    return (len(INFIX_ORDER), name)

#############################################################################
# image viewer
#############################################################################

class ImageViewer(object):
    """Show a list of images one at a time, with prev/next and zoom."""

    def __init__(self, fnames, title=''):
        self.fnames = fnames
        self.index = 0
        self.fig = None
        self.title = title

    def open(self):
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Button
        self.plt = plt
        self.fig = plt.figure(figsize=(11, 8))
        self.ax = self.fig.add_axes([0.02, 0.10, 0.96, 0.84])
        self.ax.set_xticks([])
        self.ax.set_yticks([])
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


class SkipScan(Exception):
    """'s' at a scan's prompt: discard its scores and move to the next."""


def prompt_valid(prompt, valid, viewer=None, skip=False):
    while True:
        if viewer is not None:
            viewer.refresh()
        try:
            val = input(prompt).strip()
        except EOFError:
            raise KeyboardInterrupt
        if val in valid:
            return val
        if skip and val == 's':
            raise SkipScan
        print('  invalid input, please enter one of: {0}'.format(
            ', '.join(valid)))


def prompt_qu(ftype, viewer=None):
    label = FTYPE_LABEL.get(ftype, ftype)
    return int(prompt_valid('  rate {0}? [0|1|2|3]: '.format(label),
                            ['0', '1', '2', '3'], viewer, skip=True))


def prompt_qc(viewer=None):
    return int(prompt_valid('  accept? [1|0]: ', ['1', '0'], viewer,
                            skip=True))


def prompt_notes(required, viewer=None):
    # notes are required where the reviewer failed a scan with no artifact
    # scored 3, so the reason for the fail is recorded
    while True:
        if viewer is not None:
            viewer.refresh()
        label = 'required' if required else 'optional'
        try:
            val = input('  notes ({0}): '.format(label)).strip()
        except EOFError:
            raise KeyboardInterrupt
        if val == 's':
            raise SkipScan
        val = re.sub(r',', ';', val)
        val = re.sub(r'[\'"]', '', val)
        if val:
            try:
                float(val)
            except ValueError:
                return val
            print('  numeric notes are invalid')
            continue
        if required:
            print('  notes are required when a scan fails with no artifact '
                  'scored {0}'.format(SEVERE_SCORE))
            continue
        return val

#############################################################################
# review
#############################################################################

def review_scan(scan, ftypes, verbose=0):
    # returns a review dict, or None if the scan was skipped
    images = get_images(scan)
    if not images:
        print('{0}: WARNING: no images for {1}, skipping (was rawqc run?)'
              .format(APPNAME, scan['filename']))
        return None
    print('')
    print('{0}  [{1}]  {2} images'.format(
        scan['filename'], scan['series_type'], len(images)))
    viewer = ImageViewer(images, title=scan['series_type'])
    try:
        viewer.open()
    except Exception as e:
        print('{0}: WARNING: could not open the image viewer ({1})'.format(
            APPNAME, e))
        viewer = None
    try:
        scores = {}
        for ftype in ftypes:
            scores[ftype] = prompt_qu(ftype, viewer)
        severe = [f for f in ftypes if scores[f] >= SEVERE_SCORE]
        if severe:
            # a scan carrying a severe artifact does not pass; the severity
            # governs and the verdict is not asked for
            qc = 0
            print('  -> failed automatically: {0} scored {1}'.format(
                ', '.join(FTYPE_LABEL.get(f, f) for f in severe),
                SEVERE_SCORE))
        else:
            qc = prompt_qc(viewer)
        notes = prompt_notes(qc == 0 and not severe, viewer)
    finally:
        if viewer is not None:
            viewer.close()
    return {'scan': scan, 'QC': qc, 'scores': scores, 'notes': notes}


def print_scan_summary(rev):
    scored = ' '.join(
        '{0}={1}'.format(f, v) for f, v in sorted(rev['scores'].items())
        if v > 0)
    print('  --- {0} {1}{2}{3}'.format(
        'accept' if rev['QC'] else 'REJECT',
        os.path.basename(rev['scan']['filename']),
        '  ' + scored if scored else '',
        '  "' + rev['notes'] + '"' if rev['notes'] else ''))

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


def lock_stem(output_rootdir, pid, sid):
    return os.path.join(output_rootdir, pid, sid,
                        '{0}_{1}_desc-manualrawqc'.format(pid, sid))


def save_lock_path(output_rootdir, pid, sid):
    return lock_stem(output_rootdir, pid, sid) + '_save.lock'


def review_lock_path(output_rootdir, pid, sid, modality):
    return '{0}_{1}_review.lock'.format(
        lock_stem(output_rootdir, pid, sid), modality)


def take_review_locks(args, pid, sid, modalities):
    """{modality: content} of the review locks taken, and
    {modality: holder's content} of those held by another reviewer."""
    os.makedirs(os.path.join(args.output_rootdir, pid, sid), exist_ok=True)
    taken, held = {}, {}
    for modality in modalities:
        content = lock_content(args.reviewer_id)
        ok, holder = take_lock(
            review_lock_path(args.output_rootdir, pid, sid, modality),
            content, args.lock_hours * 3600, args.override_lock)
        if ok:
            taken[modality] = content
        else:
            held[modality] = holder
    return taken, held


def release_review_locks(args, pid, sid, taken):
    for modality, content in taken.items():
        release_lock(review_lock_path(args.output_rootdir, pid, sid, modality),
                     content)

#############################################################################
# output
#############################################################################

def read_detail_rows(detail_path):
    # every existing row, normalized; the tool replaces only its own
    if not os.path.isfile(detail_path):
        return []
    rows = []
    try:
        with open(detail_path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                rows.append(copyqc.normalize_existing_row(row))
    except (IOError, OSError):
        print('{0}: WARNING: failed to read {1}'.format(APPNAME, detail_path))
    return rows


def count_reviewers(detail_path, filename):
    # how many reviewers have scored this scan, either era
    revs = set()
    for row in read_detail_rows(detail_path):
        if row.get('filename') == filename and row.get('reviewer_id'):
            revs.add(row['reviewer_id'])
    return len(revs)


def build_review_row(pid, sid, rev, reviewer_id, started, completed,
                     qu_columns):
    scan = rev['scan']
    row = {
        'participant_id': pid,
        'session_id':     sid,
        'filename':       scan['filename'],
        'series_type':    scan['series_type'],
        'reviewer_id':    reviewer_id,
        'source':         copyqc.SOURCE_REVIEWER,
        'QC':             rev['QC'],
        'qu_source':      copyqc.QU_RECORDED,
        'notes':          None if copyqc.is_empty(rev['notes']) else rev['notes'],
        'started':        started,
        'completed':      completed,
        'duration':       None,
        'revdate':        completed,
        'series_uid':     scan['series_uid'],
    }
    # an artifact that was rated is written even when clean, as 0; one that
    # was not asked for is left empty, which is what qu_source marks
    for col in qu_columns:
        row[col] = None
    for ftype, val in rev['scores'].items():
        col = QU_COLUMN_MAP.get(ftype)
        if col in qu_columns:
            row[col] = val
    return row


def write_status(status_path, pid, sid, errcode, errmsg):
    # the copy script's status columns, with this tool named as the writer
    now = datetime.datetime.now()
    procdate = now.strftime('%Y-%m-%dT%H:%M:%S.') + \
        '{0:03d}'.format(now.microsecond // 1000)
    header = ['participant_id', 'session_id', 'procdate', 'errcode', 'errmsg']
    values = [pid, sid, procdate, str(errcode), errmsg or '']
    with open(status_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        f.write('\t'.join(values) + '\n')
    # pure python with no compiled app, so there is no app_version to record
    cmig_sing_utils.add_status_versions(status_path, app_name=APPNAME)


def write_session_reviews(pid, sid, reviews, args, started, completed,
                          qu_columns):
    outdir, detail_path, summary_path, status_path = copyqc.get_output_paths(
        args.output_rootdir, pid, sid)
    os.makedirs(outdir, exist_ok=True)
    save_lock = save_lock_path(args.output_rootdir, pid, sid)
    lock = acquire_save_lock(save_lock, args.reviewer_id)
    try:
        write_locked(pid, sid, reviews, args, started, completed, qu_columns)
    finally:
        release_lock(save_lock, lock)
    return detail_path


def write_locked(pid, sid, reviews, args, started, completed, qu_columns):
    # the read-replace-write of the session files, under the save lock
    outdir, detail_path, summary_path, status_path = copyqc.get_output_paths(
        args.output_rootdir, pid, sid)
    rows = read_detail_rows(copyqc.get_existing_detail_path(
        args.output_rootdir, pid, sid))
    reviewed = set(r['scan']['filename'] for r in reviews)
    # drop this reviewer's previous rows for the scans just reviewed
    rows = [r for r in rows
            if not (r.get('reviewer_id') == args.reviewer_id and
                    r.get('source') == copyqc.SOURCE_REVIEWER and
                    r.get('filename') in reviewed)]
    for rev in reviews:
        rows.append(build_review_row(pid, sid, rev, args.reviewer_id,
                                     started, completed, qu_columns))

    order = []
    by_scan = {}
    for row in rows:
        key = (row.get('filename'), row.get('series_type'))
        if key not in by_scan:
            by_scan[key] = []
            order.append(key)
        by_scan[key].append(row)
    uids = {}
    for row in rows:
        if row.get('series_uid'):
            uids[(row.get('filename'), row.get('series_type'))] = \
                row['series_uid']
    summary_rows = [
        copyqc.compute_consensus(pid, sid, key, by_scan[key], uids.get(key))
        for key in order]

    copyqc.write_detail(detail_path, rows)
    copyqc.write_summary(summary_path, summary_rows)
    copyqc.remove_old_output(args.output_rootdir, pid, sid)
    # the status procdate is what cmig_mkjobs_bids_pipeline.py reads to
    # notice an updated review, so it moves whenever a review is written
    write_status(status_path, pid, sid, 0, '')


def select_scans(scans, detail_path, args):
    # the scans still short of max_nrev reviewers
    todo = []
    for scan in scans:
        if not args.forceflag:
            nrev = count_reviewers(detail_path, scan['filename'])
            if args.max_nrev > 0 and nrev >= args.max_nrev:
                if args.verbose:
                    print('{0}: {1} already has {2} reviewers, skipping'
                          .format(APPNAME, scan['filename'], nrev))
                continue
        todo.append(scan)
    return todo


def reselect_scans(scans, args, pid, sid):
    # select_scans under the save lock, so no save is half written
    save_lock = save_lock_path(args.output_rootdir, pid, sid)
    lock = acquire_save_lock(save_lock, args.reviewer_id)
    try:
        return select_scans(scans, copyqc.get_existing_detail_path(
            args.output_rootdir, pid, sid), args)
    finally:
        release_lock(save_lock, lock)


def review_session(pid, sid, args, qu_columns):
    scans = get_scans(args.rawqc_rootdir, pid, sid, args.modality,
                      args.verbose)
    if not scans:
        print('{0}: WARNING: no scans found for {1} {2}, skipping'.format(
            APPNAME, pid, sid))
        return 0

    detail_path = copyqc.get_existing_detail_path(args.output_rootdir, pid, sid)
    todo = select_scans(scans, detail_path, args)
    if not todo:
        if args.verbose:
            print('{0}: {1} {2}: nothing to review'.format(APPNAME, pid, sid))
        return 0

    # lock the modalities to review; those another reviewer holds are skipped
    modalities = sorted(set(s['modality'] for s in todo))
    taken, held = take_review_locks(args, pid, sid, modalities)
    try:
        for modality in sorted(held):
            print('{0}: {1} {2} {3}: under review by {4}, skipping'.format(
                APPNAME, pid, sid, modality, describe_lock(held[modality])))
        # again under the locks, as another reviewer may have saved meanwhile
        try:
            todo = reselect_scans(
                [s for s in todo if s['modality'] in taken], args, pid, sid)
        except LockError as e:
            print('{0}: {1} {2}: {3}, skipping'.format(APPNAME, pid, sid, e))
            return 0
        if not todo:
            if args.verbose:
                print('{0}: {1} {2}: nothing to review'.format(
                    APPNAME, pid, sid))
            return 0
        return review_scans(pid, sid, args, qu_columns, todo, taken)
    finally:
        release_review_locks(args, pid, sid, taken)


def review_scans(pid, sid, args, qu_columns, todo, taken):
    # the prompts and viewers for the session's scans, with the review locks
    # held; each scan is confirmed, repeated or skipped on its own, and saved
    # as soon as it is confirmed, so a stop keeps every confirmed scan
    print('')
    print('=' * 70)
    print('{0} {1}: {2} scans to review'.format(pid, sid, len(todo)))
    print('=' * 70)
    prompt_key('press Enter to begin... ', ['\n'])
    print('  s at any prompt skips this scan, discarding its scores')

    nsaved = 0
    for scan in todo:
        rev = None
        while True:
            t0 = datetime.datetime.now()
            try:
                rev = review_scan(scan, args.ftypes[
                    copyqc.modality_group(scan['series_type'])], args.verbose)
                if rev is None:
                    break
                print_scan_summary(rev)
                choice = prompt_valid('confirm, repeat, or skip? [c|r|s]: ',
                                      ['c', 'r', 's'])
            except SkipScan:
                choice = 's'
            if choice == 'r':
                print('  discarded, repeating...')
                continue
            if choice == 's':
                print('  discarded, skipping...')
                rev = None
            break
        if rev is None:
            continue
        t1 = datetime.datetime.now()
        modality = scan['modality']
        path = review_lock_path(args.output_rootdir, pid, sid, modality)
        if read_lock(path) != taken.get(modality):
            print('{0}: WARNING: the {1} review lock was taken over by '
                  '{2}'.format(APPNAME, modality,
                               describe_lock(read_lock(path))))
        saved = False
        while True:
            try:
                fpath = write_session_reviews(
                    pid, sid, [rev], args, fmt_time(t0), fmt_time(t1),
                    qu_columns)
                saved = True
                break
            except LockError as e:
                print('{0}: ERROR: not saved: {1}'.format(APPNAME, e))
                if prompt_valid('retry the save? [y|n]: ', ['y', 'n']) == 'n':
                    break
        if saved:
            nsaved += 1
            args.nsaved += 1
            if args.verbose:
                print('  saved to {0}'.format(fpath))
    print('{0}: {1} {2}: wrote {3} of {4} reviews'.format(
        APPNAME, pid, sid, nsaved, len(todo)))
    return nsaved


def fmt_time(dt):
    return dt.strftime('%Y-%m-%dT%H:%M:%S')

#############################################################################
# main
#############################################################################

def parse_args():
    parser = argparse.ArgumentParser(
        description='Manual raw QC review of BIDS data, using the images '
                    'from cmig_rawqc_bids.')
    parser.add_argument('--rawqc_rootdir', required=True,
                        help='root directory of cmig_rawqc_bids output')
    parser.add_argument('--output_rootdir', required=True,
                        help='root directory of manual rawQC output')
    parser.add_argument('--reviewer_id', required=True,
                        help='reviewer initials (e.g. DH)')
    parser.add_argument('--participant_id', nargs='+', default=None,
                        help='optional participant filter; one or more ids')
    parser.add_argument('--session_id', default=None,
                        help='optional session filter (e.g. ses-00A)')
    parser.add_argument('--selected_file', default=None,
                        help='tsv with participant_id and session_id columns')
    parser.add_argument('--modality', nargs='+', default=MODALITIES,
                        choices=MODALITIES, help='modalities to review')
    parser.add_argument('--sMRI_ftype_list', nargs='*', default=None,
                        help='artifacts rated on structural scans')
    parser.add_argument('--dMRI_ftype_list', nargs='*', default=None,
                        help='artifacts rated on diffusion scans')
    parser.add_argument('--fMRI_ftype_list', nargs='*', default=None,
                        help='artifacts rated on functional scans')
    parser.add_argument('--max_nrev', type=int, default=2,
                        help='maximum reviewers per scan; 0 = no limit')
    parser.add_argument('--rand_order_flag', type=int, default=0,
                        choices=[0, 1], help='randomize the session order')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1],
                        help='display status messages')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='review scans that already have enough reviews')
    parser.add_argument('--lock_hours', type=float, default=4,
                        help='hours after which another reviewer\'s review '
                             'lock expires and is taken over (default: 4)')
    parser.add_argument('--override_lock', type=int, default=0, choices=[0, 1],
                        help='take over review locks held by other reviewers')
    return parser.parse_args()


def check_ftypes(args):
    # the artifact lists in use, dropping any the shared schema cannot hold
    ftypes = {}
    qu_columns = list(copyqc.QU_COLUMNS)
    for group in ('sMRI', 'dMRI', 'fMRI'):
        supplied = getattr(args, '{0}_ftype_list'.format(group))
        wanted = supplied if supplied is not None else FTYPE_LISTS[group]
        keep = []
        for ftype in wanted:
            col = QU_COLUMN_MAP.get(ftype)
            if col is None:
                print('{0}: WARNING: no column for artifact {1}, dropping '
                      'it from the {2} list'.format(APPNAME, ftype, group))
                continue
            if col not in qu_columns:
                print('{0}: WARNING: {1} has no {2} column in the shared '
                      'schema, so {3} cannot be recorded; add it to '
                      'cmig_copy_manual_rawqc_bids_abcd.py to enable it'
                      .format(APPNAME, group, col, ftype))
                continue
            keep.append(ftype)
        ftypes[group] = keep
    return ftypes, qu_columns


def extend_scored_issues(ftypes):
    # compute_consensus only averages the artifacts its modality scores, so
    # an artifact this tool rates but the copy script does not parse (sMRI
    # sus) would be dropped from the summary; widen the set for this process
    for group, wanted in ftypes.items():
        issues = list(copyqc.QU_ISSUES.get(group, []))
        for ftype in wanted:
            if ftype not in issues:
                issues.append(ftype)
        copyqc.QU_ISSUES[group] = issues


def run(args):
    if not os.path.isdir(args.rawqc_rootdir):
        print('{0}: ERROR: rawqc_rootdir {1} not found'.format(
            APPNAME, args.rawqc_rootdir))
        return 1
    if not os.environ.get('DISPLAY'):
        print('{0}: WARNING: DISPLAY is not set, so no images can be shown; '
              'reconnect with ssh -X'.format(APPNAME))

    args.ftypes, qu_columns = check_ftypes(args)
    extend_scored_issues(args.ftypes)
    if args.verbose:
        for group in ('sMRI', 'dMRI', 'fMRI'):
            print('{0}: {1} artifacts: {2}'.format(
                APPNAME, group, ' '.join(args.ftypes[group]) or '(none)'))

    sessions = get_sessions(args.rawqc_rootdir, args.participant_id,
                            args.session_id, args.selected_file)
    if not sessions:
        print('{0}: ERROR: no sessions found in {1}'.format(
            APPNAME, args.rawqc_rootdir))
        return 1
    if args.rand_order_flag:
        random.shuffle(sessions)
    print('{0}: {1} sessions to review'.format(APPNAME, len(sessions)))

    # counted at each save, so a stop reports the scans saved before it
    args.nsaved = 0
    try:
        for pid, sid in sessions:
            review_session(pid, sid, args, qu_columns)
    except KeyboardInterrupt:
        print('')
        print('{0}: stopped; {1} reviews written'.format(APPNAME, args.nsaved))
        return 0
    print('{0}: {1} reviews written'.format(APPNAME, args.nsaved))
    return 0


def main():
    sys.exit(run(parse_args()))


if __name__ == '__main__':
    main()
