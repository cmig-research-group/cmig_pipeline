#!/usr/bin/env python
"""
cmig_compile_bids_tabdat.py

Compile the per-session tabulated data under one or more BIDS derivatives
trees into one TSV per table, in the layout of the release's
rawdata/phenotype/<table>.tsv.

Public function
---------------
compile_bids_tabdat(rootdirs, outdir, phenotype_dir,
                    sublist=None, seslist=None, id_events_file=None,
                    series_outdir=None, tables='all', forceflag=False,
                    verbose=True)

  With id_events_file, an id_events release csv, only the visits it lists
  (its participant_id and session_id columns) are compiled.

  Reads every sub-*/ses-*/[<modality>/]*.tsv under each rootdir, except the
  *_status.tsv and the *_corrmat.tsv (ROI-by-ROI matrices, released as
  concatenated mat files rather than tables). A task's _behavioral.tsv,
  which carries no participant_id or session_id (they are taken from its
  path), is read from func/ when there is one (selected runs with a scan),
  else from beh/; never both.

  Session-level files, one row per session, are joined on participant_id
  and session_id. Each column goes to the table whose release json lists it.
  A file no column of which the json lists is a table of its own, named for
  its columns' shared prefix, e.g. mr_y_tfmri__mid__var__gp, unless a table
  that name begins with already holds columns of its family. Otherwise a
  column no json lists goes to the longest table name it begins with,
  followed by '_' or '__', among the release tables and the new ones, a new
  table being named for the shared prefix of a file's columns that begin
  with no release table name; these are placed after every file is read, so
  a new table made from one file takes the matching columns of every other.
  The new 8.0 modality-level tables (NEW_TABLES), e.g. mr_y_qc__raw__fmri,
  take part in this routing as the release tables do, so that
  mr_y_qc__raw__fmri__fm__* goes to mr_y_qc__raw__fmri, as
  mr_y_qc__raw__dmri__fm__* goes to mr_y_qc__raw__dmri.

  A tfMRI var or gpnet corr column's __resid or __task component (from the
  desc-residual and desc-taskpred files) is left out of its table name, so
  that, e.g., mr_y_tfmri__mid__var__resid__aseg__* goes to
  mr_y_tfmri__mid__var__aseg and mr_y_tfmri__mid__corr__gpnet__aseg__task__*
  to mr_y_tfmri__mid__corr__gpnet__aseg, with the plain file's columns.

  Series-level files, one row per scan, are those carrying a filename column
  (rawqc's per-series TSVs, pc's desc-series_pc, pcqc's desc-series_pcqc,
  the manual raw QC consensus desc-series_manualrawqc). Their rows, of every
  modality, are joined on participant_id, session_id and filename into one
  table per series table name (series_table), from the file's column names:
  mr_y_qc__raw__series for rawqc, pc and pcqc, mr_y_qc__post__series for
  postqc, mr_y_qc__raw__man__series for the manual raw QC consensus. The
  files' entity columns (acq, acq_time, dir, modality, rec, run, task) are
  kept as columns. The release has no series-level tables.

  The manual QC per-reviewer files are routed by file name:
  desc-reviews_manualrawqc (or the older desc-manualrawqc_qc) to
  mr_y_qc__raw__man__series__review, keyed on participant_id, session_id,
  filename and reviewer_id; desc-<modality>_manualpostqc (or the older
  desc-manualpostqc<modality>_qc) to mr_y_qc__post__man__<modality>__review,
  keyed on participant_id, session_id and reviewer_id. A per-reviewer file
  without its table's column names (an old desc-manualpostqc_qc, or a
  desc-manualrawqc_qc with the short names), written before the rename and
  rewritten by the tools that write them, is skipped and reported.

  The series and per-reviewer tables are written to series_outdir, default
  outdir.

  tables selects the tables written: 'session' only the session-level ones
  (into outdir), 'series' only the series-level and per-reviewer ones (into
  series_outdir), 'all' (default) both. Every file is still read; a table
  of the other kind is neither written nor reported, nor are its conflicts
  counted.

  Each table is written with the release's column order first and new
  columns after, rows sorted by participant_id and session_id (and
  filename). Values are written as read, less the quoting and padding of
  the MATLAB-written TSVs. A column two files give different values for is
  a conflict: the first value is kept and errcode is set.

  Writes <outdir>/compile_tabdat_report.txt: the tables requested;
  series_outdir, when not outdir; per table, its rows, the release columns
  present and missing, and the new columns; the columns no rule could
  place; every conflict; and, with id_events_file, the visits under the
  rootdirs that it does not list.

  A file listed but gone when it is read (removed by a module rerun on the
  session while the compile runs) is skipped with a warning and listed in
  the report; the other read errors stop the compile.

  Returns (written, errcode, errmsg), where written maps each table name to
  its output path.

Created:  09/28/2026 by Don Hagler
Prev mod: 09/29/2026 by Don Hagler
Prev mod: 10/05/2026 by Don Hagler
Prev mod: 10/09/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/13/2026 by Don Hagler
Prev mod: 11/15/2026 by Don Hagler
Prev mod: 11/18/2026 by Don Hagler
Last mod: 11/21/2026 by Don Hagler
"""
from __future__ import print_function
import argparse
import csv
import glob
import json
import os
import re
import sys

SCRIPTNAME = os.path.basename(__file__)

TABLE_PREFIX = 'mr_y_'
SERIES_TABLE = 'mr_y_qc__raw__series'
# new 8.0 modality-level tables, routed like the release tables
NEW_TABLES = ['mr_y_qc__raw__fmri']
KEYS = ['participant_id', 'session_id']
SERIES_KEY = 'filename'
REPORT_NAME = 'compile_tabdat_report.txt'
# the tables argument: which kind of table is written
TABLE_CHOICES = ['all', 'session', 'series']
NA = 'n/a'
EMPTY = ('', NA)

def _msg(verbose, text):
    if verbose:
        print('{}: {}'.format(SCRIPTNAME, text))


def load_release_tables(phenotype_dir):
    """Return (col2tab, tab2cols) from the release's per-table json."""
    col2tab = {}
    tab2cols = {}
    for fname in sorted(glob.glob(os.path.join(phenotype_dir,
                                               TABLE_PREFIX + '*.json'))):
        table = os.path.splitext(os.path.basename(fname))[0]
        with open(fname) as fid:
            entries = json.load(fid)
        cols = [k for k in entries if k != 'MeasurementToolMetadata']
        tab2cols[table] = cols
        for col in cols:
            col2tab[col] = table
    return col2tab, tab2cols


BEHAV_SUFFIX = '_behavioral.tsv'
# not tabulated data: module status and ROI-by-ROI matrices
SKIP_SUFFIXES = ('_status.tsv', '_corrmat.tsv')

# the manual QC per-reviewer files, by file name, then by the older names the
#   review tools convert on their next save: (pattern, table, keys); the
#   table's '{}' takes the pattern's group, the modality (desc-session is the
#   consensus, not a per-reviewer file)
REVIEW_KEY = 'reviewer_id'
REVIEW_FILES = [
    (r'_desc-reviews_manualrawqc\.tsv$', 'mr_y_qc__raw__man__series__review',
     KEYS + [SERIES_KEY, REVIEW_KEY]),
    (r'_desc-(?!session_)([a-z0-9]+)_manualpostqc\.tsv$',
     'mr_y_qc__post__man__{}__review', KEYS + [REVIEW_KEY]),
    (r'_desc-manualrawqc_qc\.tsv$', 'mr_y_qc__raw__man__series__review',
     KEYS + [SERIES_KEY, REVIEW_KEY]),
    (r'_desc-manualpostqc([a-z0-9]+)_qc\.tsv$',
     'mr_y_qc__post__man__{}__review', KEYS + [REVIEW_KEY]),
    (r'_desc-manualpostqc_qc\.tsv$', None, None),
]


def read_id_events(fname):
    """The (participant_id, session_id) of every visit in an id_events
    release csv, which lists every kind of visit, imaging or not."""
    visits = set()
    with open(fname) as fid:
        reader = csv.DictReader(fid)
        if not reader.fieldnames or not set(KEYS) <= set(reader.fieldnames):
            raise ValueError('{} has no {} columns'.format(
                fname, ' and '.join(KEYS)))
        for row in reader:
            visits.add((row['participant_id'].strip(),
                        row['session_id'].strip()))
    return visits


def find_session_files(rootdir, sublist=None, seslist=None, visits=None,
                       left_out=None):
    """Every tabulated TSV of every session under rootdir.

    With visits, a set of (participant_id, session_id), a session not in it
    is skipped and appended to left_out as (participant_id, session_id,
    rootdir).

    A task's _behavioral.tsv may be in func/ (runs with a scan) and beh/
    (without); the func/ one is the session's, the beh/ one only when there
    is no func/ one, so the two never meet in one row.
    """
    files = []
    for sesdir in sorted(glob.glob(os.path.join(rootdir, 'sub-*', 'ses-*'))):
        pid = os.path.basename(os.path.dirname(sesdir))
        sid = os.path.basename(sesdir)
        if sublist and pid not in sublist:
            continue
        if seslist and sid not in seslist:
            continue
        if visits is not None and (pid, sid) not in visits:
            if left_out is not None:
                left_out.append((pid, sid, rootdir))
            continue
        fnames = sorted(glob.glob(os.path.join(sesdir, '*.tsv')) +
                        glob.glob(os.path.join(sesdir, '*', '*.tsv')))
        func_behav = set(os.path.basename(f) for f in fnames
                         if f.endswith(BEHAV_SUFFIX) and
                         os.path.basename(os.path.dirname(f)) == 'func')
        for fname in fnames:
            base = os.path.basename(fname)
            if base.endswith(SKIP_SUFFIXES):
                continue
            if (base.endswith(BEHAV_SUFFIX) and base in func_behav and
                    os.path.basename(os.path.dirname(fname)) == 'beh'):
                continue
            files.append((pid, sid, fname))
    return files


def read_tsv(fname):
    """Header and rows of a TSV, less MATLAB's quoting and padding."""
    with open(fname) as fid:
        rows = list(csv.reader(fid, delimiter='\t'))
    if not rows:
        return [], []
    header = [h.strip().strip('"') for h in rows[0]]
    body = [[v.strip().strip('"') for v in row] for row in rows[1:] if row]
    return header, body


# a component left out of the table name: tfMRI var and gpnet corr of the
#   residual and task-predicted time series share the plain tables
TABLE_DROP_COMPONENTS = [
    re.compile(r'^(mr_y_tfmri__[a-z]+__var)__(?:resid|task)__'),
    re.compile(r'^(mr_y_tfmri__[a-z]+__corr__gpnet(?:__aseg)?)__(?:resid|task)__'),
]


def table_col(col):
    """The name a column's table is chosen by: col, less a component the
    table name leaves out (TABLE_DROP_COMPONENTS)."""
    for pattern in TABLE_DROP_COMPONENTS:
        m = pattern.match(col)
        if m:
            return m.group(1) + '__' + col[m.end():]
    return col


def shared_prefix_table(cols):
    """A table name from the shared prefix of a file's columns, or None.

    Each column's statistic suffix, a last part after a single '_' (_mean,
    _sum), comes off first, so that a table-level column such as
    <table>_mean names the table itself rather than cutting it short. If the
    shared prefix is then one of those stems, and others extend it, it is the
    table; otherwise it is cut back to its last '__'. It must keep at least two components after
    mr_y_, e.g. mr_y_tfmri__mid__var__gp.
    """
    cols = [table_col(c) for c in cols if c.startswith(TABLE_PREFIX)]
    if not cols:
        return None
    stems = [re.sub(r'(?<=[^_])_[A-Za-z0-9]+$', '', c) for c in cols]
    prefix = os.path.commonprefix(stems)
    if '__' not in prefix:
        return None
    if prefix in stems and any(t.startswith(prefix + '__') for t in stems):
        table = prefix
    elif prefix.endswith('__'):
        table = prefix[:-2]
    else:
        table = prefix[:prefix.rfind('__')]
    if table.count('__') < 1:
        return None
    return table


def series_table(cols):
    """The series table of a series file, from its column names.

    The first mr_y_ column up to its first __series (followed by '_' or the
    end of the name) names the table, e.g. mr_y_qc__raw__series (rawqc, pc,
    pcqc), mr_y_qc__post__series (postqc) or mr_y_qc__raw__man__series
    (manual raw QC). A column with no __series gives its first two parts, the
    family, and the table <family>__series, as for the old
    mr_y_qc__post__aut__<modality>__ names of postqc.
    """
    for col in cols:
        if not col.startswith(TABLE_PREFIX):
            continue
        m = re.match(r'^(.+?__series)(?:_|$)', col)
        if m:
            return m.group(1)
        if col.count('__') >= 2:
            return '__'.join(col.split('__')[:2]) + '__series'
    return SERIES_TABLE


def review_table(fname, cols):
    """(table, keys, reason) of a manual QC per-reviewer file, by its name:
    (None, None, None) for any other file; (None, None, reason) for one to
    skip, without its table's column names."""
    base = os.path.basename(fname)
    for pattern, table, keys in REVIEW_FILES:
        m = re.search(pattern, base)
        if not m:
            continue
        if table is None:
            return None, None, 'old per-reviewer file'
        table = table.format(*m.groups())
        if not set(keys) <= set(cols) or \
                not any(c.startswith(table + '_') for c in cols):
            return None, None, 'no {} columns (old names)'.format(table)
        return table, keys, None
    return None, None, None


def longest_table(col, tables):
    """The longest of tables that col begins with, then '_', or None.

    tables is sorted longest first. '_' covers both a component, as in
    <table>__<roi>__lh_mean, and a table-level suffix, as in <table>_mean.
    """
    for table in tables:
        if col.startswith(table + '_'):
            return table
    return None


class Router(object):
    """Assigns session-level columns to tables.

    A column the release json lists goes to its table at once. The rest wait
    until every file is read (defer), then go to the longest release or new
    table they begin with (place).
    """

    def __init__(self, col2tab, tab2cols):
        self.col2tab = col2tab
        self.tab2cols = tab2cols
        self.release = set(tab2cols) | set(NEW_TABLES)
        self.new = set()

    def split(self, cols):
        """{col: table} for the listed columns, and the unlisted ones."""
        placed = {}
        rest = []
        for col in cols:
            table = self.col2tab.get(col)
            if table:
                placed[col] = table
            else:
                rest.append(col)
        return placed, rest

    def note_file(self, rest):
        """Make a new table from the shared prefix of a file's unlisted
        columns that begin with no release table name."""
        tables = sorted(self.release, key=len, reverse=True)
        orphans = [c for c in rest
                   if longest_table(table_col(c), tables) is None]
        if orphans:
            table = shared_prefix_table(orphans)
            if table:
                self.new.add(table)

    def vaguer(self, table, cols):
        """True when table is vaguer than the release table one of cols
        begins with, as the shared prefix of a file mixing families is, e.g.
        mr_y_qc__raw for mr_y_qc__raw__dmri__fm__* and __fmri__fm__*."""
        tables = sorted(self.release, key=len, reverse=True)
        for col in cols:
            rel = longest_table(table_col(col), tables)
            if rel is not None and not (table == rel or
                                        table.startswith(rel + '__')):
                return True
        return False

    def blocked(self, table, table_cols):
        """True when a release or new table that table begins with already
        holds columns of table's family, e.g. mr_y_qc__raw__dmri holding
        mr_y_qc__raw__dmri__fm__*; table_cols maps each new table to the
        columns placed in it."""
        for other in self.release | self.new:
            if other == table or not table.startswith(other + '__'):
                continue
            cols = self.tab2cols.get(other, []) + table_cols.get(other, [])
            if any(table_col(c).startswith(table + '_') for c in cols):
                return True
        return False

    def place(self, col):
        tables = sorted(self.release | self.new, key=len, reverse=True)
        return longest_table(table_col(col), tables)


class Tables(object):
    """Rows of each table, keyed on the table's key columns."""

    def __init__(self):
        self.rows = {}        # table -> {key tuple: {col: value}}
        self.cols = {}        # table -> [col, ...] in first-seen order
        self.conflicts = []   # (table, key, col, kept, other, fname)

    def add(self, table, key, values, fname):
        trows = self.rows.setdefault(table, {})
        tcols = self.cols.setdefault(table, [])
        row = trows.setdefault(key, {})
        for col, val in values.items():
            if col not in row:
                if col not in tcols:
                    tcols.append(col)
                row[col] = val
                continue
            old = row[col]
            if old in EMPTY:
                row[col] = val
            elif val not in EMPTY and val != old:
                self.conflicts.append((table, key, col, old, val, fname))


def compile_bids_tabdat(rootdirs, outdir, phenotype_dir,
                        sublist=None, seslist=None, id_events_file=None,
                        series_outdir=None, tables='all', forceflag=False,
                        verbose=True):
    """Compile per-session TSVs into one TSV per table.

    Parameters
    ----------
    rootdirs : list of str
        BIDS derivatives trees, each holding sub-*/ses-*/, e.g. a
        cmig_analysis, cmig_postqc or cmig_rawqc tree.
    outdir : str
        Output directory, created if absent.
    phenotype_dir : str
        Directory of the release's per-table json.
    sublist, seslist : list of str, optional
        Participants and sessions to include; default all.
    id_events_file : str, optional
        id_events release csv; only the visits it lists (participant_id,
        session_id) are compiled; default all.
    series_outdir : str, optional
        Output directory of the series and per-reviewer tables, created if
        absent; default outdir.
    tables : str
        'all' (default), 'session' (only the session-level tables) or
        'series' (only the series-level and per-reviewer tables).
    forceflag : bool
        Overwrite an existing <outdir>/<table>.tsv (or
        <series_outdir>/<table>.tsv).
    verbose : bool
        Print progress.

    Returns
    -------
    written : dict {table: path}
    errcode : int    0=success 1=error
    errmsg : str
    """
    written = {}
    if tables not in TABLE_CHOICES:
        return written, 1, 'tables must be one of {}: {}'.format(
            ', '.join(TABLE_CHOICES), tables)
    # the kind of table requested (tables below is the compiled Tables)
    kind = tables
    if isinstance(rootdirs, str):
        rootdirs = [rootdirs]
    for rootdir in rootdirs:
        if not os.path.isdir(rootdir):
            return written, 1, 'rootdir not found: {}'.format(rootdir)
    if not os.path.isdir(phenotype_dir):
        return written, 1, 'phenotype_dir not found: {}'.format(phenotype_dir)
    col2tab, tab2cols = load_release_tables(phenotype_dir)
    if not tab2cols:
        return written, 1, 'no {}*.json in {}'.format(TABLE_PREFIX,
                                                       phenotype_dir)
    visits = None
    left_out = []      # (pid, sid, rootdir) of the visits not in id_events_file
    if id_events_file:
        if not os.path.isfile(id_events_file):
            return written, 1, 'id_events_file not found: {}'.format(
                id_events_file)
        try:
            visits = read_id_events(id_events_file)
        except ValueError as e:
            return written, 1, str(e)
        _msg(verbose, '{} visits in {}'.format(len(visits), id_events_file))
    router = Router(col2tab, tab2cols)
    tables = Tables()
    deferred = []      # (key, {col: value}, fname, candidate table)
    candidates = set() # tables of their own, files no json lists
    unplaced = {}      # col -> first file
    skipped = []       # (fname, reason)
    vanished = []      # files listed but gone when read
    series_tables = {} # series and per-reviewer tables -> key columns
    nfiles = 0

    for rootdir in rootdirs:
        files = find_session_files(rootdir, sublist, seslist, visits,
                                   left_out)
        _msg(verbose, '{} tsv files under {}'.format(len(files), rootdir))
        for pid, sid, fname in files:
            try:
                header, body = read_tsv(fname)
            except FileNotFoundError:
                # removed after it was listed, e.g. a raw QC rerun's
                #   errors.tsv, which the module removes when all scans pass
                vanished.append(fname)
                print('{}: WARNING: {} no longer exists, skipped'.format(
                    SCRIPTNAME, fname))
                continue
            if not header or not body:
                skipped.append((fname, 'empty'))
                continue
            # a _behavioral.tsv carries no ids; they come from its path
            if fname.endswith(BEHAV_SUFFIX) and not set(KEYS) & set(header):
                header = KEYS + header
                body = [[pid, sid] + r for r in body]
            if header[:2] != KEYS:
                skipped.append((fname, 'no participant_id, session_id'))
                continue
            rtable, rkeys, reason = review_table(fname, header)
            if reason:
                skipped.append((fname, reason))
                continue
            nfiles += 1
            rows = [dict(zip(header, r)) for r in body]
            for row in rows:
                if row.get('participant_id') != pid or \
                        row.get('session_id') != sid:
                    return written, 1, '{} has {} {}, its path {} {}'.format(
                        fname, row.get('participant_id'),
                        row.get('session_id'), pid, sid)
            cols = [c for c in header if c not in KEYS]
            if rtable:
                series_tables[rtable] = rkeys
                for row in rows:
                    values = dict((c, row[c]) for c in header
                                  if c not in rkeys)
                    tables.add(rtable, tuple(row[k] for k in rkeys), values,
                               fname)
                continue
            if SERIES_KEY in header:
                stable = series_table(cols)
                series_tables[stable] = KEYS + [SERIES_KEY]
                for row in rows:
                    values = dict((c, row[c]) for c in cols
                                  if c != SERIES_KEY)
                    tables.add(stable, (pid, sid, row[SERIES_KEY]), values,
                               fname)
                continue
            if len(rows) > 1:
                skipped.append((fname, '{} rows and no {} column'.format(
                    len(rows), SERIES_KEY)))
                continue
            placed, rest = router.split(cols)
            # a file the release does not list at all may be a table of its
            #   own, named for its shared prefix, e.g.
            #   mr_y_rsfmri__corr__gpnet__trim5 rather than the release's
            #   mr_y_rsfmri__corr__gpnet; decided once every file is read
            if not placed and rest:
                table = shared_prefix_table(rest)
                if table and all(table_col(c).startswith(table + '_')
                                 for c in rest) \
                        and not router.vaguer(table, rest):
                    candidates.add(table)
                    deferred.append(((pid, sid),
                                     dict((c, rows[0][c]) for c in rest),
                                     fname, table))
                    continue
            by_table = {}
            for col, table in placed.items():
                by_table.setdefault(table, {})[col] = rows[0][col]
            for table, values in by_table.items():
                tables.add(table, (pid, sid), values, fname)
            if rest:
                router.note_file(rest)
                deferred.append(((pid, sid),
                                 dict((c, rows[0][c]) for c in rest), fname,
                                 None))

    # a candidate is a table of its own unless a table it begins with already
    #   holds columns of its family, placed from the other files
    table_cols = {}
    for key, values, fname, cand in deferred:
        if cand is None:
            for col in values:
                t = router.place(col)
                if t is not None:
                    table_cols.setdefault(t, []).append(col)
    own = set(t for t in candidates if not router.blocked(t, table_cols))
    router.new |= own

    # the unlisted columns, now that every new table is known
    placement = {}
    for key, values, fname, cand in deferred:
        by_table = {}
        for col, val in values.items():
            if col not in placement:
                placement[col] = cand if cand in own else router.place(col)
            table = placement[col]
            if table is None:
                unplaced.setdefault(col, fname)
                continue
            by_table.setdefault(table, {})[col] = val
        for table, values in by_table.items():
            tables.add(table, key, values, fname)

    if visits is not None:
        _msg(verbose, '{} visits not in {}, left out'.format(
            len(left_out), id_events_file))
    _msg(verbose, '{} files read, {} tables'.format(nfiles, len(tables.rows)))
    if vanished:
        print('{}: WARNING: {} files gone when read, skipped'.format(
            SCRIPTNAME, len(vanished)))
    if not series_outdir:
        series_outdir = outdir
    # the tables requested, of the kind tables names
    wanted = sorted(t for t in tables.rows
                    if kind in ('all', 'series' if t in series_tables
                                else 'session'))
    os.makedirs(outdir, exist_ok=True)
    if any(t in series_tables for t in wanted):
        os.makedirs(series_outdir, exist_ok=True)
    report = []
    for table in wanted:
        tdir = series_outdir if table in series_tables else outdir
        path = os.path.join(tdir, table + '.tsv')
        keys = series_tables.get(table, KEYS)
        release = tab2cols.get(table, [])
        present = set(tables.cols[table])
        cols = [c for c in release if c in present]
        new = sorted(c for c in present if c not in set(release))
        missing = [c for c in release if c not in present]
        report.append((table, len(tables.rows[table]), len(cols), len(release),
                       missing, new))
        if os.path.exists(path) and not forceflag:
            _msg(verbose, '{} exists, skipping (forceflag=False)'.format(path))
            continue
        with open(path, 'w') as fid:
            writer = csv.writer(fid, delimiter='\t', lineterminator='\n')
            writer.writerow(keys + cols + new)
            for key in sorted(tables.rows[table]):
                row = tables.rows[table][key]
                # n/a for a column the row lacks or an empty value, since a
                #   BIDS tsv cell is never empty
                writer.writerow(list(key) +
                                [row.get(c) or NA for c in cols + new])
        written[table] = path

    # the conflicts and unplaced (session-level) columns of the tables
    #   requested
    wanted_set = set(wanted)
    conflicts = [c for c in tables.conflicts if c[0] in wanted_set]
    if kind == 'series':
        unplaced = {}
    write_report(os.path.join(outdir, REPORT_NAME), rootdirs, nfiles, report,
                 unplaced, skipped, conflicts, id_events_file,
                 left_out, series_outdir if series_outdir != outdir else None,
                 kind, vanished)
    _msg(verbose, 'wrote {} tables and {}'.format(len(written), REPORT_NAME))
    errcode = 0
    errmsg = ''
    if conflicts:
        errcode = 1
        errmsg = '{} conflicting values; see {}'.format(
            len(conflicts), REPORT_NAME)
    return written, errcode, errmsg


def write_report(path, rootdirs, nfiles, report, unplaced, skipped,
                 conflicts, id_events_file=None, left_out=None,
                 series_outdir=None, kind='all', vanished=None):
    with open(path, 'w') as fid:
        fid.write('{} report\n'.format(SCRIPTNAME))
        for rootdir in rootdirs:
            fid.write('  root {}\n'.format(rootdir))
        if id_events_file:
            fid.write('  id_events {}\n'.format(id_events_file))
        fid.write('  tables {}\n'.format(
            {'all': 'all (session, series and per-reviewer)',
             'session': 'session (session-level only)',
             'series': 'series (series-level and per-reviewer only)'}[kind]))
        if series_outdir and kind != 'session':
            fid.write('  series tables in {}\n'.format(series_outdir))
        fid.write('  {} files, {} tables\n'.format(nfiles, len(report)))
        fid.write('  {} files gone when read\n\n'.format(len(vanished or [])))
        fid.write('{:<60s} {:>6s} {:>9s} {:>8s} {:>5s}\n'.format(
            'table', 'rows', 'release', 'missing', 'new'))
        for table, nrows, npresent, nrelease, missing, new in report:
            fid.write('{:<60s} {:>6d} {:>4d}/{:<4d} {:>8d} {:>5d}\n'.format(
                table, nrows, npresent, nrelease, len(missing), len(new)))
        fid.write('\nMISSING RELEASE COLUMNS, per table\n')
        for table, _, _, _, missing, _ in report:
            for col in missing:
                fid.write('  {}\t{}\n'.format(table, col))
        fid.write('\nNEW COLUMNS, per table\n')
        for table, _, _, _, _, new in report:
            for col in new:
                fid.write('  {}\t{}\n'.format(table, col))
        fid.write('\nCOLUMNS NO RULE PLACED ({})\n'.format(len(unplaced)))
        for col in sorted(unplaced):
            fid.write('  {}\t{}\n'.format(col, os.path.basename(unplaced[col])))
        fid.write('\nFILES SKIPPED ({})\n'.format(len(skipped)))
        for fname, reason in skipped:
            fid.write('  {}\t{}\n'.format(fname, reason))
        fid.write('\nFILES GONE WHEN READ, SKIPPED ({})\n'.format(
            len(vanished or [])))
        for fname in vanished or []:
            fid.write('  {}\n'.format(fname))
        fid.write('\nCONFLICTING VALUES ({})\n'.format(len(conflicts)))
        for table, key, col, old, val, fname in conflicts:
            fid.write('  {}\t{}\t{}\tkept {}\tother {}\t{}\n'.format(
                table, ' '.join(key), col, old, val, os.path.basename(fname)))
        if id_events_file:
            fid.write('\nVISITS NOT IN THE ID_EVENTS FILE, LEFT OUT ({})\n'.format(
                len(left_out)))
            for pid, sid, rootdir in sorted(left_out):
                fid.write('  {}\t{}\t{}\n'.format(pid, sid, rootdir))


def build_parser():
    parser = argparse.ArgumentParser(
        description='Compile the per-session TSVs under BIDS derivatives '
                    'trees into one TSV per table, in the release layout.')
    parser.add_argument('rootdirs', nargs='+',
                        help='BIDS derivatives trees holding sub-*/ses-*/')
    parser.add_argument('--outdir', required=True,
                        help='Output directory')
    parser.add_argument('--phenotype_dir', required=True,
                        help='Release per-table json directory '
                             '(rawdata/phenotype of the release)')
    parser.add_argument('--sublist', nargs='+', default=None,
                        help='Participants to include (default: all)')
    parser.add_argument('--seslist', nargs='+', default=None,
                        help='Sessions to include (default: all)')
    parser.add_argument('--id_events_file', default=None,
                        help='id_events release csv; only the visits it lists '
                             '(participant_id, session_id) are compiled '
                             '(default: all)')
    parser.add_argument('--series_outdir', default=None,
                        help='Output directory of the series-level tables, '
                             '<...>__series, and the manual QC per-reviewer '
                             'tables, <...>__review (default: --outdir)')
    parser.add_argument('--tables', default='all', choices=TABLE_CHOICES,
                        help='Tables to write: session, the session-level '
                             'tables (into --outdir); series, the '
                             'series-level and per-reviewer tables (into '
                             '--series_outdir); all, both (default: all)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing tables (default: 0)')
    parser.add_argument('--verbose', type=int, default=1, choices=[0, 1],
                        help='Display status messages (default: 1)')
    return parser


def parse_args():
    return build_parser().parse_args()


def run(args):
    written, errcode, errmsg = compile_bids_tabdat(
        args.rootdirs, args.outdir, phenotype_dir=args.phenotype_dir,
        sublist=args.sublist, seslist=args.seslist,
        id_events_file=args.id_events_file,
        series_outdir=args.series_outdir, tables=args.tables,
        forceflag=bool(args.forceflag),
        verbose=bool(args.verbose))
    if errcode:
        print('{}: ERROR: {}'.format(SCRIPTNAME, errmsg), file=sys.stderr)
    return errcode, errmsg


def main():
    errcode, _ = run(parse_args())
    sys.exit(errcode)


if __name__ == '__main__':
    main()
