#!/usr/bin/env python3
"""Write a status TSV and print the matching block, as cmig_write_status.m does.

The python modules that do their own work rather than calling a MATLAB app --
cmig_fsurf_bids.py and cmig_pc_bids_abcd.py -- write the same status row and
print the same block. Keeping one copy here leaves a two-way parity with
cmig_write_status.m instead of three implementations that drift.

The printed block is read by cmig_check_batchdir.py, which anchors on the
outdir line and counts four fixed offsets forward for procdate, errcode and
errmsg, then looks for the finished pair directly after them. So module and
startdate print ABOVE outdir, and nothing may be inserted between outdir and
errmsg.
"""

import datetime
import os

STATUS_PREFIX = 'cmig_write_status'
# an empty cell is written as n/a (BIDS); read_status maps it back to ''
NA = 'n/a'


def _now():
    now = datetime.datetime.now()
    return now.strftime('%Y-%m-%dT%H:%M:%S.') + '{0:03d}'.format(now.microsecond // 1000)


def status_fpath(outdir, participant_id, session_id, module):
    """The status file cmig_write_status.m would write, same naming."""
    return os.path.join(outdir, '{0}_{1}_{2}_status.tsv'.format(
        participant_id, session_id, module))


def cmig_write_status(participant_id, session_id, outdir, module, errcode, errmsg,
                      extra=None, startdate=None, procflag=1):
    """Write {participant_id}_{session_id}_{module}_status.tsv and print the block.

    module carries the desc- entity, as in cmig_write_status.m, since it becomes
    part of the filename that cmig_mkjobs_bids_pipeline.py globs for.

    extra is a dict whose keys become columns after errmsg. An empty errmsg or
    extra value is written as n/a; the printed block keeps the errmsg as given.
    startdate is the module's entry time, so a duration needs no differencing
    against the module before it.
    procflag 0 with the status file already present means the module skipped its
    work, so the file is left as it is: rewriting it would re-stamp procdate for
    an invocation that computed nothing and discard the version columns
    cmig_sing_utils.add_status_versions appended. The file-exists half matters --
    outputs can be present with no status file, and writing nothing there would
    leave mkjobs never seeing done.

    Returns the status file path, whether or not it was written.
    """
    fpath = status_fpath(outdir, participant_id, session_id, module)
    if not procflag and os.path.isfile(fpath):
        print('{0}: {1} skipped, leaving {2}'.format(
            STATUS_PREFIX, module, os.path.basename(fpath)))
        return fpath

    if not isinstance(errmsg, str):
        errmsg = ''
    for char in ('\r\n', '\n', '\r', '\t'):
        errmsg = errmsg.replace(char, ' ')
    procdate = _now()

    names, values = [], []
    if extra:
        for key in extra:
            val = extra[key]
            if isinstance(val, bool):
                val = str(int(val))
            elif val is None:
                val = NA
            else:
                val = str(val) or NA
            names.append(key)
            values.append(val)

    os.makedirs(outdir, exist_ok=True)
    header = ['participant_id', 'session_id', 'procdate', 'errcode', 'errmsg']
    row = [participant_id, session_id, procdate, str(errcode), errmsg or NA]
    if startdate is not None:
        header.append('startdate')
        row.append(startdate or NA)
    with open(fpath, 'w') as f:
        f.write('\t'.join(header + names) + '\n')
        f.write('\t'.join(row + values) + '\n')

    print('')
    print('{0}: module = {1}'.format(STATUS_PREFIX, module))
    if startdate is not None:
        print('{0}: startdate = {1}'.format(STATUS_PREFIX, startdate))
    print('{0}: outdir = {1}'.format(STATUS_PREFIX, outdir))
    print('{0}: procdate = {1}'.format(STATUS_PREFIX, procdate))
    print('{0}: errcode = {1}'.format(STATUS_PREFIX, errcode))
    print('{0}: errmsg = {1}'.format(STATUS_PREFIX, errmsg))
    return fpath


def module_start():
    """A startdate stamp, to be taken at module entry."""
    return _now()


def read_status(fpath):
    """Return the status row as a dict, or None if absent or unreadable.
    An n/a cell reads as '', so callers test for an empty errmsg as before."""
    if not os.path.isfile(fpath):
        return None
    try:
        with open(fpath) as f:
            lines = f.read().splitlines()
    except (IOError, OSError):
        return None
    if len(lines) < 2:
        return None
    return dict((k, '' if v == NA else v) for k, v in
                zip(lines[0].split('\t'), lines[1].split('\t')))


def status_time(fpath, field):
    """The named timestamp from a status file, as a datetime, or None."""
    status = read_status(fpath)
    if status is None:
        return None
    stamp = (status.get(field) or '').strip()
    if not stamp:
        return None
    for fmt in ('%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S'):
        try:
            return datetime.datetime.strptime(stamp, fmt)
        except ValueError:
            continue
    return None


def inputs_newer(status_path, input_paths,
                 field='startdate', input_field='procdate', ok_only=False):
    """The first input that ran after the output that read it, as (input
    path, its time, the output's time), or None.

    Each input is another module's output, so its own status file records when
    it ran. procdate is read rather than the data file's mtime, which a copy
    or a restore of the tree moves. The output side takes startdate, so an
    input rewritten while the run was in progress counts. No tolerance: a
    spurious hit costs one rerun, where a window hides a real one. ok_only
    skips inputs whose errcode is not 0: a failed rerun of an input left
    nothing new for the output to read.
    """
    started = status_time(status_path, field)
    if started is None:
        return None
    for fpath in input_paths:
        if ok_only:
            status = read_status(fpath)
            if status is None or (status.get('errcode') or '').strip() != '0':
                continue
        ran = status_time(fpath, input_field)
        if ran is not None and ran > started:
            return fpath, ran, started
    return None
