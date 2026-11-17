#!/usr/bin/env python3
"""
Generate cluster batch job scripts for concatenating cmig_bids per-session
concat-format mat files across subject-sessions.

One job per target set, where a target set is a datatype and, with expand_flag,
one measure x smoothing combination within it. The concatenation function takes
arbitrary measure, smoothing and layer lists; this script decides what any
given run asks for.

With several analysis root dirs, their sessions are combined into one set of
files; with an id_events release csv, only the visits it lists are kept.

Created:  08/28/2026 by Don Hagler
Prev mod: 08/28/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/14/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 11/04/2026 by Don Hagler
Prev mod: 11/13/2026 by Don Hagler
Prev mod: 11/14/2026 by Don Hagler
Last mod: 11/15/2026 by Don Hagler
"""

import argparse
import os
import sys

import cmig_write_scripts

APPNAME = 'cmig_concat_bids_sessions'
# the compiled app, and its image, that run every module, given its name first
DISPATCH_APPNAME = 'cmig_run_bids_module'

DEFAULT_BATCHDIR = os.path.join(os.path.expanduser('~'), 'batchdirs')
DEFAULT_APPVER = 'mmps_dev'
DEFAULT_APPDIR = '/usr/pubsw/packages/MMPS/apps'
DEFAULT_PUBSDIR = '/usr/pubsw'
DEFAULT_MRTDIR = '/usr/pubsw/packages/matlab/R2023a'

# datatypes, and which lists vary within each. Only the datatypes whose measures
# are enumerated per measure x smoothing can be expanded into separate jobs.
DATATYPE_CONFIG = {
    'smri_measures': {'measlist_arg': 'smri_measlist', 'expandable': True},
    'smri_surfaces': {'measlist_arg': None, 'expandable': False},
    'dmri_dti': {'measlist_arg': 'DTI_measlist', 'expandable': True},
    'dmri_isdti': {'measlist_arg': 'isDTI_measlist', 'expandable': True},
    'dmri_rsi': {'measlist_arg': 'RSI_measlist', 'expandable': True},
    'tfmri_stats': {'measlist_arg': None, 'expandable': False},
    'fmri_corr': {'measlist_arg': None, 'expandable': False},
    'fmri_tseries': {'measlist_arg': None, 'expandable': False},
}

# the 7.0 released set; the concatenation function imposes no scope of its own
DEFAULT_SMRI_MEASLIST = ['thk', 'area', 'sulc', 't1_gwc', 't2_gwc',
                         't1_gm', 't2_gm', 't1_wm', 't2_wm']
DEFAULT_DTI_MEASLIST = ['FA', 'MD', 'LD', 'TD']
DEFAULT_ISDTI_MEASLIST = ['FA', 'MD', 'LD', 'TD']
DEFAULT_RSI_MEASLIST = ['RNI', 'RND', 'RNT', 'HNT', 'FNI']
DEFAULT_SMOOTHLIST = [0, 16, 256, 1000]
DEFAULT_LAYERLIST = ['gm', 'wm']
DEFAULT_SURFLIST = ['white', 'pial']
DEFAULT_HEMILIST = ['lh', 'rh']
DEFAULT_TASKLIST = ['mid', 'nback', 'sst', 'rest']
DEFAULT_ATLASLIST = ['gpaseg']


def _write_vol_info(datatype, analysis_rootdir, output_rootdir, seslist,
                    id_events_file, ico_order, forceflag, compiled_flag, appver,
                    appdir, pubsdir, mrtdir):
    """Run cmig_concat_bids_sessions with vol_info_flag 1, before any job is
    written; return False if the datatype has no data (its invocation status
    has errcode 0 and errmsg "no <datatype> data in ...")."""
    import cmig_run_bids_module
    print('writing vol_info for {}...'.format(datatype))
    errcode, errmsg = cmig_run_bids_module.run_module(
        'cmig_concat_bids_sessions',
        datatype=datatype,
        analysis_rootdir=analysis_rootdir,
        output_rootdir=output_rootdir,
        seslist=list(seslist) if seslist else None,
        id_events_file=id_events_file or None,
        ico_order=ico_order,
        vol_info_flag=1,
        forceflag=int(bool(forceflag)),
        compiled_flag=compiled_flag,
        appver=appver, appdir=appdir, pubsdir=pubsdir, mrtdir=mrtdir)
    if errcode:
        sys.exit('Error: failed to write vol_info for {}: {}'.format(
            datatype, errmsg))
    # the invocation status carries an errmsg with errcode 0 only for a
    #   no-data skip ("no <datatype> data in ..."), vol_info read back from an
    #   earlier run included
    return not errmsg


def _build_job_list(datatype, measlist, smoothlist, expand_flag):
    """Return a list of (job_name, per-job parms) tuples."""
    cfg = DATATYPE_CONFIG[datatype]
    if not (expand_flag and cfg['expandable']):
        parms = {}
        if cfg['measlist_arg'] and measlist:
            parms['measlist'] = measlist
        parms['smoothlist'] = smoothlist
        return [('concat_{}'.format(datatype), parms)]
    jobs = []
    for meas in measlist:
        for smooth in smoothlist:
            name = 'concat_{}_{}_sm{}'.format(datatype, meas, smooth)
            jobs.append((name, {'measlist': [meas], 'smoothlist': [smooth]}))
    return jobs


def cmig_mkjobs_bids_concat(
    analysis_rootdir,
    output_rootdir,
    datatypelist=None,
    batchdir=DEFAULT_BATCHDIR,
    study='',
    seslist=None,
    id_events_file=None,
    smri_measlist=None,
    DTI_measlist=None,
    isDTI_measlist=None,
    RSI_measlist=None,
    smoothlist=None,
    layerlist=None,
    surflist=None,
    hemilist=None,
    tasklist=None,
    atlaslist=None,
    ico_order=5,
    expand_flag=False,
    forceflag=False,
    compiled_flag=1,
    appver=DEFAULT_APPVER,
    appdir=DEFAULT_APPDIR,
    pubsdir=DEFAULT_PUBSDIR,
    mrtdir=DEFAULT_MRTDIR,
):
    """Create a batch job directory of concatenation scripts.

    Parameters
    ----------
    analysis_rootdir : str or list of str
        root of the per-session cmig_analysis tree, or several roots whose
        sessions are combined
    output_rootdir : str
        root directory for the concatenated output
    datatypelist : list of str, optional
        datatypes to concatenate; default is all of them
    id_events_file : str, optional
        id_events release csv; only the visits it lists (participant_id,
        session_id) are kept; default all
    expand_flag : bool
        if True, one job per measure x smoothing for the expandable datatypes

    Returns
    -------
    jobdir : str
        the batch job directory that was written
    """
    # one root is passed as a string, as before; several as a cell array
    if isinstance(analysis_rootdir, (list, tuple)):
        analysis_rootdir = list(analysis_rootdir)
        if len(analysis_rootdir) == 1:
            analysis_rootdir = analysis_rootdir[0]
    if datatypelist is None:
        datatypelist = sorted(DATATYPE_CONFIG)
    for datatype in datatypelist:
        if datatype not in DATATYPE_CONFIG:
            raise ValueError('unknown datatype: {}'.format(datatype))

    if smri_measlist is None:
        smri_measlist = list(DEFAULT_SMRI_MEASLIST)
    if DTI_measlist is None:
        DTI_measlist = list(DEFAULT_DTI_MEASLIST)
    if isDTI_measlist is None:
        isDTI_measlist = list(DEFAULT_ISDTI_MEASLIST)
    if RSI_measlist is None:
        RSI_measlist = list(DEFAULT_RSI_MEASLIST)
    if smoothlist is None:
        smoothlist = list(DEFAULT_SMOOTHLIST)
    if layerlist is None:
        layerlist = list(DEFAULT_LAYERLIST)
    if surflist is None:
        surflist = list(DEFAULT_SURFLIST)
    if hemilist is None:
        hemilist = list(DEFAULT_HEMILIST)
    if tasklist is None:
        tasklist = list(DEFAULT_TASKLIST)
    if atlaslist is None:
        atlaslist = list(DEFAULT_ATLASLIST)

    jobdir_name = '{}_concat'.format(study) if study else 'concat'
    jobdir = os.path.join(batchdir, jobdir_name)
    if not os.path.isdir(jobdir):
        os.makedirs(jobdir)

    measlist_by_datatype = {
        'smri_measures': smri_measlist,
        'dmri_dti': DTI_measlist,
        'dmri_isdti': isDTI_measlist,
        'dmri_rsi': RSI_measlist,
    }

    scriptlist = []
    for datatype in datatypelist:
        # every job of a datatype writes into the same output folders (the
        #   release 7.0 layout; fmri_corr and fmri_tseries share
        #   corrmat/rsfmri and corrmat/tfmri), so the session list is fixed
        #   once here, in each folder's vol_info, which the jobs, and the
        #   datatype sharing the folder, then read back
        if not _write_vol_info(datatype, analysis_rootdir, output_rootdir,
                               seslist, id_events_file, ico_order, forceflag,
                               compiled_flag, appver, appdir, pubsdir, mrtdir):
            print('WARNING: no {} data, no jobs written'.format(datatype))
            continue
        measlist = measlist_by_datatype.get(datatype)
        jobs = _build_job_list(datatype, measlist, smoothlist, expand_flag)
        for job_name, job_parms in jobs:
            parms = {
                'seslist': seslist if seslist else [],
                'id_events_file': id_events_file,
                'layerlist': layerlist,
                'surflist': surflist,
                'hemilist': hemilist,
                'tasklist': tasklist,
                'atlaslist': atlaslist,
                'ico_order': ico_order,
                'forceflag': int(bool(forceflag)),
            }
            parms.update(job_parms)
            # an empty value is written as a bare space, which the shell drops,
            #   leaving the tag paired with the next tag's name; the matlab
            #   defaults for these are empty anyway, so drop them here
            parms = {k: v for k, v in parms.items()
                     if not (v is None or (isinstance(v, (list, tuple)) and not v))}
            # the session list comes from the vol_info written above
            tags = sorted(parms)
            args = [datatype, analysis_rootdir, output_rootdir]
            if compiled_flag == 0:
                fname = os.path.join(jobdir, job_name + '.m')
                cmig_write_scripts.cmig_write_matlab_script(
                    fname, APPNAME, args, tags=tags, parms=parms)
            else:
                fname = os.path.join(jobdir, job_name + '.csh')
                cmig_write_scripts.cmig_write_runtime_script(
                    fname, DISPATCH_APPNAME, [APPNAME] + args, mrtdir, appdir,
                    tags=tags, parms=parms, appver=appver)
            scriptlist.append(os.path.basename(fname))
            print('wrote {}'.format(fname))

    fname_list = os.path.join(jobdir, 'scriptlist.txt')
    with open(fname_list, 'w') as fid:
        for name in scriptlist:
            fid.write('{}\n'.format(os.path.splitext(name)[0]))
    print('{} jobs in {}'.format(len(scriptlist), jobdir))
    return jobdir


def parse_args():
    parser = argparse.ArgumentParser(
        description='create batch jobs to concatenate cmig_bids per-session '
                    'mat files across subject-sessions')
    parser.add_argument('--analysis_rootdir', required=True, nargs='+',
                        help='root of the per-session cmig_analysis tree, or '
                             'several roots whose sessions are combined')
    parser.add_argument('--output_rootdir', required=True,
                        help='root directory for the concatenated output')
    parser.add_argument('--datatypelist', nargs='+', default=None,
                        help='datatypes to concatenate (default: all): {}'.format(
                            ' '.join(sorted(DATATYPE_CONFIG))))
    parser.add_argument('--batchdir', default=DEFAULT_BATCHDIR)
    parser.add_argument('--study', default='')
    parser.add_argument('--seslist', nargs='+', default=None,
                        help='session ids to keep, e.g. ses-00A ses-02A')
    parser.add_argument('--id_events_file', default=None,
                        help='id_events release csv; only the visits it lists '
                             '(participant_id, session_id) are kept '
                             '(default: all)')
    parser.add_argument('--smri_measlist', nargs='+', default=None)
    parser.add_argument('--DTI_measlist', nargs='+', default=None)
    parser.add_argument('--isDTI_measlist', nargs='+', default=None)
    parser.add_argument('--RSI_measlist', nargs='+', default=None)
    parser.add_argument('--smoothlist', nargs='+', type=int, default=None)
    parser.add_argument('--layerlist', nargs='+', default=None)
    parser.add_argument('--surflist', nargs='+', default=None)
    parser.add_argument('--hemilist', nargs='+', default=None)
    parser.add_argument('--tasklist', nargs='+', default=None)
    parser.add_argument('--atlaslist', nargs='+', default=None)
    parser.add_argument('--ico_order', type=int, default=5)
    parser.add_argument('--expand_flag', type=int, default=0,
                        help='one job per measure x smoothing (default: 0)')
    parser.add_argument('--forceflag', type=int, default=0)
    parser.add_argument('--compiled_flag', type=int, default=1, choices=[0, 1],
                        help='0 = uncompiled matlab .m, 1 = compiled runtime csh '
                             '(default: 1)')
    parser.add_argument('--appver', default=DEFAULT_APPVER)
    parser.add_argument('--appdir', default=DEFAULT_APPDIR)
    parser.add_argument('--pubsdir', default=DEFAULT_PUBSDIR)
    parser.add_argument('--mrtdir', default=DEFAULT_MRTDIR)
    return parser.parse_args()


def main():
    args = parse_args()
    cmig_mkjobs_bids_concat(**vars(args))
    return 0


if __name__ == '__main__':
    sys.exit(main())
