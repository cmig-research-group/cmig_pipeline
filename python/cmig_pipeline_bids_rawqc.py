#!/usr/bin/env python3
"""
Pipeline: run raw QC for all modalities (anat, func, dwi, fmap).

Unlike other pipelines, each modality runs independently — a failure in one
does not prevent the others from running. An overall status file is written
at {output_rootdir}/{pid}/{sid}/{pid}_{sid}_desc-rawqc_status.tsv, with
errcode 1 where any modality, pc or pcqc failed;
cmig_mkjobs_bids_pipeline.py checks it after the per-modality status files. With --pc_flag 1,
protocol compliance is written into --pc_rootdir (default: output_rootdir).

Created:  06/11/2026 by Don Hagler
Prev mod: 06/22/2026 by Don Hagler
Prev mod: 06/25/2026 by Don Hagler
Prev mod: 07/24/2026 by Don Hagler
Prev mod: 08/12/2026 by Don Hagler
Prev mod: 10/13/2026 by Don Hagler
Prev mod: 10/19/2026 by Don Hagler
Prev mod: 10/23/2026 by Don Hagler
Prev mod: 10/26/2026 by Don Hagler
Prev mod: 10/29/2026 by Don Hagler
Prev mod: 11/06/2026 by Don Hagler
Prev mod: 11/11/2026 by Don Hagler
Prev mod: 11/12/2026 by Don Hagler
Last mod: 11/20/2026 by Don Hagler
"""

import argparse
import os
import sys
from argparse import Namespace

import cmig_run_bids_module
import cmig_sing_utils
import cmig_write_status as cmig_status

APPNAME = 'cmig_pipeline_bids_rawqc'

DEFAULT_APPVER = 'mmps_dev'
DEFAULT_APPDIR = '/usr/pubsw/packages/MMPS/apps'
DEFAULT_PUBSDIR = '/usr/pubsw'
DEFAULT_MRTDIR = '/usr/pubsw/packages/matlab/R2023a'
DEFAULT_MRTDIR_TARGET = '/opt/matlabruntime/R2023a'

MODALITIES = ['anat', 'func', 'dwi', 'fmap']


def _demarcate(label):
    print('\n' + '=' * 60)
    print(label)
    print('=' * 60)


def _write_status(outdir, pid, sid, errcode, errmsg, startdate):
    fpath = cmig_status.cmig_write_status(
        pid, sid, outdir, 'desc-rawqc', errcode, errmsg, startdate=startdate)
    cmig_sing_utils.add_status_versions(fpath, app_name=APPNAME)


def build_parser():
    parser = argparse.ArgumentParser(
        description='Run raw BIDS QC for all modalities (anat, func, dwi, fmap).'
    )
    parser.add_argument('--participant_id', required=True,
                        help='Participant ID (e.g., sub-XXXXXXXX)')
    parser.add_argument('--session_id', required=True,
                        help='Session ID (e.g., ses-00A)')
    parser.add_argument('--raw_rootdir', required=True,
                        help='Root directory of raw BIDS data (read-only input)')
    parser.add_argument('--output_rootdir', required=True,
                        help='Root directory for QC output (created if absent)')
    parser.add_argument('--cleanup_flag', type=int, default=1, choices=[0, 1],
                        help='Remove temporary files (default: 1)')
    parser.add_argument('--export_flag', type=int, default=1, choices=[0, 1],
                        help='Create summary images for manual review (default: 1)')
    parser.add_argument('--motion_radius', type=int, default=50,
                        help='Head radius (mm) for rotation-to-distance conversion, func+dwi (default: 50)')
    parser.add_argument('--pc_flag', type=int, default=0, choices=[0, 1],
                        help='Run protocol compliance check (cmig_pc_bids_abcd) and the pcqc merge after rawqc (default: 0)')
    parser.add_argument('--exclude_noncompliant_flag', type=int, default=0, choices=[0, 1],
                        help='With pc_flag, eprime__scan__run__valid_count also leaves out runs with compliant 0 in the rawdata _scans.tsv (default: 0)')
    parser.add_argument('--pc_rootdir', default=None,
                        help='Root directory of protocol compliance output '
                             '(default: the rawqc output_rootdir)')
    parser.add_argument('--manualqc_rootdir', default=None,
                        help='Root directory of manual raw QC output '
                             '(default: the rawqc output_rootdir)')
    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1],
                        help='Display status messages (default: 0)')
    parser.add_argument('--forceflag', type=int, default=0, choices=[0, 1],
                        help='Overwrite existing output (default: 0)')
    parser.add_argument('--compiled_flag', type=int, default=1, choices=[0, 1],
                        help='0 = MATLAB code, 1 = compiled app (default: 1)')
    parser.add_argument('--reproducible_flag', type=int, default=1, choices=[0, 1],
                        help=cmig_run_bids_module.REPRODUCIBLE_HELP)
    parser.add_argument('--appver', default=DEFAULT_APPVER,
                        help='Container image version tag (default: %(default)s)')
    parser.add_argument('--appdir', default=DEFAULT_APPDIR,
                        help='Directory containing the compiled app (default: %(default)s)')
    parser.add_argument('--pubsdir', default=DEFAULT_PUBSDIR,
                        help='Public software directory bound to /usr/pubsw (default: %(default)s)')
    parser.add_argument('--mrtdir', default=DEFAULT_MRTDIR,
                        help='MATLAB runtime directory on host (default: %(default)s)')
    parser.add_argument('--mrtdir_target', default=DEFAULT_MRTDIR_TARGET,
                        help='MATLAB runtime mount point inside container (default: %(default)s)')
    return parser


def parse_args():
    return build_parser().parse_args()


def run(args):
    startdate = cmig_status.module_start()
    pid = args.participant_id
    sid = args.session_id

    if not os.path.isdir(args.raw_rootdir):
        errmsg = 'raw_rootdir not found: {}'.format(args.raw_rootdir)
        print('ERROR: {}'.format(errmsg), file=sys.stderr)
        return 1, errmsg

    infra = dict(
        compiled_flag=args.compiled_flag,
        reproducible_flag=args.reproducible_flag,
        appver=args.appver,
        appdir=args.appdir,
        pubsdir=args.pubsdir,
        mrtdir=args.mrtdir,
        mrtdir_target=args.mrtdir_target,
    )

    pc_rootdir = args.pc_rootdir or args.output_rootdir
    results = {}

    # the modality-specific options, each as the earlier per-modality wrapper
    #   passed it
    extra = {
        'anat': {},
        'func': dict(motion_radius=args.motion_radius, method='fslmean'),
        'dwi': dict(motion_radius=args.motion_radius),
        'fmap': {},
    }
    for modality in MODALITIES:
        name = 'cmig_rawqc_bids_' + modality
        raw_subdir = os.path.join(args.raw_rootdir, pid, sid, modality)
        if not os.path.isdir(raw_subdir):
            if args.verbose:
                print('{}: skipping (no raw/{}/{}/{})'.format(name, pid, sid, modality))
            continue
        _demarcate(name)
        errcode, errmsg = cmig_run_bids_module.run_module(
            'cmig_rawqc_bids',
            participant_id=pid,
            session_id=sid,
            raw_rootdir=args.raw_rootdir,
            output_rootdir=args.output_rootdir,
            modality=modality,
            cleanup_flag=args.cleanup_flag,
            export_flag=args.export_flag,
            verbose=args.verbose,
            forceflag=args.forceflag,
            **dict(extra[modality], **infra)
        )
        results[modality] = (errcode, errmsg)

    # protocol compliance (optional)
    if args.pc_flag:
        import cmig_pc_bids_abcd
        _demarcate('cmig_pc_bids_abcd')
        errcode, errmsg = cmig_pc_bids_abcd.run(Namespace(
            participant_id=pid,
            session_id=sid,
            raw_rootdir=args.raw_rootdir,
            output_rootdir=pc_rootdir,
            verbose=args.verbose,
            forceflag=args.forceflag,
            exclude_noncompliant_flag=args.exclude_noncompliant_flag,
        ))
        results['pc'] = (errcode, errmsg)

    # merge compliance with QC, last, reading the pc and rawqc output this run
    # has just written. It runs with pc rather than on a flag of its own: the
    # merge needs the compliance output and nothing else produces it
    if args.pc_flag:
        import cmig_check_event_pcqc
        _demarcate('cmig_check_event_pcqc')
        errcode, errmsg = cmig_check_event_pcqc.run(Namespace(
            participant_id=pid,
            session_id=sid,
            pc_rootdir=pc_rootdir,
            rawqc_rootdir=args.output_rootdir,
            manualqc_rootdir=(args.manualqc_rootdir or args.output_rootdir),
            output_rootdir=args.output_rootdir,
            verbose=args.verbose,
            forceflag=args.forceflag,
        ))
        results['pcqc'] = (errcode, errmsg)

    # summary
    _demarcate('summary')
    errmsg_parts = []
    for label in (list(MODALITIES) +
                  (['pc', 'pcqc'] if args.pc_flag else [])):
        if label not in results:
            print('{}: skipped'.format(label))
            continue
        ec, em = results[label]
        status_str = 'ok' if ec == 0 else 'ERROR (errcode={})'.format(ec)
        print('{}: {}{}'.format(label, status_str, ': ' + em if em else ''))
        if em:
            errmsg_parts.append('{}: {}'.format(label, em))

    # fail if pc, pcqc or any rawqc modality failed (every modality has run)
    pc_failed = args.pc_flag and 'pc' in results and results['pc'][0] != 0
    pcqc_failed = (args.pc_flag and 'pcqc' in results and
                   results['pcqc'][0] != 0)
    any_rawqc_failed = any(ec != 0 for k, (ec, _) in results.items()
                           if k in MODALITIES)
    combined_errcode = 1 if (pc_failed or pcqc_failed or any_rawqc_failed) else 0

    combined_errmsg = '; '.join(errmsg_parts) if errmsg_parts else None

    # write overall status file
    outdir = os.path.join(args.output_rootdir, pid, sid)
    _write_status(outdir, pid, sid, combined_errcode, combined_errmsg, startdate)

    return combined_errcode, combined_errmsg


def main():
    errcode, errmsg = run(parse_args())
    if errmsg:
        print(errmsg, file=sys.stderr)
    sys.exit(errcode)


if __name__ == '__main__':
    main()
