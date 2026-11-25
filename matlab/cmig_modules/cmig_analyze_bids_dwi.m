function [outdir,errcode,errmsg] = cmig_analyze_bids_dwi(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_analyze_bids_dwi(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Analyze diffusion MRI data from BIDS processed derivatives
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%   fsurf_rootdir: BIDS freesurfer derivatives root directory
%   output_rootdir: BIDS output root directory
%
% Optional Input ('key',value):
%   'dti_flag': [0|1] whether to run DTI tensor analysis
%     {default = 1}
%   'rsi_flag': [0|1] whether to run RSI analysis
%     {default = 1}
%   'fiber_flag': [0|1] whether to run AtlasTrack fiber analysis
%     set to 0 automatically if atlas fiber dir not found
%     {default = 1}
%   'aseg_flag': [0|1] whether to extract aseg ROI measures
%     {default = 1}
%   'aseg_disp_flag': [0|1] whether to compute MD-dispersion-weighted
%     averages for aseg ROIs (rather than plain means); the inner-shell
%     DTI MD volume supplies the weighting values for all modalities
%     (dti, dti_inner, rsi); disabled automatically if that MD is missing
%     {default = 1}
%   'aseg_dispfact': multiplicative factor applied to aseg dispersion values
%     {default = 4.7}
%   'aseg_disp_suffix': suffix attached to dispersion-weighted output files
%     {default = 'dwtd'}
%   'aseg_disp_roicodes': aseg ROI codes corresponding to aseg_dispvec
%     {default = 30-element list matching aseg_dispvec}
%   'aseg_dispvec': per-ROI MD dispersion values (MAD estimates); snapshot
%     of MMPS atlases/DTI_ROI_Dispersion/DTI_aseg_dispersion.mat
%     {default = 30-element list matching aseg_disp_roicodes}
%   'xcg_flag': [0|1] exclude CSF and gray matter from fiber ROIs
%     {default = 1}
%   'thresh_prob': probability threshold applied to atlas-derived fiber
%     ROIs when extracting fiber-tract measures
%     {default = 0.08}
%   'inner_dti_flag': [0|1] whether to run inner-shell DTI tensor analysis
%     {default = 1}
%   'inner_bval_max': maximum b-value used for inner-shell DTI fit
%     {default = 1100}
%   'dt_outfix_inner': output file suffix for inner-shell DTI
%     {default = 'isdti'}
%   'gwcsurf_outfix_inner': gwcsurf output stem suffix for inner-shell DTI
%     {default = 'gwcsurf_inner'}
%   'dt_outfix_full': output file suffix for full-shell DTI
%     {default = 'dti_full'}
%   'rsi_outfix': output file suffix for RSI
%     {default = 'res_hind_free'}
%   'dt_measlist_full': cell array of DTI measures (full-shell)
%     {default = {'FA','MD','LD','TD'}}
%   'dt_scalefacts_full': scale factors for DTI measures
%     {default = [1,1000,1000,1000]}
%   'reregister_flag': [0|1] register the dMRI to the T1 rather than using the
%     matrix supplied in the sidecar; also done when no matrix is supplied
%     {default = 0}
%   'wmbc_flag': [0|1] white matter bias correction of the RSI input
%     {default = 1}
%   'wmbc_vol_target': target white matter intensity for the correction
%     a fixed value, so the scale is the same for every subject, which the
%     non-normalized measures depend on
%     {default = 1000}
%   'wmbc_image_range': range the corrected image is clipped to
%     {default = [0,25000]}
%   'wmbc_wm_codes': ribbon.mgz codes for the white matter mask
%     {default = [2,41]}
%   'wmbc_gm_codes': ribbon.mgz codes for the gray matter mask
%     {default = [3,42]}
%   'wmbc_excl_codes': aseg codes excluded from the white matter mask,
%     the subcortical gray and ventricle labels the ribbon includes
%   'rsi_measlist': cell array of RSI measures
%     {default = {'N0','N0s2','ND','NDs2','NT','NIf','NTs2'}}
%   'motion_radius': head radius (mm) for rotation-to-distance conversion
%     {default = 50}
%   'fibers': fiber numbers to track from the fiber atlas
%     {default = [101:110,115:123,133:138,141:150]}
%   'trans_fibers': fiber numbers to transform and analyze, including the
%     subdivisions and combinations that dti_AtlasTrack derives from 'fibers'
%     {default = [101:110,115:123,133:138,141:150,1014,1024,2000:2004]}
%   'fiber_atlasdir': full path of atlas fiber directory
%     if empty, set from $CMIG_PIPELINE_ATLASDIR/DTI_Atlas/AllSubjects
%     {default = []}
%   'atlasdir': full path of atlas directory
%     if empty, set from $CMIG_PIPELINE_ATLASDIR
%     {default = []}
%   'tabdir': directory containing ROI abbreviation lookup table CSV files
%     if empty, uses $CMIG_PIPELINE_DIR/parms/tabulated
%     {default = []}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     will have form of {output_rootdir}/{participant_id}/{session_id}/dwi
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string
%
% Created:  06/08/2026 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 06/09/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 07/14/2026 by Don Hagler
% Prev Mod: 07/22/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/20/2026 by Don Hagler
% Prev Mod: 08/24/2026 by Don Hagler
% Prev Mod: 08/25/2026 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Prev Mod: 09/03/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/28/2026 by Don Hagler
% Prev Mod: 10/08/2026 by Don Hagler
% Prev Mod: 10/17/2026 by Don Hagler
% Prev Mod: 10/19/2026 by Don Hagler
% Prev Mod: 10/20/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/12/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,5), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

% check parameters
parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin);

% locate input data
[parms,errcode,errmsg] = check_input(parms);
if errcode
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate);
  return;
end

% check whether processing required
[parms,procflag] = check_output(parms);

% process data
if procflag
  % load DWI data once for all model fits
  [vol_all,qmat_all,bvals_all,M,load_errcode,load_errmsg] = load_dwi_data(parms);
  if load_errcode
    errcode = 1; errmsg = cmig_append_errmsg(errmsg,load_errmsg);
    cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate);
    return;
  end

  % register to the T1 before any step consumes the matrix, so DTI, RSI, the
  %   aseg sampling and gwcsurf all use the same one
  [parms,reg_errcode,reg_errmsg] = reregister_to_T1(parms,vol_all,bvals_all,M);
  if reg_errcode
    errcode = 1; errmsg = cmig_append_errmsg(errmsg,reg_errmsg);
    cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate);
    return;
  end

  % white matter bias correction, once, for all three fits
  % NOTE: all three fit the same corrected volume, so one average b=0 and one
  %   brain mask serve them all. On failure everything falls back to vol_all
  [vol_fit,wmbc_errcode,wmbc_errmsg] = wmbc_correct(parms,vol_all,bvals_all,M);
  if wmbc_errcode
    fprintf('%s: WARNING: white matter bias correction failed: %s\n',...
      mfilename,wmbc_errmsg);
    fprintf('%s: WARNING: falling back to uncorrected data\n',mfilename);
    vol_fit = vol_all;
    errmsg = cmig_append_errmsg(errmsg,sprintf('white matter bias correction failed: %s',wmbc_errmsg));
  end
  parms.wmbc_applied_flag = ~wmbc_errcode && parms.wmbc_flag;
  % the non-normalized RSI measures are not meaningful on uncorrected data,
  %   whether the correction failed or was never asked for
  if ~parms.wmbc_applied_flag
    parms.rsi_measlist = normalized_measlist(parms.rsi_measlist);
  end

  % the brain mask, from the volume the fits actually see, exported as
  %   _desc-brain_mask.nii.gz beside the coefficients; the fits, coefficients and measures
  %   themselves are not masked
  [parms,mask_errcode,mask_errmsg] = get_shared_mask(parms,vol_fit,bvals_all,M);
  if mask_errcode
    fprintf('%s: WARNING: %s\n',mfilename,mask_errmsg);
    fprintf('%s: WARNING: no brain mask will be exported\n',mfilename);
    errmsg = cmig_append_errmsg(errmsg,mask_errmsg);
  end
  % NOTE: an all-ones volmask fits every voxel: dti_fit_tensor and rsi_fit_MFOD
  %   fit only inside their (dilated) mask, and dilate_brain_mask leaves an
  %   all-ones mask unchanged; given no volmask, they would build their own
  parms.fitmask = ones(size(vol_fit,1),size(vol_fit,2),size(vol_fit,3));

  % run DTI tensor fitting (full shell)
  dti_errcode = 0;
  if parms.dti_flag
    [dti_errcode,dti_errmsg] = calc_dti(parms,vol_fit,qmat_all,bvals_all,M,'dti');
    if dti_errcode
      fprintf('%s: WARNING: DTI fitting failed: %s\n',mfilename,dti_errmsg);
    end
    errmsg = cmig_append_errmsg(errmsg,dti_errmsg);
  end

  % run inner-shell DTI tensor fitting
  inner_dti_errcode = 0;
  if parms.inner_dti_flag
    [inner_dti_errcode,inner_dti_errmsg] = calc_dti(parms,vol_fit,qmat_all,bvals_all,M,'dti_inner');
    errmsg = cmig_append_errmsg(errmsg,inner_dti_errmsg);
  end

  % run RSI fitting
  rsi_errcode = 0;
  if parms.rsi_flag
    [rsi_errcode,rsi_errmsg] = calc_rsi(parms,vol_fit,qmat_all,bvals_all,M);
    if rsi_errcode
      fprintf('%s: WARNING: RSI fitting failed: %s\n',mfilename,rsi_errmsg);
    end
    errmsg = cmig_append_errmsg(errmsg,rsi_errmsg);
  end
  clear vol_fit;

  % fail if all enabled analyses failed
  if parms.dti_flag && parms.rsi_flag && dti_errcode && rsi_errcode
    if ~parms.inner_dti_flag || inner_dti_errcode
      errcode = 1;
      errmsg = cmig_append_errmsg(errmsg,'DTI and RSI failed');
      cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate);
      return;
    else
      fprintf('%s: WARNING: DTI full-shell and RSI failed but inner-shell DTI succeeded\n',mfilename);
    end
  end

  % run atlas fiber tracking (full-shell DTI only)
  if parms.fiber_flag && ~dti_errcode
    errmsg = cmig_append_errmsg(errmsg,atlas_track(parms));
    % export fiber-derived volumes so they survive tmpdir cleanup
    errmsg = cmig_append_errmsg(errmsg,export_fseg(parms));
    errmsg = cmig_append_errmsg(errmsg,export_fiber_maps(parms));
  end

  % extract ROI measures from model fits
  if ~dti_errcode
    errmsg = cmig_append_errmsg(errmsg,analyze_dwi_meas(parms,'dti'));
  end
  if parms.inner_dti_flag && ~inner_dti_errcode
    errmsg = cmig_append_errmsg(errmsg,analyze_dwi_meas(parms,'dti_inner'));
  end
  if ~rsi_errcode
    errmsg = cmig_append_errmsg(errmsg,analyze_dwi_meas(parms,'rsi'));
  end

  % gwcsurf analysis, sphere concat, and GIfTI conversion — per modality
  if parms.gwcsurf_flag
    % the surface coverage of the dMRI series, shared by every measure
    [parms,cov_errcode,cov_errmsg] = make_coverage(parms,vol_all,M);
    if cov_errcode
      fprintf('%s: WARNING: %s\n',mfilename,cov_errmsg);
      fprintf('%s: WARNING: smoothing within the cortex label, ROIs not restricted to the coverage\n',mfilename);
    end
    errmsg = cmig_append_errmsg(errmsg,cov_errmsg);
    if ~dti_errcode
      errmsg = cmig_append_errmsg(errmsg,analyze_gwcsurf(parms,'dti'));
      [errcode,tmp_errmsg] = concat_dwi_measures(parms,'dti');
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end
      errmsg = cmig_append_errmsg(errmsg,convert_dwi_measures(parms,'dti'));
    end
    if parms.inner_dti_flag && ~inner_dti_errcode
      errmsg = cmig_append_errmsg(errmsg,analyze_gwcsurf(parms,'dti_inner'));
      [errcode,tmp_errmsg] = concat_dwi_measures(parms,'dti_inner');
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end
      errmsg = cmig_append_errmsg(errmsg,convert_dwi_measures(parms,'dti_inner'));
    end
    if ~rsi_errcode
      errmsg = cmig_append_errmsg(errmsg,analyze_gwcsurf(parms,'rsi'));
      [errcode,tmp_errmsg] = concat_dwi_measures(parms,'rsi');
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end
      errmsg = cmig_append_errmsg(errmsg,convert_dwi_measures(parms,'rsi'));
    end
  end

  % write BIDS TSV tabulated data
  errmsg = cmig_append_errmsg(errmsg,write_tabdat(parms));

  % write dmri motion QC tabulated data
  cmig_write_tabdat_dmri_motion('fnames_data',parms.fnames_data,...
    'series_type',parms.series_type,'motion_radius',parms.motion_radius,...
    'outdir',parms.outdir,'outstem',parms.outstem,...
    'participant_id',parms.participant_id,'session_id',parms.session_id,...
    'forceflag',parms.forceflag);
end

% cleanup tmpdirs
if parms.cleanup_flag
  errmsg = cmig_append_errmsg(errmsg,cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir));
end

% a tolerated failure from a resource error (e.g. out of memory) fails the
%   module, so the job is rerun
if ~errcode && cmig_is_resource_error(errmsg), errcode = 1; end

% write status file
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate,procflag);

outdir = parms.outdir;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,options)
  options = cmig_parse_deployed_args(options);
  parms = mmil_args2parms(options,{...
    'participant_id',participant_id,[],...
    'session_id',session_id,[],...
    'proc_rootdir',proc_rootdir,[],...
    'fsurf_rootdir',fsurf_rootdir,[],...
    'output_rootdir',output_rootdir,[],...
  ...
    'dti_flag',true,[false true],...
    'rsi_flag',true,[false true],...
    'inner_dti_flag',true,[false true],...
    'inner_bval_max',1100,[100 Inf],...
    'dt_outfix_inner','isdti',[],...
    'gwcsurf_outfix_inner','gwcsurf_inner',[],...
    'fiber_flag',true,[false true],...
    'aseg_flag',true,[false true],...
    'aseg_disp_flag',true,[false true],...
    'aseg_dispfact',4.7,[1e-6,1e6],...
    'aseg_disp_suffix','dwtd',[],...
    'aseg_disp_roicodes',[2 4 5 7 8 10 11 12 13 14 15 16 17 18 24 26 28 41 43 44 46 47 49 50 51 52 53 54 58 60],[],...
    'aseg_dispvec',[0.07614465 0.61517701 0.37268076 0.07707907 0.16776379 0.09967775 0.16673911 0.04700424 0.06402047 0.51972693 0.85810881 0.35213536 0.20591846 0.19182216 0.59725660 0.07338498 0.20371076 0.07614465 0.61517701 0.37268076 0.07707907 0.16776379 0.09967775 0.16673911 0.04700424 0.06402047 0.20591846 0.19182216 0.07338498 0.20371076],[],...
    'xcg_flag',true,[false true],...
    'thresh_prob',0.08,[0,1],...
    'dt_outfix_full','dti_full',[],...
    'rsi_outfix','res_hind_free',[],...
    'dt_measlist_full',{'FA','MD','LD','TD'},[],...
    'dt_scalefacts_full',[1,1000,1000,1000],[],...
    'rsi_measlist',{'N0','N0s2','ND','NDs2','NT','NIf','NTs2',...
                    'F0','F0s2','FD','FDs2','FT','FTs2','If'},[],...
    'reregister_flag',false,[false true],...
    'wmbc_flag',true,[false true],...
    'volmask',[],[],...
    'wmbc_applied_flag',false,[false true],...
    'rsi_sentinel','N0',[],...
    'wmbc_vol_target',1000,[1,1e6],...
    'wmbc_image_range',[0,25000],[],...
    'wmbc_wm_codes',[2 41],[],...
    'wmbc_gm_codes',[3 42],[],...
    'wmbc_excl_codes',[4 5 10 11 12 13 14 15 17 18 24 26 28 43 44 49 50 51 52 53 54 58 60],[],...
    'motion_radius',50,[1,100],...
    'fibers',[101:110,115:123,133:138,141:150],[],...
    'trans_fibers',[101:110,115:123,133:138,141:150,1014,1024,2000:2004],[],...
    'fiber_atlasdir',[],[],...
    'atlasdir',[],[],...
    'tabdir',[],[],...
    'gwcsurf_flag',true,[false true],...
    'ico_order',5,[0 7],...
    'aparc_infix',{'aparc','aparc.a2009s'},[],...
    'gwcsurf_outfix','gwcsurf',[],...
    'sphere_flag',true,[false true],...
    'sphsmoothsteps',[0,16,256,1000],[],...
    'coverage_thresh',0.9,[0,1],...
    'fnames_coverage_label',[],[],...
    'fnames_coverage_mask',[],[],...
    'gwmask_flag',true,[false true],...
    'pmin',0.8,[0,10],...
    'pmax',2.0,[0,10],...
    'pstep',0.2,[0.001,10],...
    'tukey_flag',true,[false true],...
    'tukey_fact',0.5,[1e-10,1e10],...
    'cleanup_flag',true,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'modality','dwi',[],...
    'series_type','desc-preproc_dwi',[],...
    'status_extra',[],[],...
  });

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];

  % set output directory
  parms.outdir = sprintf('%s/%s/%s/%s',parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);

  % enumerate DWI NIfTI files from proc_rootdir
  [parms.fnames_data,parms.fnames_info] = get_dwi_data(parms);
  if isempty(parms.fnames_data)
    errcode = 1;
    errmsg = sprintf('no DWI files found in %s/%s/%s/%s',...
      parms.proc_rootdir,parms.participant_id,parms.session_id,parms.modality);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % derive paired bvec/bval paths
  nscans = length(parms.fnames_data);
  parms.fnames_bval = cell(nscans,1);
  parms.fnames_bvec = cell(nscans,1);
  for s = 1:nscans
    parms.fnames_bval{s} = regexprep(parms.fnames_data{s},'\.nii.+','.bval');
    parms.fnames_bvec{s} = regexprep(parms.fnames_data{s},'\.nii.+','.bvec');
    if ~exist(parms.fnames_bval{s},'file')
      errcode = 1;
      errmsg = sprintf('bval file not found: %s',parms.fnames_bval{s});
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    if ~exist(parms.fnames_bvec{s},'file')
      errcode = 1;
      errmsg = sprintf('bvec file not found: %s',parms.fnames_bvec{s});
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
  end

  % check FreeSurfer (always required for ROI extraction)
  [parms.fspath,parms.fs_version,parms.fs_build,parms.fs_complete_flag,errcode,errmsg] = ...
    cmig_check_bids_fsurf(parms.participant_id,parms.session_id,parms.fsurf_rootdir);
  parms.status_extra = struct('fs_version',parms.fs_version,...
    'fs_build',parms.fs_build,'fs_complete_flag',parms.fs_complete_flag);
  if errcode
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  % volume-only recon is fatal: surface sampling requires a complete recon
  if ~parms.fs_complete_flag
    errcode = 1;
    errmsg = sprintf('not processed: freesurfer recon %s is volume-only; surface sampling requires a complete recon',parms.fspath);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % read RegistrationMatrixT1w from first valid JSON sidecar
  % RegistrationMatrixT1w is T1-to-dwi, as cmig_proc_bids_func writes it and
  %   cmig_check_regT1_epi reads it
  % NOTE: reg_supplied_flag records whether a matrix was actually found, which
  %   identity alone cannot tell: it is the initial value, it is what ABCD
  %   sidecars legitimately contain, and it is what a missing or unparseable
  %   field leaves behind. Only the flag distinguishes registered with no
  %   transform needed from never registered
  parms.M_T1_to_DTI = eye(4);
  parms.reg_supplied_flag = 0;
  for s = 1:nscans
    if isempty(parms.fnames_info{s}) || ~exist(parms.fnames_info{s},'file')
      continue;
    end
    try
      json = jsondecode(fileread(parms.fnames_info{s}));
      M = cmig_get_reg_T1w(json);
      if ~isempty(M)
        if isnumeric(M) && isequal(size(M),[4,4])
          parms.M_T1_to_DTI = M;
          parms.reg_supplied_flag = 1;
          break;
        end
      end
    catch
    end
  end
  if ~parms.reg_supplied_flag
    fprintf('%s: WARNING: no registration matrix in the sidecars for %s %s\n',...
      mfilename,parms.participant_id,parms.session_id);
  end

  % check atlas fiber directory; disable fiber tracking if absent
  if isempty(parms.fiber_atlasdir)
    parms.fiber_atlasdir = sprintf('%s/DTI_Atlas/AllSubjects',cmig_getenv('CMIG_PIPELINE_ATLASDIR'));
  end
  if ~exist(parms.fiber_atlasdir,'dir')
    fprintf('%s: WARNING: atlas fiber dir not found: %s; fiber_flag set to 0\n',...
      mfilename,parms.fiber_atlasdir);
    errmsg = sprintf('atlas fiber dir not found: %s',parms.fiber_atlasdir);
    parms.fiber_flag = 0;
  end

  % set atlasdir
  if isempty(parms.atlasdir)
    parms.atlasdir = cmig_getenv('CMIG_PIPELINE_ATLASDIR');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_data,fnames_info] = get_dwi_data(parms)
  % the proc dwi named by series_type (desc-preproc_dwi); if a session has
  %   none, the plain _dwi name of the released mproc derivatives
  [fnames_data,fnames_info] = cmig_get_bids_data(...
    parms.participant_id,parms.session_id,parms.proc_rootdir,...
    parms.modality,parms.series_type);
  if ~isempty(fnames_data)
    fprintf('%s: using %s files\n',mfilename,parms.series_type);
    return;
  end
  old_type = regexprep(parms.series_type,'^desc-[^_]+_','');
  if strcmp(old_type,parms.series_type), return; end
  [fnames_data,fnames_info] = cmig_get_bids_data(...
    parms.participant_id,parms.session_id,parms.proc_rootdir,...
    parms.modality,old_type);
  if isempty(fnames_data), return; end
  % keep only {stem}_{old_type}.nii[.gz] with no desc- entity
  pat = sprintf('_%s\\.nii(\\.gz)?$',old_type);
  idx = find(~cellfun(@isempty,regexp(fnames_data,pat,'once')) & ...
    cellfun(@isempty,regexp(fnames_data,'_desc-','once')));
  fnames_data = fnames_data(idx);
  fnames_info = fnames_info(idx);
  if ~isempty(fnames_data)
    fprintf('%s: no %s files; using %s files\n',...
      mfilename,parms.series_type,old_type);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnames = expected_outputs(parms)
  % the full paths of every output the enabled analyses write, named as their
  %   writers name them: export_coef, export_meas, export_meas_resT1,
  %   export_fseg, export_fiber_maps, write_tabdat, convert_dwi_measures,
  %   cmig_concat_dwi_measures and cmig_write_tabdat_dmri_motion
  stem = sprintf('%s/%s',parms.outdir,parms.outstem);
  fnames = {};
  layers = {'gm','wm','gwc'};
  fwhms = round(1.25 * sqrt(parms.sphsmoothsteps));
  % the RSI measures: the normalized ones alone where the bias correction is
  %   off, or failed in the run that wrote the coefficients
  rsi_measlist = parms.rsi_measlist;
  if ~parms.wmbc_flag
    rsi_measlist = normalized_measlist(rsi_measlist);
  else
    fname_info = sprintf('%s_rsicoef.json',stem);
    if exist(fname_info,'file')
      try
        info = jsondecode(fileread(fname_info));
        if ~mmil_getfield(info,'WMBiasCorrection',true)
          rsi_measlist = normalized_measlist(rsi_measlist);
        end
      catch
      end
    end
  end
  families = {};
  if parms.dti_flag
    families{end+1} = {'dti',parms.dt_measlist_full};
  end
  if parms.inner_dti_flag
    families{end+1} = {'isdti',parms.dt_measlist_full};
  end
  if parms.rsi_flag && ~isempty(rsi_measlist)
    ext = cellfun(@get_rsi_extmeas,rsi_measlist,'UniformOutput',false);
    families{end+1} = {'rsi',ext};
  end
  if ~isempty(families)
    fnames{end+1} = sprintf('%s_desc-brain_mask.nii.gz',stem);
    fnames{end+1} = sprintf('%s_desc-brain_mask.json',stem);
    fnames{end+1} = sprintf('%s_dwiref.nii.gz',stem);
    fnames{end+1} = sprintf('%s_dwiref.json',stem);
    % the coverage mask paired with the ico7 GIfTI files (make_coverage)
    if parms.gwcsurf_flag
      for hemi={'L','R'}
        fnames{end+1} = sprintf('%s_hemi-%s_atlas-ico7_mask.func.gii',stem,hemi{1});
        fnames{end+1} = sprintf('%s_hemi-%s_atlas-ico7_mask.json',stem,hemi{1});
      end
    end
  end
  for f=1:length(families)
    family = families{f}{1};
    measlist = families{f}{2};
    % the coefficients and the measure volumes
    fnames{end+1} = sprintf('%s_%scoef.nii.gz',stem,family);
    fnames{end+1} = sprintf('%s_%scoef.json',stem,family);
    if strcmp(family,'rsi')
      vols = measlist;
    else
      vols = get_dti_export_meas();
    end
    for m=1:length(vols)
      fnames{end+1} = sprintf('%s_%s%s.nii.gz',stem,family,vols{m});
      fnames{end+1} = sprintf('%s_%s%s.json',stem,family,vols{m});
    end
    % the ROI tables
    for m=1:length(measlist)
      meas = measlist{m};
      if parms.aseg_flag
        fnames{end+1} = sprintf('%s_atlas-aseg_%s%s.tsv',stem,family,meas);
      end
      if parms.fiber_flag
        fnames{end+1} = sprintf('%s_atlas-at_%s%s.tsv',stem,family,meas);
      end
      if parms.gwcsurf_flag
        for l=1:length(layers)
          for atlas={'dsk','dst'}
            fnames{end+1} = sprintf('%s_atlas-%s_%s%s%s.tsv',...
              stem,atlas{1},family,meas,layers{l});
          end
          for k=1:length(fwhms)
            for hemi={'L','R'}
              fnames{end+1} = sprintf('%s_hemi-%s_atlas-ico7_desc-fwhm%d_%s%s%s.func.gii',...
                stem,hemi{1},fwhms(k),family,meas,layers{l});
              fnames{end+1} = sprintf('%s_hemi-%s_atlas-ico7_desc-fwhm%d_%s%s%s.json',...
                stem,hemi{1},fwhms(k),family,meas,layers{l});
              fnames{end+1} = sprintf('%s_hemi-%s_atlas-ico%d_desc-fwhm%d_%s%s%s.mat',...
                stem,hemi{1},parms.ico_order,fwhms(k),family,meas,layers{l});
            end
          end
        end
      end
    end
  end
  if parms.rsi_flag && ismember('ND',rsi_measlist)
    fnames{end+1} = sprintf('%s_space-T1w_rsiRND.nii.gz',stem);
    fnames{end+1} = sprintf('%s_space-T1w_rsiRND.json',stem);
  end
  if parms.fiber_flag && parms.dti_flag
    fnames{end+1} = sprintf('%s_space-T1w_atlas-at_dseg.nii.gz',stem);
    fnames{end+1} = sprintf('%s_space-T1w_atlas-at_dseg.json',stem);
    fnames{end+1} = sprintf('%s_atlas-at_probseg.nii.gz',stem);
    fnames{end+1} = sprintf('%s_atlas-at_probseg.json',stem);
  end
  fnames{end+1} = sprintf('%s_motionqc.tsv',stem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;

  % set output file stem and temporary directory
  parms.outstem = sprintf('%s_%s',parms.participant_id,parms.session_id);
  parms.tmpdir = sprintf('%s/analyze_dwi',parms.outdir);

  % set intermediate file stems (no BIDS entities — these are internal tmpdir files)
  parms.dt_fstem_full  = sprintf('%s/%s_%s',parms.tmpdir,parms.outstem,parms.dt_outfix_full);
  parms.dt_fstem_inner = sprintf('%s/%s_%s',parms.tmpdir,parms.outstem,parms.dt_outfix_inner);
  parms.rsi_fstem      = sprintf('%s/%s_%s',parms.tmpdir,parms.outstem,parms.rsi_outfix);
  parms.at_outdir      = sprintf('%s/AtlasTrack',parms.tmpdir);

  % check whether processing is required
  if ~exist(parms.outdir,'dir') || parms.forceflag
    procflag = 1;
  else
    % every output the enabled analyses write: one missing means a prior run
    %   was incomplete, an output was deleted, or new outputs were added, and
    %   the writers skip the outputs that exist
    fnames_check = expected_outputs(parms);
    for fi = 1:length(fnames_check)
      if ~exist(fnames_check{fi},'file')
        procflag = 1;
        break;
      end
    end
  end

  % create temporary directory (and output directory)
  if procflag, mmil_mkdir(parms.tmpdir); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol_all,qmat_all,bvals_all,M,errcode,errmsg] = load_dwi_data(parms)
  vol_all = []; qmat_all = []; bvals_all = []; M = [];
  errcode = 0; errmsg = [];
  for s = 1:length(parms.fnames_data)
    fname_mgz = sprintf('%s/run%02d_dwi.mgz',parms.tmpdir,s);
    if ~exist(fname_mgz,'file') || parms.forceflag
      try
        fs_mri_convert(parms.fnames_data{s},fname_mgz);
      catch err
        errcode = 1;
        errmsg = sprintf('fs_mri_convert failed for %s: %s',parms.fnames_data{s},err.message);
        return;
      end
    end
    [vol_s,M_s] = fs_load_mgh(fname_mgz);
    if isempty(vol_s)
      errcode = 1;
      errmsg = sprintf('failed to load %s',fname_mgz);
      return;
    end
    if isempty(M), M = M_s; end
    bvec = load(parms.fnames_bvec{s});
    if size(bvec,1) ~= 3, bvec = bvec'; end
    bval = load(parms.fnames_bval{s});
    bval = bval(:);
    % qmat is in RAS world, as MMPS's qmat is
    if mmps_bvec_flag(parms,s)
      % the release mproc bvecs are MMPS's RAS qmat with x negated
      %   (mmil_write_bfiles), not BIDS bvecs for their image axes
      fprintf('%s: %s: MMPS bvecs (registration_matrix_T1 in the sidecar); x negated\n',...
        mfilename,parms.fnames_bvec{s});
      bvec(1,:) = -bvec(1,:);
      qmat = bvec';
    else
      % BIDS bvecs are relative to the image axes, with x negated when the
      %   affine's 3x3 has a positive determinant (BIDS/FSL convention)
      R = M_s(1:3,1:3);
      R = bsxfun(@rdivide,R,sqrt(sum(R.^2,1)));
      if det(R) > 0, bvec(1,:) = -bvec(1,:); end
      qmat = (R*bvec)';
    end
    vol_all = cat(4,vol_all,double(vol_s));
    qmat_all = cat(1,qmat_all,qmat);
    bvals_all = cat(1,bvals_all,bval);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flag = mmps_bvec_flag(parms,s)
  % true for a scan of the release mproc derivatives: its sidecar has the
  %   field name used before 08/2026, registration_matrix_T1, which only
  %   those derivatives carry (cmig_get_reg_T1w)
  flag = false;
  if s > length(parms.fnames_info) || isempty(parms.fnames_info{s}) || ...
     ~exist(parms.fnames_info{s},'file')
    return;
  end
  try
    json = jsondecode(fileread(parms.fnames_info{s}));
  catch
    return;
  end
  flag = isfield(json,'registration_matrix_T1') && ...
    ~isfield(json,'RegistrationMatrixT1w');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = calc_dti(parms, vol_all, qmat_all, bvals_all, M, modality)
  errcode = 0; errmsg = [];
  dti_inner = strcmp(modality,'dti_inner');
  if dti_inner
    fprintf('%s: running inner-shell DTI tensor fitting (max_bval=%d)...\n',...
      mfilename,parms.inner_bval_max);
    fstem = parms.dt_fstem_inner;
    max_bval = parms.inner_bval_max;
    family = 'isdti';
  else
    fprintf('%s: running DTI tensor fitting...\n',mfilename);
    fstem = parms.dt_fstem_full;
    max_bval = Inf;
    family = 'dti';
  end
  % fit DTI tensor
  try
    DTfit = dti_fit_tensor(vol_all,qmat_all,'bvals',bvals_all,'M',M,...
      'max_bval',max_bval,'volmask',parms.fitmask);
  catch err
    errcode = 1; errmsg = sprintf('dti_fit_tensor failed: %s',err.message);
    if dti_inner, fprintf('%s: WARNING: %s\n',mfilename,errmsg); end
    return;
  end
  if isempty(DTfit)
    errcode = 1; errmsg = 'dti_fit_tensor returned empty result';
    if dti_inner, fprintf('%s: WARNING: %s\n',mfilename,errmsg); end
    return;
  end
  % NOTE: cast to single before the measures are computed, so the shipped
  %   measures are exactly what a consumer recomputes from the exported
  %   coefficients, which are written as single. The cast moves the measures by
  %   a few ULP at float32 and nothing more
  DTfit.volDT = single(DTfit.volDT);
  % calculate DTI measures (FA, MD, LD, TD, b0, V0)
  try
    DTmeas = dti_calc_DTmeas(DTfit,0);
  catch err
    errcode = 1; errmsg = sprintf('dti_calc_DTmeas failed: %s',err.message);
    if dti_inner, fprintf('%s: WARNING: %s\n',mfilename,errmsg); end
    return;
  end
  % save MGZ files (b0, FA, V0, MD, LD, TD)
  try
    dti_convert_DTmeas(DTmeas,fstem,'M',M,...
      'measlist',{'b0','FA','V0','MD','LD','TD'},'forceflag',parms.forceflag);
  catch err
    errcode = 1; errmsg = sprintf('dti_convert_DTmeas failed: %s',err.message);
    if dti_inner, fprintf('%s: WARNING: %s\n',mfilename,errmsg); end
    return;
  end
  % export the coefficients and the measure volumes to the BIDS output
  errmsg = cmig_append_errmsg(errmsg,export_coef(parms,DTfit,family));
  errmsg = cmig_append_errmsg(errmsg,export_meas(parms,fstem,get_dti_export_meas(),family));
  if dti_inner
    fprintf('%s: inner-shell DTI tensor fitting done\n',mfilename);
  else
    fprintf('%s: DTI tensor fitting done\n',mfilename);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_b0,errcode,errmsg] = get_avg_b0(parms,vol_all,bvals_all,M,outdir)
  % average of the b=0 frames, the T2-weighted image both the registration and
  %   the bias correction work from
  fname_b0 = []; errcode = 0; errmsg = [];
  % b <= 10 is b=0, as in dti_fit_tensor (b0_thresh), so the Philips
  %   b=0.001 frames count
  idx_b0 = find(bvals_all<=10);
  if isempty(idx_b0)
    errcode = 1; errmsg = 'no b=0 frames found'; return;
  end
  fname_b0 = sprintf('%s/avg_b0.mgz',outdir);
  if ~exist(fname_b0,'file') || parms.forceflag
    mmil_mkdir(outdir);
    fs_save_mgh(mean(vol_all(:,:,:,idx_b0),4),fname_b0,M);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = get_shared_mask(parms,vol_fit,bvals_all,M)
  % one brain mask for every fit, from the average b=0 of the volume they are
  %   given. Passed as volmask so dti_fit_tensor and rsi_fit_MFOD skip their own
  %   mask creation; an empty mask leaves each to build its own, as before
  % NOTE: the b=0 and the mask are named by whether the correction was applied,
  %   so a rerun that falls back cannot pick up the corrected pair cached by an
  %   earlier run that succeeded, or the reverse
  errcode = 0; errmsg = [];
  parms.volmask = [];
  if parms.wmbc_applied_flag
    maskdir = sprintf('%s/mask_wmbc',parms.tmpdir);
  else
    maskdir = sprintf('%s/mask',parms.tmpdir);
  end
  [fname_b0,errcode,errmsg] = get_avg_b0(parms,vol_fit,bvals_all,M,maskdir);
  if errcode, return; end
  fname_mask = sprintf('%s/brainmask.mgz',maskdir);
  if ~exist(fname_mask,'file') || parms.forceflag
    try
      cmig_create_brainmask(fname_b0,'fname_mask',fname_mask,...
        'binarize_flag',true,'forceflag',parms.forceflag);
    catch err
      errcode = 1;
      errmsg = sprintf('cmig_create_brainmask failed: %s',err.message);
      return;
    end
  end
  parms.volmask = fs_load_mgh(fname_mask);
  if isempty(parms.volmask)
    errcode = 1; errmsg = sprintf('failed to load %s',fname_mask);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function measlist = normalized_measlist(measlist)
  % the RSI measures that survive without the bias correction, whether it
  %   failed or was never asked for. The normalized ones are ratios, so a
  %   multiplicative field cancels; the rest are not meaningful on uncorrected
  %   data. Internal names starting with N are the normalized ones: N0, ND, NT,
  %   NIf and their s{N} frames
  keep = false(1,length(measlist));
  for m=1:length(measlist)
    keep(m) = strncmp(measlist{m},'N',1);
  end
  if ~all(keep)
    fprintf('%s: WARNING: dropping non-normalized RSI measures: %s\n',...
      mfilename,strjoin(measlist(~keep),' '));
  end
  measlist = measlist(keep);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = reregister_to_T1(parms,vol_all,bvals_all,M)
  % register the dMRI to the T1 and replace parms.M_T1_to_DTI with the result
  % NOTE: run once, here, so every step downstream - wmbc, the aseg sampling and
  %   gwcsurf - uses the same matrix. A step registering for itself is how the
  %   steps come to disagree
  errcode = 0; errmsg = [];
  if ~parms.reregister_flag && parms.reg_supplied_flag, return; end
  if ~parms.reregister_flag
    fprintf('%s: NOTE: registering to the T1, since no matrix was supplied\n',...
      mfilename);
  end

  regdir = sprintf('%s/reg_T1',parms.tmpdir);
  fname_reg = sprintf('%s/regT2toT1.mat',regdir);
  if exist(fname_reg,'file') && ~parms.forceflag
    reginfo = load(fname_reg);
    parms.M_T1_to_DTI = reginfo.M_T1_to_T2;
    return;
  end
  mmil_mkdir(regdir);

  fname_T1 = sprintf('%s/mri/nu.mgz',parms.fspath);
  if ~exist(fname_T1,'file')
    errcode = 1; errmsg = sprintf('T1 for registration not found: %s',fname_T1);
    return;
  end
  [fname_b0,errcode,errmsg] = get_avg_b0(parms,vol_all,bvals_all,M,regdir);
  if errcode, return; end
  try
    fname_bm = cmig_create_aseg_brainmask(parms.fspath,regdir,...
      parms.verbose,parms.forceflag);
  catch err
    errcode = 1;
    errmsg = sprintf('cmig_create_aseg_brainmask failed: %s',err.message);
    return;
  end

  p = [];
  p.fname_T1_mask = fname_bm;
  p.outdir = regdir;
  p.atlasdir = parms.atlasdir;
  p.cleanup_flag = parms.cleanup_flag;
  p.forceflag = parms.forceflag;
  a = mmil_parms2args(p);
  try
    M_T1_to_T2 = mmil_reg_T2toT1(fname_b0,fname_T1,a{:});
  catch err
    errcode = 1; errmsg = sprintf('mmil_reg_T2toT1 failed: %s',err.message);
    return;
  end
  if isempty(M_T1_to_T2)
    errcode = 1; errmsg = 'mmil_reg_T2toT1 returned no matrix'; return;
  end
  save(fname_reg,'M_T1_to_T2','fname_T1','fname_b0');
  parms.M_T1_to_DTI = M_T1_to_T2;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol_out,errcode,errmsg] = wmbc_correct(parms,vol_all,bvals_all,M)
  % white matter bias correction of the DWI series, for the RSI fit
  % NOTE: a single multiplicative field, estimated once from the average b=0
  %   frame and applied to every frame, so the diffusion weighting is untouched.
  %   Follows MMIL_wmCorr_DTI: average b0, mmil_wm_corr with the T1 supplied,
  %   ratio of corrected to original, multiply through
  vol_out = vol_all; errcode = 0; errmsg = [];
  if ~parms.wmbc_flag
    fprintf('%s: NOTE: white matter bias correction disabled\n',mfilename);
    return;
  end

  % NOTE: no cast on the way back in. fs_load_mgh returns double, keepsingle
  %   defaulting to 0, holding the float32 values the file stores -- which is
  %   what the product below is rounded to, so a rerun reading this file fits
  %   the same values as the run that wrote it
  fname_out = sprintf('%s/%s_wmbc.mgz',parms.tmpdir,parms.outstem);
  if exist(fname_out,'file') && ~parms.forceflag
    vol_out = fs_load_mgh(fname_out);
    if isempty(vol_out)
      errcode = 1; errmsg = sprintf('failed to load %s',fname_out); return;
    end
    return;
  end

  % NOTE: everything wmbc writes lives in one subdirectory, removed whole on
  %   forceflag, so a rerun cannot mix new masks with an old correction
  wmbc_dir = sprintf('%s/wmbc',parms.tmpdir);
  if parms.forceflag && exist(wmbc_dir,'dir')
    [st,msg] = mmil_unix(sprintf('rm -rf %s',wmbc_dir));
    if st
      errcode = 1; errmsg = sprintf('failed to remove %s: %s',wmbc_dir,msg); return;
    end
  end
  mmil_mkdir(wmbc_dir);

  % T1 in FreeSurfer conformed space, so the masks need no resampling
  fname_T1 = sprintf('%s/mri/nu.mgz',parms.fspath);
  if ~exist(fname_T1,'file')
    errcode = 1; errmsg = sprintf('T1 for bias correction not found: %s',fname_T1);
    return;
  end

  % segmentation masks from the recon
  [fname_wm,fname_bm,fname_gm,errcode,errmsg] = create_wmbc_masks(parms,wmbc_dir);
  if errcode, return; end

  % average the b=0 frames
  [fname_b0,errcode,errmsg] = get_avg_b0(parms,vol_all,bvals_all,M,wmbc_dir);
  if errcode, return; end

  % estimate the correction on the average b0
  % NOTE: forceflag stays false here whatever ours is. mmil_wm_corr passes it to
  %   mmil_wmseg, which returns early only when the wm and bm files exist AND
  %   forceflag is off; with it on, mmil_wmseg would segment the T1 itself and
  %   overwrite the masks written above. Recomputation is forced by removing
  %   wmbc_dir instead
  p = [];
  p.fname_T1 = fname_T1;
  p.fname_out = sprintf('%s/avg_b0_wmbc.mgz',wmbc_dir);
  p.fname_wm_T1 = fname_wm;
  p.fname_bm_T1 = fname_bm;
  p.fname_gm_T1 = fname_gm;
  p.fname_wm = sprintf('%s/wm_b0.mgz',wmbc_dir);
  p.fname_bm = sprintf('%s/bm_b0.mgz',wmbc_dir);
  p.fname_gm = sprintf('%s/gm_b0.mgz',wmbc_dir);
  % NOTE: a fixed vol_target, so white matter lands on the same intensity for
  %   every subject. Left empty, mmil_corr_wm_bias_spsm falls back to
  %   median(vol(vol_wm>0)), which rescales each subject to its own median and
  %   leaves the non-normalized measures incomparable between subjects
  p.vol_target = parms.wmbc_vol_target;
  p.image_range = parms.wmbc_image_range;
  % NOTE: the registration comes from the sidecar, read in check_input, so
  %   mmil_wm_corr skips registering and every step uses the same matrix. The b0
  %   is in dwi space, so M_T1_to_DTI is the M_T1_to_T2 it wants
  p.M_T1_to_T2 = parms.M_T1_to_DTI;
  % NOTE: passed explicitly, so mmil_reg_T2toT1 uses the same atlas
  %   directory as this module
  p.atlasdir = parms.atlasdir;
  p.tmpdir = sprintf('%s/tmp_wm_corr',wmbc_dir);
  p.cleanup_flag = parms.cleanup_flag;
  p.verbose = parms.verbose;
  p.forceflag = false;
  a = mmil_parms2args(p);
  try
    fname_b0_corr = mmil_wm_corr(fname_b0,a{:});
  catch err
    errcode = 1; errmsg = sprintf('mmil_wm_corr failed: %s',err.message); return;
  end
  if isempty(fname_b0_corr) || ~exist(fname_b0_corr,'file')
    errcode = 1; errmsg = 'mmil_wm_corr produced no output'; return;
  end

  % the bias field is the ratio of corrected to original
  vol_b0 = fs_load_mgh(fname_b0);
  vol_b0_corr = fs_load_mgh(fname_b0_corr);
  vol_ratio = double(vol_b0_corr) ./ double(vol_b0);
  vol_ratio(~isfinite(vol_ratio)) = 0;
  fs_save_mgh(vol_ratio,sprintf('%s/bias_ratio.mgz',wmbc_dir),M);

  % apply to every frame
  % NOTE: rounded through float32, which is all mgz stores, so that this run
  %   and any rerun reading the file back fit identical values. The rounding is
  %   for the round trip and for nothing else, which is why vol_all and the two
  %   paths that return it uncorrected are left alone
  % NOTE: rounded, not cast. Making the array single would put the
  %   least-squares solve in single, and rsi_fit_MFOD and dti_fit_tensor both
  %   preallocate their output as double, so the result would be widened back
  %   again with nothing gained
  vol_out = double(single(bsxfun(@times,vol_all,vol_ratio)));
  fs_save_mgh(vol_out,fname_out,M);
  fprintf('%s: white matter bias correction done\n',mfilename);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_wm,fname_bm,fname_gm,errcode,errmsg] = create_wmbc_masks(parms,outdir)
  % wm, brain and gm masks from the recon, in FreeSurfer conformed space
  % NOTE: written at the paths mmil_wm_corr hands to mmil_wmseg as fname_wm and
  %   fname_bm, which makes mmil_wmseg return early and use these rather than
  %   segmenting the T1 itself
  fname_wm = []; fname_bm = []; fname_gm = []; errcode = 0; errmsg = [];
  fname_ribbon = sprintf('%s/mri/ribbon.mgz',parms.fspath);
  fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
  fnames_seg = {fname_ribbon,fname_aseg};
  for f=1:length(fnames_seg)
    if ~exist(fnames_seg{f},'file')
      errcode = 1; errmsg = sprintf('segmentation not found: %s',fnames_seg{f});
      return;
    end
  end

  % wm from the cortical ribbon, which is bounded by the white surface and so
  %   excludes the cortex and brainstem a binarized wm.mgz includes, minus the
  %   aseg subcortical gray and ventricle labels, which the ribbon does include
  fname_wm = sprintf('%s/wm_T1.mgz',outdir);
  if ~exist(fname_wm,'file') || parms.forceflag
    [vol_wm,M_wm] = fs_aseg_mask(fname_ribbon,'aseg_codes',parms.wmbc_wm_codes);
    vol_excl = fs_aseg_mask(fname_aseg,'aseg_codes',parms.wmbc_excl_codes);
    vol_wm(vol_excl>0) = 0;
    if ~any(vol_wm(:))
      errcode = 1; errmsg = 'white matter mask is empty'; return;
    end
    fs_save_mgh(vol_wm,fname_wm,M_wm);
  end

  % gm: cortex only, from the ribbon
  fname_gm = sprintf('%s/gm_T1.mgz',outdir);
  if ~exist(fname_gm,'file') || parms.forceflag
    fs_aseg_mask(fname_ribbon,'fname_mask',fname_gm,...
      'aseg_codes',parms.wmbc_gm_codes,'forceflag',parms.forceflag);
  end

  % bm: from aseg, which is tighter than brainmask.mgz
  try
    fname_bm = cmig_create_aseg_brainmask(parms.fspath,outdir,...
      parms.verbose,parms.forceflag);
  catch err
    errcode = 1;
    errmsg = sprintf('cmig_create_aseg_brainmask failed: %s',err.message);
    return;
  end
  if isempty(fname_bm) || ~exist(fname_bm,'file')
    errcode = 1; errmsg = 'aseg brain mask not created'; return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = calc_rsi(parms, vol_all, qmat_all, bvals_all, M)
  errcode = 0; errmsg = [];
  fprintf('%s: running RSI fitting...\n',mfilename);

  % fit RSI multi-compartment FOD model
  % NOTE: these match the old pipeline's RSI configuration, DTI_RSI_* in
  %   cmig_proc_dMRIv2_full. num_ADC_trans sets the number of transverse ADC
  %   scales, so it decides which scale the s2 measures name: with 2 the
  %   values are [0 0.9e-3] and frame 2 is the hindered end, which is what
  %   HNT means. The rsi_fit_MFOD default of 5 puts frame 2 at 0.225e-3
  try
    MFfit = rsi_fit_MFOD(vol_all,qmat_all,'bvals',bvals_all,'M',M,...
      'num_ADC_trans',2,...
      'iso_restricted_flag',0,'iso_hindered_flag',0,'iso_free_flag',1,...
      'volmask',parms.fitmask);
  catch err
    errcode = 1; errmsg = sprintf('rsi_fit_MFOD failed: %s',err.message); return;
  end
  if isempty(MFfit)
    errcode = 1; errmsg = 'rsi_fit_MFOD returned empty result'; return;
  end

  % NOTE: cast to single before the measures are computed, for the same reason
  %   as in calc_dti: the shipped measures then match what a consumer
  %   recomputes from the exported coefficients
  MFfit.volMF = single(MFfit.volMF);
  % calculate RSI measures
  try
    MFmeas = rsi_calc_MFmeas(MFfit,'mask_flag',false);
  catch err
    errcode = 1; errmsg = sprintf('rsi_calc_MFmeas failed: %s',err.message); return;
  end

  % save MGZ files (exclude s{N} measures — they are frames of multi-frame base MGZs)
  base_measlist = {};
  for m = 1:length(parms.rsi_measlist)
    if isempty(regexp(parms.rsi_measlist{m},'s\d+$','once'))
      base_measlist{end+1} = parms.rsi_measlist{m};
    end
  end
  try
    rsi_convert_MFmeas(MFmeas,parms.rsi_fstem,'M',M,...
      'measlist',base_measlist,'forceflag',parms.forceflag);
  catch err
    errcode = 1; errmsg = sprintf('rsi_convert_MFmeas failed: %s',err.message); return;
  end
  % export the coefficients and the measure volumes to the BIDS output
  errmsg = cmig_append_errmsg(errmsg,export_coef(parms,MFfit,'rsi'));
  errmsg = cmig_append_errmsg(errmsg,export_meas(parms,parms.rsi_fstem,base_measlist,'rsi'));
  % the restricted directional measure on the T1w grid, for manual QC
  errmsg = cmig_append_errmsg(errmsg,export_meas_resT1(parms,parms.rsi_fstem,'ND','rsi'));
  fprintf('%s: RSI fitting done\n',mfilename);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = atlas_track(parms)
  errmsg = [];
  fprintf('%s: running atlas fiber tracking...\n',mfilename);
  fname_FA = sprintf('%s_FA.mgz',parms.dt_fstem_full);
  fname_V0 = sprintf('%s_V0.mgz',parms.dt_fstem_full);
  % nu.mgz, as for the DTI-to-T1 registration: T1.mgz carries a filter that
  %   flattens the white matter intensities
  fname_T1 = sprintf('%s/%s_%s/mri/nu.mgz',...
    parms.fsurf_rootdir,parms.participant_id,parms.session_id);
  if ~exist(fname_FA,'file')
    tmp_errmsg = sprintf('FA file not found: %s; skipping atlas tracking',fname_FA);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  if ~exist(fname_T1,'file')
    tmp_errmsg = sprintf('T1 file not found: %s; skipping atlas tracking',fname_T1);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  try
    dti_AtlasTrack(...
      'fname_T1',fname_T1,...
      'fname_FA',fname_FA,...
      'fname_V0',fname_V0,...
      'M_T1_to_DTI',parms.M_T1_to_DTI,...
      'fibers',parms.fibers,...
      'fiber_atlasdir',parms.fiber_atlasdir,...
      'atlasdir',parms.atlasdir,...
      'outdir',parms.at_outdir,...
      'create_paths_flag',0,...
      'forceflag',parms.forceflag);
  catch err
    tmp_errmsg = sprintf('dti_AtlasTrack failed: %s',err.message);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  end
  % dti_AtlasTrack only writes the base fiber files (fiber_N_prob_countatlas);
  % create the xcg-masked, prob-thresholded variants
  % (fiber_N_prob_countatlas_xcg_pthresh%0.2f) that mmil_fiber_analysis reads,
  % plus the fseg fiber segmentation volume at T1 resolution
  fiber_dir = sprintf('%s/fiber_maps',parms.at_outdir);
  if exist(fiber_dir,'dir')
    fname_aseg = sprintf('%s/%s_%s/mri/aseg.mgz',...
      parms.fsurf_rootdir,parms.participant_id,parms.session_id);
    if ~exist(fname_aseg,'file')
      tmp_errmsg = sprintf('aseg not found (%s); skipping xcg fiber masking',fname_aseg);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    else
      [M_T1,volsz_T1] = fs_read_header(fname_T1);
      volsz_T1 = volsz_T1(1:3);
      try
        dti_transform_fibers(fiber_dir,...
          'outdir',parms.at_outdir,...
          'fiber_outdir',fiber_dir,...
          'atlas_flag',2,...
          'fibers',parms.trans_fibers,...
          'xcg_flag',1,...
          'thresh_prob',parms.thresh_prob,...
          'fname_aseg',fname_aseg,...
          'M_T1_to_DTI',parms.M_T1_to_DTI,...
          'M_T1',M_T1,...
          'volsz_T1',volsz_T1,...
          'fname_FA',fname_FA,...
          'fseg_flag',1,...
          'fseg_resT1flag',1,...
          'fname_fseg','fseg_resT1_xcg.mgz',...
          'create_paths_flag',0,...
          'forceflag',parms.forceflag);
      catch err
        tmp_errmsg = sprintf('dti_transform_fibers failed: %s',err.message);
        fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
        errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      end
    end
  end
  fprintf('%s: atlas fiber tracking done\n',mfilename);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function extmeas = get_rsi_extmeas(meas)
  % RSI internal → external measure name
  % NOTE: the first seven are the ABCD 7.0 released names; the rest are new.
  %   Prefix is compartment, R restricted, H hindered, F free; N means
  %   normalized, and s2 is frame 2 of the base measure
  rsi_int = {'N0','ND','NT','NIf','NTs2','N0s2','NDs2',...
             'F0','FD','FT','If','FTs2','F0s2','FDs2'};
  rsi_ext = {'RNI','RND','RNT','FNI','HNT','HNI','HND',...
             'RI','RD','RT','FI','HT','HI','HD'};
  idx_ext = find(strcmp(rsi_int,meas));
  if ~isempty(idx_ext)
    extmeas = rsi_ext{idx_ext};
  else
    extmeas = meas;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function measlist = get_dti_export_meas()
  % the DTI measure volumes written by dti_convert_DTmeas, in the order it is
  %   asked for them in calc_dti
  measlist = {'b0','FA','V0','MD','LD','TD'};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = export_coef(parms,fit,family)
  errmsg = [];
  % export the model fit coefficients to the BIDS output as a multi-frame
  %   nii.gz, with a sidecar json describing the model and naming each frame
  % NOTE: nothing in the pipeline reads these back. The measures still come
  %   from the in-memory fit and are still written as mgz for the ROI and
  %   surface analysis; this is a copy for external consumers
  % family is the BIDS family the file is named for: dti, isdti or rsi
  if strcmp(family,'rsi')
    vol = fit.volMF;
  else
    vol = fit.volDT;
  end
  suffix = sprintf('%scoef',family);
  fname_out = sprintf('%s/%s_%s.nii.gz',parms.outdir,parms.outstem,suffix);
  fname_json = sprintf('%s/%s_%s.json',parms.outdir,parms.outstem,suffix);
  % the mask and the reference b=0 travel with the coefficients and serve both
  %   coefficient files, so they are written whether or not this one is current
  errmsg = cmig_append_errmsg(errmsg,export_coef_ref(parms,fit));
  if exist(fname_out,'file') && exist(fname_json,'file') && ~parms.forceflag
    return;
  end
  try
    cmig_write_nifti(single(vol),fname_out,fit.M);
  catch err
    tmp_errmsg = sprintf('failed to write %s: %s',fname_out,err.message);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  errmsg = cmig_append_errmsg(errmsg,write_json(coef_info(parms,fit,family),fname_json));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = export_coef_ref(parms,fit)
  errmsg = [];
  % the brain mask and the average b=0 that the coefficients are read against
  % NOTE: separate files rather than extra frames, so each stays homogeneous:
  %   the mask is binary and the b=0 is intensity. One pair serves both
  %   coefficient files, since every fit is given the same volume and mask
  % NOTE: the brain mask from get_shared_mask, not the fit's mask, which is
  %   all ones; nothing written is masked with it
  % NOTE: single rather than an integer type, even though the mask is binary.
  %   MATLAB refuses an elementwise product of an integer array with a single
  %   one, which is what dti_calc_DTmeas and rsi_calc_MFmeas do with it, so an
  %   integer mask would force every consumer to cast. A binary volume
  %   compresses to almost nothing either way
  fname_mask = sprintf('%s/%s_desc-brain_mask.nii.gz',parms.outdir,parms.outstem);
  if isempty(parms.volmask)
    fprintf('%s: WARNING: no brain mask to write to %s\n',mfilename,fname_mask);
  elseif ~exist(fname_mask,'file') || parms.forceflag
    try
      cmig_write_nifti(single(parms.volmask>0),fname_mask,fit.M);
    catch err
      tmp_errmsg = sprintf('failed to write %s: %s',fname_mask,err.message);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    end
  end
  if exist(fname_mask,'file') && (~exist(json_name(fname_mask),'file') || parms.forceflag)
    info = sidecar_info(sprintf('Brain mask from mri_synthstrip on the %s',b0_desc(parms)),...
      proc_uris(parms),dwi_ref(parms));
    info.Type = 'Brain';
    info.SkullStripped = false;
    errmsg = cmig_append_errmsg(errmsg,write_json(info,json_name(fname_mask)));
  end
  fname_ref = sprintf('%s/%s_dwiref.nii.gz',parms.outdir,parms.outstem);
  if ~exist(fname_ref,'file') || parms.forceflag
    try
      cmig_write_nifti(single(fit.volb0),fname_ref,fit.M);
    catch err
      tmp_errmsg = sprintf('failed to write %s: %s',fname_ref,err.message);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    end
  end
  if exist(fname_ref,'file') && (~exist(json_name(fname_ref),'file') || parms.forceflag)
    info = sidecar_info(sprintf('The %s, the reference image for the model fit coefficients',b0_desc(parms)),...
      proc_uris(parms),dwi_ref(parms));
    info.SkullStripped = false;
    errmsg = cmig_append_errmsg(errmsg,write_json(info,json_name(fname_ref)));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = coef_info(parms,fit,family)
  % sidecar contents for a coefficient file, in CamelCase
  info = [];
  if strcmp(family,'rsi')
    info.ModelName = 'RSI';
    info.Description = ['restriction spectrum imaging multi-compartment '...
      'fiber orientation density fit; frames are spherical harmonic '...
      'coefficients, blocked by size scale, followed by the isotropic terms'];
    info.NumberOfCoefficients = fit.nb;
    info.SHOrder = fit.SH_order;
    % NOTE: derived from the vertex count rather than hardcoded. The ico order
    %   is fixed inside create_rsi_matrix, and a change there would otherwise
    %   leave a recorded constant silently wrong while the shapes still agreed
    info.IcoOrder = get_ico_order(size(fit.icoverts,1));
    info.Lambda = fit.lambda;
    info.ADCLong = fit.ADC_long;
    info.ADCHindered = fit.ADC_hindered;
    info.ADCFree = fit.ADC_free;
    info.ADCTransVals = fit.ADC_trans_vals(:)';
    info.NumADCTrans = fit.num_ADC_trans;
    info.ADCIsoVals = fit.ADC_iso_vals(:)';
    info.NumADCIso = fit.num_ADC_iso;
    info.IsoFreeFlag = fit.iso_free_flag;
    info.FrameNames = rsi_frame_names(fit);
  else
    info.ModelName = 'DTI';
    info.Description = ['log-linear diffusion tensor fit; the first six '...
      'frames are the tensor components in the order xx, yy, zz, xy, xz, yz '...
      'and the seventh is the natural log of the fitted b=0 signal'];
    info.NumberOfCoefficients = size(fit.volDT,4);
    info.FrameNames = {'Dxx','Dyy','Dzz','Dxy','Dxz','Dyz','logS0'};
    % NOTE: bvals is already trimmed to the frames the fit was given, by
    %   max_bval and nob0_flag, so it answers both of these directly. Do not
    %   index it with i_fit, which holds indices into the untrimmed input
    info.NumberOfFramesFit = length(fit.bvals);
    info.MaxBvalFit = max(fit.bvals);
  end
  % NOTE: qmat is in RAS world (load_dwi_data), so the tensor components and
  %   the spherical harmonic coefficients are too, whatever the axis order of
  %   the grid they are written on
  info.CoefficientFrame = 'RAS';
  info.CoefficientFrameDescription = ['tensor components and spherical '...
    'harmonic coefficients are in RAS world coordinates (x right, y '...
    'anterior, z superior), not the image axes of this NIfTI grid'];
  % the same volume, mask and reference b=0 serve every fit
  info.MaskFile = sprintf('%s_desc-brain_mask.nii.gz',parms.outstem);
  info.ReferenceFile = sprintf('%s_dwiref.nii.gz',parms.outstem);
  info.Sources = proc_uris(parms);
  info.SpatialReference = dwi_ref(parms);
  % NOTE: an uncorrected file is structurally indistinguishable from a
  %   corrected one, differing only in scale, so both have to be recorded
  info.WMBiasCorrection = parms.wmbc_applied_flag;
  info.WMBiasCorrectionTarget = parms.wmbc_vol_target;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ico_order = get_ico_order(nverts)
  % ico order from a vertex count: an icosahedral sphere of order n has
  %   10*4^n + 2 vertices
  ico_order = round(log((nverts-2)/10)/log(4));
  if 10*4^ico_order + 2 ~= nverts, ico_order = -1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function names = rsi_frame_names(fit)
  % one name per coefficient frame, C<n> for the n-th compartment, following
  %   the order create_rsi_matrix builds the forward matrix in: the
  %   transverse size scales first, each a full set of spherical harmonic
  %   coefficients, then the isotropic size scales, then isotropic restricted,
  %   hindered and free, each l=0 only. With the fit here C1 is restricted,
  %   C2 hindered and C3 free
  % NOTE: within a scale the spherical harmonic column for degree l and order m
  %   is at (l^2+l+2)/2+m, even l only, which puts l=0 at frame 1, l=2 at 2:6
  %   and l=4 at 7:15 for SH_order 4
  names = {};
  c = 0;
  for t=1:fit.num_ADC_trans
    c = c + 1;
    for l=0:2:fit.SH_order
      for m=-l:l
        names{end+1} = sprintf('C%d_l%d_m%+d',c,l,m);
      end
    end
  end
  niso = fit.num_ADC_iso + fit.iso_restricted_flag + fit.iso_hindered_flag + ...
    fit.iso_free_flag;
  for t=1:niso
    c = c + 1;
    names{end+1} = sprintf('C%d_l0_m+0',c);
  end
  if length(names) ~= fit.nb
    fprintf('%s: WARNING: %d frame names for %d coefficients\n',...
      mfilename,length(names),fit.nb);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = export_meas(parms,fstem,measlist,family,rsi_measlist)
  errmsg = [];
  % export the voxelwise measure volumes to the BIDS output as nii.gz
  % NOTE: the measures otherwise exist only as mgz inside analyze_dwi, which
  %   cleanup_flag removes at the end of the run
  % NOTE: read back from the mgz just written rather than recomputed here, so
  %   the exported volume is the same array the ROI analysis reads and no
  %   scaling or frame convention is duplicated
  if ~exist('rsi_measlist','var'), rsi_measlist = []; end
  rsi = strcmp(family,'rsi');
  if rsi && isempty(rsi_measlist), rsi_measlist = parms.rsi_measlist; end
  for m=1:length(measlist)
    meas = measlist{m};
    fname_in = sprintf('%s_%s.mgz',fstem,meas);
    if ~exist(fname_in,'file')
      tmp_errmsg = sprintf('measure volume not found: %s',fname_in);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      continue;
    end
    if rsi
      extmeas = get_rsi_extmeas(meas);
    else
      extmeas = meas;
    end
    fname_out = sprintf('%s/%s_%s%s.nii.gz',...
      parms.outdir,parms.outstem,family,extmeas);
    fname_json = sprintf('%s/%s_%s%s.json',...
      parms.outdir,parms.outstem,family,extmeas);
    fname_coef = sprintf('%s/%s_%scoef.nii.gz',parms.outdir,parms.outstem,family);
    % NOTE: an RSI base measure is checked frame by frame below
    if ~rsi && exist(fname_out,'file') && exist(fname_json,'file') && ~parms.forceflag
      continue;
    end
    [vol,M] = fs_load_mgh(fname_in);
    if isempty(vol)
      tmp_errmsg = sprintf('failed to load %s',fname_in);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      continue;
    end
    frame_names = meas_frame_names(meas,extmeas,size(vol,4),rsi,rsi_measlist);
    if rsi
      % each frame of an RSI base measure is a measure of its own (RND and
      %   HND), so each is written as a single-frame file named for it
      if isempty(frame_names), frame_names = {extmeas}; end
      for k=1:length(frame_names)
        fname_k = sprintf('%s/%s_%s%s.nii.gz',...
          parms.outdir,parms.outstem,family,frame_names{k});
        if ~exist(fname_k,'file') || parms.forceflag
          try
            cmig_write_nifti(single(vol(:,:,:,k)),fname_k,M);
          catch err
            tmp_errmsg = sprintf('failed to write %s: %s',fname_k,err.message);
            fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
            errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
            continue;
          end
        elseif exist(json_name(fname_k),'file')
          continue;
        end
        info = sidecar_info(meas_description(parms,family,frame_names{k}),...
          {out_uri(parms,fname_coef)},dwi_ref(parms));
        errmsg = cmig_append_errmsg(errmsg,write_json(info,json_name(fname_k)));
      end
      continue;
    end
    try
      cmig_write_nifti(single(vol),fname_out,M);
    catch err
      tmp_errmsg = sprintf('failed to write %s: %s',fname_out,err.message);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      continue;
    end
    info = sidecar_info(meas_description(parms,family,extmeas),...
      {out_uri(parms,fname_coef)},dwi_ref(parms));
    % a multi-frame measure needs its frames named; a single-frame one does not
    if ~isempty(frame_names)
      info.FrameNames = frame_names;
    end
    errmsg = cmig_append_errmsg(errmsg,write_json(info,fname_json));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = export_meas_resT1(parms,fstem,meas,family)
  errmsg = [];
  % export the first frame of a measure resampled to the T1w grid, as
  %   _space-T1w_<family><measure>.nii.gz, the review image of
  %   cmig_manual_postqc_bids beside the AtlasTrack dseg
  % NOTE: the reference is the recon's nu.mgz, the image atlas_track registers
  %   to the atlas (fseg_resT1_xcg.mgz on its grid), and the registration
  %   M_T1_to_DTI, as for the ROIs
  fname_in = sprintf('%s_%s.mgz',fstem,meas);
  fname_T1 = sprintf('%s/%s_%s/mri/nu.mgz',...
    parms.fsurf_rootdir,parms.participant_id,parms.session_id);
  extmeas = meas;
  if strcmp(family,'rsi'), extmeas = get_rsi_extmeas(meas); end
  fname_out = sprintf('%s/%s_space-T1w_%s%s.nii.gz',...
    parms.outdir,parms.outstem,family,extmeas);
  fname_json = json_name(fname_out);
  if exist(fname_out,'file') && exist(fname_json,'file') && ~parms.forceflag
    return;
  end
  if ~exist(fname_in,'file') || ~exist(fname_T1,'file')
    tmp_errmsg = sprintf('%s or %s not found; no T1w-space %s%s',fname_in,fname_T1,family,extmeas);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  try
    [vol,M] = fs_load_mgh(fname_in);
    [M_T1,volsz_T1] = fs_read_header(fname_T1);
    vol_res = mmil_resample_vol(vol(:,:,:,1),M,...
      'M_ref',M_T1,'nvox_ref',volsz_T1(1:3),...
      'M_reg',parms.M_T1_to_DTI,'interpm',1);
    cmig_write_nifti(single(vol_res),fname_out,M_T1);
  catch err
    tmp_errmsg = sprintf('failed to write %s: %s',fname_out,err.message);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  fname_meas = sprintf('%s/%s_%s%s.nii.gz',parms.outdir,parms.outstem,family,extmeas);
  uri_T1 = fsurf_uri(parms,'mri/nu.mgz');
  info = sidecar_info(sprintf('%s, resampled to the FreeSurfer nu.mgz grid',...
    meas_description(parms,family,extmeas)),{out_uri(parms,fname_meas),uri_T1},uri_T1);
  errmsg = cmig_append_errmsg(errmsg,write_json(info,fname_json));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function names = meas_frame_names(meas,extmeas,nframes,rsi,rsi_measlist)
  % frame names for a multi-frame measure volume, empty for a single-frame one
  % NOTE: the frames of an RSI base measure are the transverse ADC size
  %   scales, and the s{N} entries in rsi_measlist name them (frame 2 of RNI
  %   is HNI); export_meas writes each RSI frame as its own file under the
  %   name given here
  names = {};
  if nframes<2, return; end
  if ~rsi
    % DTI V0 is the principal eigenvector, one frame per axis
    if strcmp(meas,'V0'), names = {'x','y','z'}; end
    return;
  end
  names = cell(1,nframes);
  names{1} = extmeas;
  for k=2:nframes, names{k} = sprintf('%ss%d',extmeas,k); end
  for k=1:length(rsi_measlist)
    tok = regexp(rsi_measlist{k},sprintf('^%ss(\\d+)$',meas),'tokens');
    if isempty(tok), continue; end
    idx = str2double(tok{1}{1});
    if idx>=1 && idx<=nframes
      names{idx} = get_rsi_extmeas(rsi_measlist{k});
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = write_json(info,fname)
  errmsg = [];
  % NOTE: PrettyPrint was added in R2021a (MATLAB 9.10); guarded so an
  %   older MATLAB writes a single line rather than failing
  if verLessThan('matlab','9.10')
    json = jsonencode(info);
  else
    json = jsonencode(info,'PrettyPrint',true);
  end
  fid = fopen(fname,'wt');
  if fid<0
    tmp_errmsg = sprintf('failed to open %s for writing',fname);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  fprintf(fid,'%s',json);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_json = json_name(fname)
  % the sidecar of an image or surface file: its whole extension replaced
  fname_json = regexprep(fname,'(\.(func|shape|label))?\.(nii\.gz|nii|gii)$','.json');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = sidecar_info(description,sources,spatial_ref)
  % the common sidecar fields of a derivative image or surface
  info = [];
  info.Description = description;
  info.Sources = sources;
  info.SpatialReference = spatial_ref;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = out_uri(parms,fname)
  % BIDS URI of a file in this module's output dir
  [~,fstem,fext] = fileparts(fname);
  uri = sprintf('bids::%s/%s/%s/%s%s',...
    parms.participant_id,parms.session_id,parms.modality,fstem,fext);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = proc_uri(parms,fname)
  % BIDS URI of a file in proc_rootdir, under the dataset name cmig_proc
  [~,fstem,fext] = fileparts(fname);
  uri = sprintf('bids:cmig_proc:%s/%s/%s/%s%s',...
    parms.participant_id,parms.session_id,parms.modality,fstem,fext);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uris = proc_uris(parms)
  % BIDS URIs of the processed dMRI runs
  uris = cellfun(@(x) proc_uri(parms,x),parms.fnames_data(:)','UniformOutput',false);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = dwi_ref(parms)
  % the processed dMRI run whose grid the outputs on the DTI grid share
  %   (load_dwi_data takes its vox2ras from the first run)
  uri = proc_uri(parms,parms.fnames_data{1});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = fsurf_uri(parms,relpath)
  % BIDS URI of a file in the recon, under the dataset name freesurfer
  uri = sprintf('bids:freesurfer:%s_%s/%s',...
    parms.participant_id,parms.session_id,relpath);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uris = surf_uris(parms,hemi)
  % the recon surfaces the measures are sampled on and resampled with
  uris = {fsurf_uri(parms,sprintf('surf/%s.white',hemi)),...
          fsurf_uri(parms,sprintf('surf/%s.sphere.reg',hemi))};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uris = at_source_uris(parms)
  % the inputs to AtlasTrack: the full-shell DTI FA and V0 and the recon nu.mgz
  uris = {out_uri(parms,sprintf('%s_dtiFA.nii.gz',parms.outstem)),...
          out_uri(parms,sprintf('%s_dtiV0.nii.gz',parms.outstem)),...
          fsurf_uri(parms,'mri/nu.mgz')};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ref,desc] = fsaverage_ref()
  % the ico7 maps: all of fsaverage's vertices
  [ref,desc] = cmig_fsaverage_ref(7);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = b0_desc(parms)
  % the average b=0 image the brain mask and the reference image come from
  if parms.wmbc_applied_flag
    desc = 'average b=0 image of the dMRI after white matter bias correction';
  else
    desc = 'average b=0 image of the dMRI';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = meas_description(parms,family,extmeas)
  % plain description of a measure, by family (dti, isdti, rsi) and its
  %   external name
  switch family
    case 'dti'
      model = 'diffusion tensor fit to all frames';
    case 'isdti'
      model = sprintf('diffusion tensor fit to the frames with b <= %g',...
        parms.inner_bval_max);
    otherwise
      model = 'restriction spectrum imaging fit';
  end
  dti_names = {'b0','FA','V0','MD','LD','TD'};
  dti_descs = {'Fitted b=0 signal','Fractional anisotropy',...
    'Principal eigenvector (x, y, z in RAS world coordinates)',...
    'Mean diffusivity','Longitudinal diffusivity','Transverse diffusivity'};
  meas_desc = extmeas;
  if strcmp(family,'rsi')
    tok = regexp(extmeas,'^([RHF])(N?)([IDT])$','tokens','once');
    if ~isempty(tok)
      comps = struct('R','Restricted','H','Hindered','F','Free');
      types = struct('I','isotropic','D','directional','T','total');
      meas_desc = comps.(tok{1});
      if ~isempty(tok{2}), meas_desc = [meas_desc ' normalized']; end
      meas_desc = sprintf('%s %s',meas_desc,types.(tok{3}));
    end
  else
    idx = find(strcmp(dti_names,extmeas),1);
    if ~isempty(idx), meas_desc = dti_descs{idx}; end
  end
  desc = sprintf('%s, from the %s',meas_desc,model);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = layer_description(parms,layer)
  % how mmil_gwcsurf_analysis samples a layer along the white surface normals
  switch layer
    case 'gm'
      desc = sprintf('sampled %g to %g mm outside the white surface (gray matter)',...
        parms.pmin,parms.pmax);
    case 'wm'
      desc = sprintf('sampled %g to %g mm inside the white surface (white matter)',...
        parms.pmin,parms.pmax);
    otherwise
      desc = 'gray/white contrast, 2*(wm-gm)/(wm+gm)';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = fwhm_description(fwhm,coverage_flag)
  if fwhm==0
    desc = ', unsmoothed';
  else
    desc = sprintf(', smoothed on the sphere to FWHM %d mm',fwhm);
  end
  if coverage_flag
    desc = [desc ', 0 outside the coverage mask'];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = strip_gii_metadata(fname)
  errmsg = [];
  % remove the UserName and Date metadata that mri_convert writes into a GIfTI
  str = fileread(fname);
  str = regexprep(str,...
    '\s*<MD>\s*<Name><!\[CDATA\[(UserName|Date)\]\]></Name>\s*<Value>.*?</Value>\s*</MD>','');
  fid = fopen(fname,'w');
  if fid<0
    tmp_errmsg = sprintf('failed to open %s for writing',fname);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  fwrite(fid,str,'char');
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = export_fseg(parms)
  errmsg = [];
  % export the fseg fiber segmentation volume to the BIDS output
  fname_fseg = sprintf('%s/fseg_resT1_xcg.mgz',parms.at_outdir);
  if ~exist(fname_fseg,'file')
    tmp_errmsg = sprintf('fseg volume not found: %s',fname_fseg);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  fname_out = sprintf('%s/%s_space-T1w_atlas-at_dseg.nii.gz',...
    parms.outdir,parms.outstem);
  fname_json = json_name(fname_out);
  if exist(fname_out,'file') && exist(fname_json,'file') && ~parms.forceflag
    return;
  end
  % NOTE: read and written here rather than through fs_mri_convert. Cast to
  %   single, the type fs_save_mgh wrote
  [vol,M] = fs_load_mgh(fname_fseg);
  if isempty(vol)
    tmp_errmsg = sprintf('failed to load %s',fname_fseg);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  % renumber AtlasTrack's labels, 10000 + the fiber number, to the codes of
  %   CMIG_ColorLUT.txt, 11000 higher, which FreeSurfer 8.x's
  %   FreeSurferColorLUT does not reuse
  ind = vol >= 10000 & vol < 12000;
  vol(ind) = vol(ind) + 11000;
  try
    cmig_write_nifti(single(vol),fname_out,M);
  catch err
    tmp_errmsg = sprintf('failed to write %s: %s',fname_out,err.message);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  uri_T1 = fsurf_uri(parms,'mri/nu.mgz');
  sources = [at_source_uris(parms),{fsurf_uri(parms,'mri/aseg.mgz')}];
  info = sidecar_info(['AtlasTrack fiber tract segmentation on the FreeSurfer '...
    'nu.mgz grid, excluding CSF and gray matter; each label is 21000 plus the '...
    'fiber number, the codes of CMIG_ColorLUT.txt'],sources,uri_T1);
  errmsg = cmig_append_errmsg(errmsg,write_json(info,fname_json));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = export_fiber_maps(parms)
  errmsg = [];
  % export the individual AtlasTrack fiber probability maps to the BIDS
  %   output as a single 4D volume at DTI resolution, with a sidecar json
  %   naming the fiber for each frame
  fiber_dir = sprintf('%s/fiber_maps',parms.at_outdir);
  if ~exist(fiber_dir,'dir')
    tmp_errmsg = sprintf('fiber maps directory not found: %s',fiber_dir);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  fname_out = sprintf('%s/%s_atlas-at_probseg.nii.gz',parms.outdir,parms.outstem);
  fname_json = sprintf('%s/%s_atlas-at_probseg.json',parms.outdir,parms.outstem);
  if exist(fname_out,'file') && exist(fname_json,'file') && ~parms.forceflag
    return;
  end
  % load fiber legend for fiber tract names
  FiberLegend = [];
  try
    FiberLegend = DTI_MMIL_Read_Fiber_Legend();
  catch err
    tmp_errmsg = sprintf('could not load fiber legend: %s',err.message);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  end
  % load each fiber probability map and concatenate as frames
  vol_fibers = []; M = eye(4);
  fiber_numbers = []; fiber_names = {};
  for f=1:length(parms.trans_fibers)
    fnum = parms.trans_fibers(f);
    fname_fiber = sprintf('%s/fiber_%d_prob_countatlas.mat',fiber_dir,fnum);
    if ~exist(fname_fiber,'file')
      tmp_errmsg = sprintf('fiber map not found: %s',fname_fiber);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      continue;
    end
    [vol,M] = mmil_load_sparse(fname_fiber);
    vol_fibers = cat(4,vol_fibers,single(vol));
    fiber_numbers(end+1) = fnum;
    fiber_names{end+1} = get_fiber_name(fnum,FiberLegend);
  end
  if isempty(vol_fibers)
    tmp_errmsg = sprintf('no fiber maps found in %s',fiber_dir);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  % save the 4D volume
  % NOTE: written straight to nii.gz rather than via an mgz and
  %   fs_mri_convert. vol_fibers is already single
  try
    cmig_write_nifti(vol_fibers,fname_out,M);
  catch err
    tmp_errmsg = sprintf('failed to write %s: %s',fname_out,err.message);
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    return;
  end
  % write sidecar json naming the fiber for each frame
  info = sidecar_info(['AtlasTrack fiber tract probabilities on the DTI grid, '...
    'one frame per fiber'],at_source_uris(parms),dwi_ref(parms));
  info.LabelMap = fiber_names;
  info.FiberNumbers = fiber_numbers;
  info.FiberNames = fiber_names;
  errmsg = cmig_append_errmsg(errmsg,write_json(info,fname_json));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fiber_name = get_fiber_name(fnum,FiberLegend)
  % fiber tract name for a fiber code, from the fiber legend
  fiber_name = sprintf('fiber_%d',fnum);
  if isempty(FiberLegend), return; end
  legend_fibernums = cell2mat({FiberLegend.FiberNumber});
  idx = find(legend_fibernums==fnum,1);
  if ~isempty(idx)
    fiber_name = FiberLegend(idx).FiberName;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = analyze_dwi_meas(parms, modality)
  errmsg = [];
  % modality: 'dti', 'dti_inner', or 'rsi'
  dti = strcmp(modality,'dti');
  dti_inner = strcmp(modality,'dti_inner');
  fprintf('%s: running %s ROI analysis...\n',mfilename,upper(modality));
  fspath = parms.fspath;
  fiber_dir = sprintf('%s/fiber_maps',parms.at_outdir);
  if dti
    meas_list = parms.dt_measlist_full;
    fstem = parms.dt_fstem_full;
  elseif dti_inner
    meas_list = parms.dt_measlist_full;
    fstem = parms.dt_fstem_inner;
  else
    meas_list = parms.rsi_measlist;
    fstem = parms.rsi_fstem;
  end
  for m = 1:length(meas_list)
    meas = meas_list{m};
    if dti || dti_inner
      scalefact = parms.dt_scalefacts_full(m);
      base_meas = meas;
      frame_num = 1;
    else
      scalefact = 1;
      tok = regexp(meas,'^(.+)s(\d+)$','tokens','once');
      if ~isempty(tok)
        base_meas = tok{1};
        frame_num = str2double(tok{2});
      else
        base_meas = meas;
        frame_num = 1;
      end
    end
    fname_meas = sprintf('%s_%s.mgz',fstem,base_meas);
    if ~exist(fname_meas,'file')
      tmp_errmsg = sprintf('%s measure file not found: %s',upper(modality),fname_meas);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      continue;
    end
    % aseg ROI extraction
    if parms.aseg_flag
      if dti
        aseg_outstem = sprintf('%s/%s_dti_%s_%s',...
          parms.tmpdir,parms.outstem,parms.dt_outfix_full,meas);
      elseif dti_inner
        aseg_outstem = sprintf('%s/%s_dti_%s_%s',...
          parms.tmpdir,parms.outstem,parms.dt_outfix_inner,meas);
      else
        aseg_outstem = sprintf('%s/%s_rsi_%s',...
          parms.tmpdir,parms.outstem,meas);
      end
      aseg_args = {...
        'outstem',aseg_outstem,...
        'scalefact',scalefact,...
        'M_reg',parms.M_T1_to_DTI,...
        'erode_flag',0,...
        'csv_flag',1,...
        'frames',frame_num,...
        'forceflag',parms.forceflag};
      % MD-dispersion-weighted averaging (production always-on): the
      % inner-shell DTI MD supplies the weighting values for all modalities
      % NOTE: the single-exponential DTI model fits multi-b data poorly,
      %   compressing the CSF-vs-tissue MD contrast the dispersion
      %   weighting relies on, so inner-shell MD is the better weight
      if parms.aseg_disp_flag
        fname_MD = sprintf('%s_MD.mgz',parms.dt_fstem_inner);
        if exist(fname_MD,'file')
          aseg_args = [aseg_args,{...
            'fname_vals',fname_MD,...
            'dispvec',parms.aseg_dispvec,...
            'disp_roicodes',parms.aseg_disp_roicodes,...
            'disp_suffix',parms.aseg_disp_suffix,...
            'dispfact',parms.aseg_dispfact}];
        else
          tmp_errmsg = sprintf('inner-shell DTI MD not found (%s); aseg dispersion weighting disabled for %s %s',fname_MD,upper(modality),meas);
          fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
          errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
        end
      end
      try
        mmil_aseg_analysis(fname_meas,fspath,aseg_args{:});
      catch err
        tmp_errmsg = sprintf('mmil_aseg_analysis failed for %s %s: %s',upper(modality),meas,err.message);
        fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
        errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      end
    end
    % fiber tract ROI extraction
    if parms.fiber_flag && exist(fiber_dir,'dir')
      if dti
        fiber_outstem = sprintf('%s/%s_dti_%s_%s_fiber',...
          parms.tmpdir,parms.outstem,parms.dt_outfix_full,meas);
      elseif dti_inner
        fiber_outstem = sprintf('%s/%s_dti_%s_%s_fiber',...
          parms.tmpdir,parms.outstem,parms.dt_outfix_inner,meas);
      else
        fiber_outstem = sprintf('%s/%s_rsi_%s_fiber',...
          parms.tmpdir,parms.outstem,meas);
      end
      try
        mmil_fiber_analysis(fname_meas,fiber_dir,...
          'outstem',fiber_outstem,...
          'fibers',parms.trans_fibers,...
          'atlas_flag',2,...
          'xcg_flag',parms.xcg_flag,...
          'scalefact',scalefact,...
          'frames',frame_num,...
          'M_reg',[],...
          'thresh_prob',parms.thresh_prob,...
          'csv_flag',1,...
          'forceflag',parms.forceflag);
      catch err
        tmp_errmsg = sprintf('mmil_fiber_analysis failed for %s %s: %s',upper(modality),meas,err.message);
        fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
        errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      end
    end
  end
  fprintf('%s: %s ROI analysis done\n',mfilename,upper(modality));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = make_coverage(parms,vol_all,M)
  % the fraction of frames in which each voxel of the dMRI series is nonzero
  %   (field of view and motion), painted onto the gm and wm surfaces and
  %   resampled to ico7 as the measures are; a vertex is covered where both
  %   are at least coverage_thresh
  % the subject-surface mask restricts the aparc ROI values; the ico7 label,
  %   intersected with fsaverage ?h.cortex.label, restricts the sphere
  %   smoothing and is written as _hemi-{L|R}_atlas-ico7_mask.func.gii (1 covered)
  errcode = 0; errmsg = [];
  parms.fnames_coverage_label = [];
  parms.fnames_coverage_mask = [];
  hemilist = {'lh','rh'};
  hemi_bids = {'L','R'};
  fname_cov = sprintf('%s/coverage.mgz',parms.tmpdir);
  outstem = sprintf('%s/coverage_%s',parms.tmpdir,parms.gwcsurf_outfix);
  try
    if ~exist(fname_cov,'file') || parms.forceflag
      fs_save_mgh(single(mean(vol_all~=0,4)),fname_cov,M);
    end
    % the projection and parameters of analyze_gwcsurf, unsmoothed
    mmil_gwcsurf_analysis(fname_cov,parms.fspath,...
      'outstem',outstem,...
      'icolevel',7,...
      'sphere_flag',parms.sphere_flag,...
      'sphsmoothsteps',0,...
      'gwmask_flag',parms.gwmask_flag,...
      'pmin',parms.pmin,...
      'pmax',parms.pmax,...
      'pstep',parms.pstep,...
      'tukey_flag',parms.tukey_flag,...
      'tukey_fact',parms.tukey_fact,...
      'M_reg',parms.M_T1_to_DTI,...
      'csv_flag',false,...
      'forceflag',parms.forceflag);
    labelpath = sprintf('%s/subjects/fsaverage/label',getenv('FREESURFER_HOME'));
    fnames_label = cell(1,2);
    fnames_mask = cell(1,2);
    for h=1:2
      hemi = hemilist{h};
      % subject surface
      gm = fs_load_mgh(sprintf('%s_gm-%s.mgz',outstem,hemi));
      wm = fs_load_mgh(sprintf('%s_wm-%s.mgz',outstem,hemi));
      covmask = min(gm,wm) >= parms.coverage_thresh;
      fnames_mask{h} = sprintf('%s/coverage_mask-%s.mgz',parms.tmpdir,hemi);
      fs_save_mgh(single(mmil_colvec(covmask)),fnames_mask{h});
      fprintf('%s: %s coverage >= %0.2f: gm %d, wm %d, both %d of %d vertices\n',...
        mfilename,hemi,parms.coverage_thresh,nnz(gm>=parms.coverage_thresh),...
        nnz(wm>=parms.coverage_thresh),nnz(covmask),numel(covmask));
      if ~parms.sphere_flag, continue; end
      % ico7
      gm = fs_load_mgh(sprintf('%s_gm-sphere-%s.mgz',outstem,hemi));
      wm = fs_load_mgh(sprintf('%s_wm-sphere-%s.mgz',outstem,hemi));
      cortex = false(numel(gm),1);
      cortex(fs_read_label(sprintf('%s/%s.cortex.label',labelpath,hemi))) = true;
      covmask = (min(gm(:),wm(:)) >= parms.coverage_thresh) & cortex;
      fnames_label{h} = sprintf('%s/coverage_ico7-%s.label',parms.tmpdir,hemi);
      fs_write_label(find(covmask),fnames_label{h},'fsaverage');
      fprintf('%s: %s ico7 coverage label: %d of %d cortex vertices\n',...
        mfilename,hemi,nnz(covmask),nnz(cortex));
      fname_out = sprintf('%s/%s_hemi-%s_atlas-ico7_mask.func.gii',...
        parms.outdir,parms.outstem,hemi_bids{h});
      if ~exist(fname_out,'file') || parms.forceflag
        fname_mgz = sprintf('%s/coverage_ico7_mask-%s.mgz',parms.tmpdir,hemi);
        fs_save_mgh(single(covmask),fname_mgz);
        fs_mri_convert(fname_mgz,fname_out);
        errmsg = cmig_append_errmsg(errmsg,strip_gii_metadata(fname_out));
        fprintf('%s: wrote %s\n',mfilename,fname_out);
      end
      if ~exist(json_name(fname_out),'file') || parms.forceflag
        [ico_ref,ico_desc] = fsaverage_ref();
        info = sidecar_info(sprintf(['dMRI coverage on %s: 1 where '...
          'the gray and white matter samples are nonzero in at least %g of the '...
          'frames, within the fsaverage cortex label, 0 elsewhere'],...
          ico_desc,parms.coverage_thresh),...
          [proc_uris(parms),surf_uris(parms,hemi)],ico_ref);
        info.Type = 'ROI';
        errmsg = cmig_append_errmsg(errmsg,write_json(info,json_name(fname_out)));
      end
    end
  catch err
    errcode = 1;
    errmsg = cmig_append_errmsg(errmsg,sprintf('coverage failed: %s',err.message));
    return;
  end
  parms.fnames_coverage_mask = fnames_mask;
  if parms.sphere_flag
    parms.fnames_coverage_label = fnames_label;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = analyze_gwcsurf(parms, modality)
  errmsg = [];
  fprintf('%s: running gwcsurf analysis (%s)...\n',mfilename,modality);
  fspath = parms.fspath;
  % NOTE: applied whatever it is, including identity, so every step uses the one
  %   matrix. Special-casing identity to [] only created a second path that
  %   behaves differently once the matrix is not identity
  M_reg_gwc = parms.M_T1_to_DTI;
  rsi = strcmp(modality,'rsi');
  dti_inner = strcmp(modality,'dti_inner');
  if rsi
    meas_list = parms.rsi_measlist;
    fstem = parms.rsi_fstem;
    gwcsurf_outfix = parms.gwcsurf_outfix;
  elseif dti_inner
    meas_list = parms.dt_measlist_full;
    fstem = parms.dt_fstem_inner;
    gwcsurf_outfix = parms.gwcsurf_outfix_inner;
  else
    meas_list = parms.dt_measlist_full;
    fstem = parms.dt_fstem_full;
    gwcsurf_outfix = parms.gwcsurf_outfix;
  end
  for m = 1:length(meas_list)
    meas = meas_list{m};
    if rsi
      scalefact = 1;
      tok = regexp(meas,'^(.+)s(\d+)$','tokens','once');
      if ~isempty(tok)
        base_meas = tok{1};
        frame_num = str2double(tok{2});
      else
        base_meas = meas;
        frame_num = 1;
      end
    else
      scalefact = parms.dt_scalefacts_full(m);
      base_meas = meas;
      frame_num = 1;
    end
    fname_meas = sprintf('%s_%s.mgz',fstem,base_meas);
    if ~exist(fname_meas,'file')
      tmp_errmsg = sprintf('%s measure file not found for gwcsurf: %s',modality,fname_meas);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      continue;
    end
    gwc_outstem = sprintf('%s/%s_%s',parms.tmpdir,meas,gwcsurf_outfix);
    for a = 1:length(parms.aparc_infix)
      aparc = parms.aparc_infix{a};
      % the annotation itself: outfix only names the output, and
      %   mmil_gwcsurf_analysis samples ?h.aparc.annot when fnames_aparc is
      %   empty, which gave the aparc.a2009s tables Desikan values
      fnames_aparc = {sprintf('%s/label/lh.%s.annot',fspath,aparc),...
                      sprintf('%s/label/rh.%s.annot',fspath,aparc)};
      try
        % sphere-resample at ico7 (production default; matches the
        % atlas-ico7 GIfTI label); concat_dwi_measures later reduces the
        % ico7 surface to ico_order (default ico5) for tabulated output
        mmil_gwcsurf_analysis(fname_meas,fspath,...
          'outstem',gwc_outstem,...
          'outfix',aparc,...
          'fnames_aparc',fnames_aparc,...
          'frames',frame_num,...
          'scalefact',scalefact,...
          'icolevel',7,...
          'sphere_flag',parms.sphere_flag,...
          'sphsmoothsteps',parms.sphsmoothsteps,...
          'gwmask_flag',parms.gwmask_flag,...
          'pmin',parms.pmin,...
          'pmax',parms.pmax,...
          'pstep',parms.pstep,...
          'tukey_flag',parms.tukey_flag,...
          'tukey_fact',parms.tukey_fact,...
          'M_reg',M_reg_gwc,...
          'fnames_coverage_label',parms.fnames_coverage_label,...
          'fnames_coverage_mask',parms.fnames_coverage_mask,...
          'forceflag',parms.forceflag);
      catch err
        tmp_errmsg = sprintf('mmil_gwcsurf_analysis failed for %s %s %s: %s',modality,meas,aparc,err.message);
        fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
        errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      end
    end
  end
  fprintf('%s: gwcsurf analysis done (%s)\n',mfilename,modality);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = convert_dwi_measures(parms, modality)
  errmsg = [];
  % convert sphere-resampled dMRI surface maps from MGZ to GIfTI format
  fprintf('%s: converting dMRI surface maps to GIfTI (%s)...\n',mfilename,modality);
  layerlist = {'gm','wm','gwc'};
  hemilist = {'lh','rh'};
  hemi_bids = {'L','R'};
  rsi = strcmp(modality,'rsi');
  dti_inner = strcmp(modality,'dti_inner');
  if rsi
    meas_list = parms.rsi_measlist;
    gwcsurf_outfix = parms.gwcsurf_outfix;
  elseif dti_inner
    meas_list = parms.dt_measlist_full;
    gwcsurf_outfix = parms.gwcsurf_outfix_inner;
  else
    meas_list = parms.dt_measlist_full;
    gwcsurf_outfix = parms.gwcsurf_outfix;
  end
  for m = 1:length(meas_list)
    meas = meas_list{m};
    if rsi
      extmeas = get_rsi_extmeas(meas);
      bids_meas = extmeas;
      bids_prefix = 'rsi';
    elseif dti_inner
      bids_meas = meas;
      bids_prefix = 'isdti';
    else
      bids_meas = meas;
      bids_prefix = 'dti';
    end
    gwc_outstem = sprintf('%s/%s_%s',parms.tmpdir,meas,gwcsurf_outfix);
    for l = 1:length(layerlist)
      layer = layerlist{l};
      bids_suffix = sprintf('%s%s%s',bids_prefix,bids_meas,layer);
      for smfi = 1:length(parms.sphsmoothsteps)
        sphsm = parms.sphsmoothsteps(smfi);
        bids_fwhm = round(1.25 * sqrt(sphsm));
        for h = 1:length(hemilist)
          hemi = hemilist{h};
          if sphsm == 0
            fname = sprintf('%s_%s-sphere-%s.mgz',gwc_outstem,layer,hemi);
          else
            fname = sprintf('%s_%s-sphere-sm%d-%s.mgz',gwc_outstem,layer,sphsm,hemi);
          end
          if ~exist(fname,'file')
            tmp_errmsg = sprintf('missing sphere file: %s',fname);
            fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
            errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
            continue;
          end
          fname_out = sprintf('%s/%s_hemi-%s_atlas-ico7_desc-fwhm%d_%s.func.gii',...
            parms.outdir,parms.outstem,...
            hemi_bids{h},bids_fwhm,bids_suffix);
          fname_json = json_name(fname_out);
          if exist(fname_out,'file') && exist(fname_json,'file') && ~parms.forceflag
            continue;
          end
          try
            % 0 outside the coverage label (the smoothed maps already are)
            if ~isempty(parms.fnames_coverage_label)
              [vals,M] = fs_load_mgh(fname);
              cov = false(numel(vals),1);
              cov(fs_read_label(parms.fnames_coverage_label{h})) = true;
              vals(~cov) = 0;
              fname = sprintf('%s_%s-sphere-sm%d-cov-%s.mgz',...
                gwc_outstem,layer,sphsm,hemi);
              fs_save_mgh(vals,fname,M);
            end
            fs_mri_convert(fname,fname_out);
            errmsg = cmig_append_errmsg(errmsg,strip_gii_metadata(fname_out));
            fprintf('%s: wrote %s\n',mfilename,fname_out);
          catch err
            tmp_errmsg = sprintf('GIfTI conversion failed for %s: %s',fname_out,err.message);
            fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
            errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
            continue;
          end
          fname_meas = sprintf('%s/%s_%s%s.nii.gz',...
            parms.outdir,parms.outstem,bids_prefix,bids_meas);
          sources = [{out_uri(parms,fname_meas)},surf_uris(parms,hemi)];
          if ~isempty(parms.fnames_coverage_label)
            sources{end+1} = out_uri(parms,sprintf('%s/%s_hemi-%s_atlas-ico7_mask.func.gii',...
              parms.outdir,parms.outstem,hemi_bids{h}));
          end
          [ico_ref,ico_desc] = fsaverage_ref();
          info = sidecar_info(sprintf('%s, %s, resampled to %s%s',...
            meas_description(parms,bids_prefix,bids_meas),layer_description(parms,layer),...
            ico_desc,fwhm_description(bids_fwhm,~isempty(parms.fnames_coverage_label))),...
            sources,ico_ref);
          errmsg = cmig_append_errmsg(errmsg,write_json(info,fname_json));
        end
      end
    end
  end
  fprintf('%s: GIfTI conversion done (%s)\n',mfilename,modality);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = concat_dwi_measures(parms, modality)
  errcode = 0; errmsg = [];
  % the concat reads the ico7 sphere maps, which gwcsurf writes only with
  %   sphere_flag
  if ~parms.sphere_flag
    fprintf('%s: skipping surface map concat (%s): sphere_flag = 0\n',mfilename,modality);
    return;
  end
  args = {'dt_measlist_full',parms.dt_measlist_full,...
          'rsi_measlist',parms.rsi_measlist,...
          'gwcsurf_outfix',parms.gwcsurf_outfix,...
          'gwcsurf_outfix_inner',parms.gwcsurf_outfix_inner,...
          'ico_order',parms.ico_order,...
          'sphsmoothsteps',parms.sphsmoothsteps,...
          'fnames_coverage_label',parms.fnames_coverage_label,...
          'forceflag',parms.forceflag};
  [~,errcode,errmsg] = cmig_concat_dwi_measures(parms.participant_id,parms.session_id,...
    parms.output_rootdir,modality,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = write_tabdat(parms)
  errmsg = [];
  fprintf('%s: writing tabulated data...\n',mfilename);
  mmil_mkdir(parms.outdir);

  % aseg dispersion-weighted roi_data files carry a suffix built the same
  % way as mmil_aseg_analysis's set_output() (e.g. '_dwtd4.7')
  if parms.aseg_disp_flag
    aseg_disp_infix = sprintf('_%s%0.1f',parms.aseg_disp_suffix,parms.aseg_dispfact);
  else
    aseg_disp_infix = '';
  end

  % fiber roi_data files carry the xcg/thresh_prob suffix from
  % mmil_fiber_analysis (e.g. '_xcg_pthresh0.08'); build it the same way
  fiber_suffix = dti_set_fiber_suffix('xcg_flag',parms.xcg_flag,...
    'thresh_prob',parms.thresh_prob);
  if isempty(fiber_suffix)
    fiber_infix = '';
  else
    fiber_infix = ['_' fiber_suffix];
  end

  % load fiber legend for fiber tract abbreviations
  FiberLegend = [];
  if parms.fiber_flag && (parms.dti_flag || parms.rsi_flag)
    try
      FiberLegend = DTI_MMIL_Read_Fiber_Legend();
    catch err
      tmp_errmsg = sprintf('could not load fiber legend: %s',err.message);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    end
  end

  % write DTI aseg TSVs
  if parms.dti_flag && parms.aseg_flag
    for m = 1:length(parms.dt_measlist_full)
      meas = parms.dt_measlist_full{m};
      fname_mat = sprintf('%s/%s_dti_%s_%s_aseg%s_roi_data.mat',...
        parms.tmpdir,parms.outstem,parms.dt_outfix_full,meas,aseg_disp_infix);
      fname_tsv = sprintf('%s/%s_atlas-aseg_dti%s.tsv',...
        parms.outdir,parms.outstem,meas);
      colstem = sprintf('mr_y_dti__fs__%s__aseg',lower(meas));
      cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
        'participant_id',parms.participant_id,'session_id',parms.session_id,...
        'colstem',colstem,'atlas','aseg','tabdir',parms.tabdir);
    end
  end

  % write RSI aseg TSVs
  if parms.rsi_flag && parms.aseg_flag
    for m = 1:length(parms.rsi_measlist)
      meas = parms.rsi_measlist{m};
      extmeas = get_rsi_extmeas(meas);
      fname_mat = sprintf('%s/%s_rsi_%s_aseg%s_roi_data.mat',...
        parms.tmpdir,parms.outstem,meas,aseg_disp_infix);
      fname_tsv = sprintf('%s/%s_atlas-aseg_rsi%s.tsv',...
        parms.outdir,parms.outstem,extmeas);
      colstem = sprintf('mr_y_rsi__%s__aseg',lower(extmeas));
      cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
        'participant_id',parms.participant_id,'session_id',parms.session_id,...
        'colstem',colstem,'atlas','aseg','tabdir',parms.tabdir);
    end
  end

  % write DTI fiber TSVs
  if parms.dti_flag && parms.fiber_flag
    for m = 1:length(parms.dt_measlist_full)
      meas = parms.dt_measlist_full{m};
      fname_mat = sprintf('%s/%s_dti_%s_%s_fiber_fiber_countatlas_wtd%s_roi_data.mat',...
        parms.tmpdir,parms.outstem,parms.dt_outfix_full,meas,fiber_infix);
      fname_tsv = sprintf('%s/%s_atlas-at_dti%s.tsv',...
        parms.outdir,parms.outstem,meas);
      colstem = sprintf('mr_y_dti__fs__%s__at',lower(meas));
      cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
        'participant_id',parms.participant_id,'session_id',parms.session_id,...
        'colstem',colstem,'atlas','at','FiberLegend',FiberLegend);
    end
  end

  % write RSI fiber TSVs
  if parms.rsi_flag && parms.fiber_flag
    for m = 1:length(parms.rsi_measlist)
      meas = parms.rsi_measlist{m};
      extmeas = get_rsi_extmeas(meas);
      fname_mat = sprintf('%s/%s_rsi_%s_fiber_fiber_countatlas_wtd%s_roi_data.mat',...
        parms.tmpdir,parms.outstem,meas,fiber_infix);
      fname_tsv = sprintf('%s/%s_atlas-at_rsi%s.tsv',...
        parms.outdir,parms.outstem,extmeas);
      colstem = sprintf('mr_y_rsi__%s__at',lower(extmeas));
      cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
        'participant_id',parms.participant_id,'session_id',parms.session_id,...
        'colstem',colstem,'atlas','at','FiberLegend',FiberLegend);
    end
  end

  % write DTI gwcsurf TSVs
  if parms.dti_flag && parms.gwcsurf_flag
    aparc_names = {'aparc','aparc.a2009s'};
    atlas_names = {'dsk','dst'};
    layerlist = {'gm','wm','gwc'};
    for m = 1:length(parms.dt_measlist_full)
      meas = parms.dt_measlist_full{m};
      gwc_outstem = sprintf('%s/%s_%s',parms.tmpdir,meas,parms.gwcsurf_outfix);
      for a = 1:length(aparc_names)
        aparc = aparc_names{a};
        atlas = atlas_names{a};
        for l = 1:length(layerlist)
          layer = layerlist{l};
          fname_mat = sprintf('%s_%s_%s_roi_data.mat',gwc_outstem,layer,aparc);
          fname_tsv = sprintf('%s/%s_atlas-%s_dti%s%s.tsv',...
            parms.outdir,parms.outstem,atlas,meas,layer);
          colstem = sprintf('mr_y_dti__fs__%s__%s__%s',lower(meas),layer,atlas);
          cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
            'participant_id',parms.participant_id,'session_id',parms.session_id,...
            'colstem',colstem,'atlas',atlas,'tabdir',parms.tabdir,...
            'fspath',parms.fspath,'aparc_infix',aparc);
        end
      end
    end
  end

  % write RSI gwcsurf TSVs
  if parms.rsi_flag && parms.gwcsurf_flag
    aparc_names = {'aparc','aparc.a2009s'};
    atlas_names = {'dsk','dst'};
    layerlist = {'gm','wm','gwc'};
    for m = 1:length(parms.rsi_measlist)
      meas = parms.rsi_measlist{m};
      extmeas = get_rsi_extmeas(meas);
      gwc_outstem = sprintf('%s/%s_%s',parms.tmpdir,meas,parms.gwcsurf_outfix);
      for a = 1:length(aparc_names)
        aparc = aparc_names{a};
        atlas = atlas_names{a};
        for l = 1:length(layerlist)
          layer = layerlist{l};
          fname_mat = sprintf('%s_%s_%s_roi_data.mat',gwc_outstem,layer,aparc);
          fname_tsv = sprintf('%s/%s_atlas-%s_rsi%s%s.tsv',...
            parms.outdir,parms.outstem,atlas,extmeas,layer);
          colstem = sprintf('mr_y_rsi__%s__%s__%s',lower(extmeas),layer,atlas);
          cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
            'participant_id',parms.participant_id,'session_id',parms.session_id,...
            'colstem',colstem,'atlas',atlas,'tabdir',parms.tabdir,...
            'fspath',parms.fspath,'aparc_infix',aparc);
        end
      end
    end
  end

  % write inner-shell DTI aseg TSVs
  if parms.inner_dti_flag && parms.aseg_flag
    for m = 1:length(parms.dt_measlist_full)
      meas = parms.dt_measlist_full{m};
      fname_mat = sprintf('%s/%s_dti_%s_%s_aseg%s_roi_data.mat',...
        parms.tmpdir,parms.outstem,parms.dt_outfix_inner,meas,aseg_disp_infix);
      fname_tsv = sprintf('%s/%s_atlas-aseg_isdti%s.tsv',...
        parms.outdir,parms.outstem,meas);
      colstem = sprintf('mr_y_dti__is__%s__aseg',lower(meas));
      cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
        'participant_id',parms.participant_id,'session_id',parms.session_id,...
        'colstem',colstem,'atlas','aseg','tabdir',parms.tabdir);
    end
  end

  % write inner-shell DTI fiber TSVs
  if parms.inner_dti_flag && parms.fiber_flag
    for m = 1:length(parms.dt_measlist_full)
      meas = parms.dt_measlist_full{m};
      fname_mat = sprintf('%s/%s_dti_%s_%s_fiber_fiber_countatlas_wtd%s_roi_data.mat',...
        parms.tmpdir,parms.outstem,parms.dt_outfix_inner,meas,fiber_infix);
      fname_tsv = sprintf('%s/%s_atlas-at_isdti%s.tsv',...
        parms.outdir,parms.outstem,meas);
      colstem = sprintf('mr_y_dti__is__%s__at',lower(meas));
      cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
        'participant_id',parms.participant_id,'session_id',parms.session_id,...
        'colstem',colstem,'atlas','at','FiberLegend',FiberLegend);
    end
  end

  % write inner-shell DTI gwcsurf TSVs
  if parms.inner_dti_flag && parms.gwcsurf_flag
    aparc_names = {'aparc','aparc.a2009s'};
    atlas_names = {'dsk','dst'};
    layerlist = {'gm','wm','gwc'};
    for m = 1:length(parms.dt_measlist_full)
      meas = parms.dt_measlist_full{m};
      gwc_outstem = sprintf('%s/%s_%s',parms.tmpdir,meas,parms.gwcsurf_outfix_inner);
      for a = 1:length(aparc_names)
        aparc = aparc_names{a};
        atlas = atlas_names{a};
        for l = 1:length(layerlist)
          layer = layerlist{l};
          fname_mat = sprintf('%s_%s_%s_roi_data.mat',gwc_outstem,layer,aparc);
          fname_tsv = sprintf('%s/%s_atlas-%s_isdti%s%s.tsv',...
            parms.outdir,parms.outstem,atlas,meas,layer);
          colstem = sprintf('mr_y_dti__is__%s__%s__%s',lower(meas),layer,atlas);
          cmig_write_tabdat_dmri('fname_in',fname_mat,'fname_out',fname_tsv,...
            'participant_id',parms.participant_id,'session_id',parms.session_id,...
            'colstem',colstem,'atlas',atlas,'tabdir',parms.tabdir,...
            'fspath',parms.fspath,'aparc_infix',aparc);
        end
      end
    end
  end

  fprintf('%s: tabulated data done\n',mfilename);
return;

