function auto_qcinfo = autoqc_dMRI_shells(indir,varargin)
%function auto_qcinfo = autoqc_dMRI_shells(indir,[options])
%
% Purpose: Calculate automated quality control metrics from dMRI data for each b-value shell
%
% Usage:
%  autoqc_dMRI_shells(indir,'key1', value1,...);
%
% Required Input:
%  indir: full path of DTIPROC directory
%
% Optional Paramters:
%  'data': struct containing fields vol, M, qmat, bvals
%    if supplied, indir may be empty
%    optional fields vol_raw, qmat_raw, bvals_raw: raw (uncorrected) data
%      if supplied, slice censoring is detected from a fit to the raw
%      data; required when vol holds corrected data, whose slices are no
%      longer acquisition-aligned and whose censored slices are already
%      replaced; qmat_raw/bvals_raw default to qmat/bvals
%    {default = []}
%  'outdir': output directory
%    provide full path or relative to indir
%    if empty, will write output to indir
%    {default = []}
%  'outstem': output file stem
%    {default = []}
%  'tmpdir': temporary processing file output directory
%    full path or relative to indir
%    {default = 'autoQC'}
%  'infix': file name infix (e.g. [], 'corr_regT1')
%    {default = []}
%  'revflag': [0|1|2] specify whether to use non-rev or rev data
%    if revflag=0, use non-rev data
%    if revflag=1, use rev data
%      rev scans have names like 'DTI1_rev.mgz'
%    if revflag=2, use concatenated non-rev and rev data
%    {default = 0}
%  'min_ndirs': require at least this many diffusion directions
%    {default = 10}
%  'min_bval': minimum b value a scan must have to be included
%    {default = 2000}
%  'flex_flag': [0|1] DTI_flex scans included in tensor fit
%    {default = 1}
%  'meas_flag': [0|1] calculate DTI measures from tensor fit
%    {default = 0}
%  'fname_DTmeas': output mat file for DT measures (requires meas_flag)
%    for conversion and image plots by the caller
%    {default = []}
%  'fname_DTerr': output mat file for the DT rms error volume
%    for conversion and image plots by the caller
%    {default = []}
%  'cleanupflag': [0|1] remove temporary files
%    {default = 1}
%  'verbose': [0|1] display status messages and warnings
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Output:
%   qcinfo: struct containing metrics related to DTI fit error
%
% Created:  11/08/21 by Don Hagler
% Prev Mod: 01/30/24 by Don Hagler
% Prev Mod: 04/23/25 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 06/19/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 08/12/2026 by Don Hagler
% Prev Mod: 08/20/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: calculate averages for white matter, gray matter, ventricles, brain
%%       FreeSurfer as optional input?
%%       or use corseg or FSL FIRST so can be used on raw?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on MMIL_autoQC_Exam
% Created:  02/20/16 by Don Hagler
% Prev Mod: 08/21/19 by Feng Xue
% Prev Mod: 02/24/21 by Don Hagler
% Prev Mod: 08/21/2026 by Don Hagler
% Prev Mod: 08/22/2026 by Don Hagler
% Last Mod: 08/25/2026 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

auto_qcinfo = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(indir,varargin);

% load data
[data,parms] = load_data(parms);

% initialize output values
auto_qcinfo = init_output(data);

% calculate tSNR metrics  
auto_qcinfo = calc_tSNR_metrics(data,auto_qcinfo,parms);
    
% calculate FWHM metrics
auto_qcinfo = calc_FWHM_metrics(data,auto_qcinfo,parms);

% calculate DTI metrics
auto_qcinfo = calc_DTI_metrics(data,auto_qcinfo,parms);

% calculate RSI metrics
if parms.flex_flag
  auto_qcinfo = calc_RSI_metrics(data,auto_qcinfo,parms);
end

if parms.cleanupflag
  cleanup_tmpdir(parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(indir,options)
  parms = mmil_args2parms(options,{...
    'indir',indir,[],...
    'data',[],[],...
  ...
    'outdir',[],[],...
    'outstem',[],[],...
    'tmpdir','autoQC',[],...
    'infix',[],[],...
    'revflag',0,[0,1,2],...
    'min_ndirs',10,[2,1000],...
    'min_bval',2000,[],...
    'flex_flag',true,[false true],...
    'meas_flag',false,[false true],...
    'fname_DTmeas',[],[],...
    'fname_DTerr',[],[],...
    'cleanupflag',true,[false true],...
    'verbose',false,[false true],...
    'bids_flag',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden for censoring
    'censor_thresh',[],[],...
    'censor_prctile',99,[],...
    'min_fract',0.05,[],...
    'lambda',0.01,[],... % regularization term for RSI fit
  ... % hidden for subplots
    'planestrings',{'HOR','SAG','COR'},[],...
    'threeviews_plim',[1,99],[],...
    'threeviews_slice_fract',[0.5,0.6,0.5],[],...
    'threeviews_fig_size',[10 4],[],...
    'censor_fig_size',[4 4],[],...
    'fontsize',8,[],...
    'tif_dpi',300,[],...
    'normmax',[],[],...
    'relmax',1,[],...
    'colorbar_flag',true,[false true],...
    'plot_extra_flag',false,[false true],...
  ...
    'DTI_meas_list',{'FA','MD','b0'},[],...
    'RSI_meas_list',{'N0','ND','N0s2','NDs2','NIf'},[],...
  ...
    'load_data_tags',{'infix','revflag','min_ndirs','min_bval','flex_flag'},[],...
    'mri_convert_tags',{'out_orient','options','verbose','forceflag'},[],...
    'brainmask_tags',{'outdir','fname_mask','binarize_flag','verbose','forceflag'},[],...
    'info_tags',{'norm_bvals_flag','norm_slices_flag','winsorize_flag',...
                 'clip_neg_flag','erode_nvoxels','smooth_nvoxels',....
                 'censor_thresh','censor_prctile','min_fract',...
                 'sigma','thresh'},[],...
  });

  parms.outfix = parms.infix;
  if ~isempty(parms.outfix), parms.outfix = ['_' parms.outfix]; end

  % NOTE: the intermediate volumes are named from outstem, so two runs of one
  %   session, or a rerun of either, cannot read each other's files back. Their
  %   FWHM outputs inherit the name through cmig_3dFWHM, whose outstem defaults to
  %   the input stem. Falls back to 'data' for callers that supply no outstem
  parms.tmpstem = 'data';
  if ~isempty(parms.outstem), parms.tmpstem = parms.outstem; end

  if isempty(parms.outdir)
    parms.outdir = parms.indir;
  elseif mmil_isrelative(parms.outdir)
    parms.outdir = [parms.indir '/' parms.outdir];
  end;
  
  if isempty(parms.tmpdir)
    parms.tmpdir = parms.indir;
  elseif mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [parms.indir '/' parms.tmpdir];
  end;
  
  mmil_mkdir(parms.outdir);
  mmil_mkdir(parms.tmpdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = init_output(data)
  auto_qcinfo = [];
  %% todo: metrics for wm, gm, and ventricles
  auto_qcinfo.brain_FWHM = [];
  auto_qcinfo.brain_FWHM_resid = [];
  auto_qcinfo.brain_tSNR = [];
  auto_qcinfo.brain_DTmeas = [];
  auto_qcinfo.brain_DTerr = [];
  auto_qcinfo.brain_RSImeas = [];
  auto_qcinfo.brain_RSIerr = [];
  auto_qcinfo.bvals = data.bvals;
  auto_qcinfo.round_bvals = data.round_bvals;
  auto_qcinfo.unique_bvals = data.unique_bvals;
  auto_qcinfo.nbvals = data.nbvals;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,parms] = load_data(parms)
  data = [];
  
  % check data
  if ~isempty(parms.data)
    data = parms.data;
    required_fields = {'vol','M','qmat','bvals'};
    for r=1:length(required_fields)
      if ~isfield(data,required_fields{r})
        error('data struct missing field %s',required_fields{r});
      end
    end
  else
    if parms.verbose
      fprintf('%s: loading dMRI data...\n',mfilename);
    end;
    % load dMRI data
    a = mmil_parms2args(parms,parms.load_data_tags);
    [data.vol,data.M,data.qmat,data.bvals,TEs,errcode]=DTI_MMIL_Load_Data(parms.indir,a{:});
    if errcode, return; end
  end
  
  % identify unique b-values
  [data.round_bvals,data.unique_bvals,data.nbvals] = get_unique_bvals(data.bvals);
  fprintf('%s: %d unique b-values: %s\n',...
    mfilename,length(data.unique_bvals),sprintf('%d ',data.unique_bvals));
  % check mask
  if ~isfield(data,'vol_mask')
    % create brain mask
    [data.vol_mask,parms.fname_mask] = create_brainmask(data.vol,data.M,data.round_bvals,parms);
  else
    % save mask file
    % NOTE: nifti, written only so cmig_3dFWHM can hand it to 3dFWHMx.
    %   Nothing reads it back
    parms.fname_mask = sprintf('%s/%s%s_mask.nii',parms.tmpdir,parms.tmpstem,parms.outfix);
    cmig_write_nifti(single(data.vol_mask),parms.fname_mask,data.M);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = calc_DTI_metrics(data,auto_qcinfo,parms)
  % calculate diffusion tensor fit
  DTfit = calc_DTfit(data,parms);
  if isempty(DTfit), return; end
  % override mask
  DTfit.volmask = data.vol_mask;
  % calculate residual error from diffusion tensor fit
  DTerr = calc_DTerr(DTfit,data,parms);
  % save rms error volume for conversion and image plots by caller
  if ~isempty(parms.fname_DTerr)
    save_DTerr(DTerr,data,parms);
  end
  % calculate average DTI err within brain mask
  %% todo: gm, wm, and ventricles masks
  auto_qcinfo.brain_DTerr = ...
    calc_brain_DTerr(DTerr,data,parms);
  % plot images of residual error
  if parms.bids_flag
    plot_err_images(DTerr,data,parms,parms.outstem,'DTI');
    if isempty(parms.outstem)
      fname_censor = sprintf('%s/DTI_censor.jpg',parms.outdir);
    else
      fname_censor = sprintf('%s/%s_desc-DTI_censor.jpg',parms.outdir,parms.outstem);
    end
  else
    if isempty(parms.outstem)
      outstem = 'DTI';
    else
      outstem = [parms.outstem '_DTI'];
    end
    plot_err_images(DTerr,data,parms,outstem);
    fname_censor = sprintf('%s/%s_censor.jpg',parms.outdir,outstem);
  end
  % calculate measures from diffusion tensor fit (FA, MD, b0)
  if parms.meas_flag
    DTmeas = calc_DTmeas(DTfit,parms);
    % save DT measures for conversion and image plots by caller
    if ~isempty(parms.fname_DTmeas)
      save_DTmeas(DTmeas,data,parms);
    end
    % calculate average DTI measures within brain mask
    %% todo: gm, wm, and ventricles masks
    auto_qcinfo.brain_DTmeas = ...
      calc_brain_DTmeas(DTmeas,data,parms);
  end
  % calculate censoring, from raw data if supplied
  [cen_fit,cen_err,cen_data] = set_censor_input(DTfit,DTerr,data,parms,'DTI');
  [auto_qcinfo.DTI_censor_info,censor_data] = calc_censor_info(cen_fit,cen_err,cen_data,parms);
  % spatial smoothness of the fit residual, on the same data the censoring
  % detection used, so raw when the caller supplied it. This is the only
  % measure of noise smoothness alone: the data arm includes anatomy, which
  % contributes roughly a factor of two, and unlike the data arm it does not
  % vary with subject anatomy
  if ~isempty(cen_err) && isfield(cen_err,'vol') && ~isempty(cen_err.vol)
    infix = '_dterr';
    if isfield(data,'vol_raw') && ~isempty(data.vol_raw), infix = '_raw_dterr'; end
    auto_qcinfo.brain_FWHM_resid = calc_shell_fwhm(cen_err.vol,...
      cen_data.round_bvals,cen_data.unique_bvals,cen_data.M,infix,parms);
  end
  % plot images of censor data
  plot_censor_images(auto_qcinfo.DTI_censor_info,censor_data,parms,fname_censor);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = calc_RSI_metrics(data,auto_qcinfo,parms)
  % calculate diffusion tensor fit
  RSIfit = calc_RSIfit(data,parms);
  if isempty(RSIfit), return; end
  % override mask
  RSIfit.volmask = data.vol_mask;
  % calculate residual error from diffusion tensor fit
  RSIerr = calc_RSIerr(RSIfit,data,parms);
  % calculate average DTI err within brain mask
  %% todo: gm, wm, and ventricles masks
  auto_qcinfo.brain_RSIerr = ...
    calc_brain_RSIerr(RSIerr,data,parms);
  % plot images of residual error
  if parms.bids_flag
    plot_err_images(RSIerr,data,parms,parms.outstem,'RSI');
    if isempty(parms.outstem)
      fname_censor = sprintf('%s/RSI_censor.jpg',parms.outdir);
    else
      fname_censor = sprintf('%s/%s_desc-RSI_censor.jpg',parms.outdir,parms.outstem);
    end
  else
    if isempty(parms.outstem)
      outstem = 'RSI';
    else
      outstem = [parms.outstem '_RSI'];
    end
    plot_err_images(RSIerr,data,parms,outstem);
    fname_censor = sprintf('%s/%s_censor.jpg',parms.outdir,outstem);
  end
  % calculate measures from RSI fit
  if parms.meas_flag
    RSImeas = calc_RSImeas(RSIfit,parms);
    % calculate average RSI measures within brain mask
    %% todo: gm, wm, and ventricles masks
    auto_qcinfo.brain_RSImeas = ...
      calc_brain_RSImeas(RSImeas,data,parms);
  end
  % calculate censoring, from raw data if supplied
  [cen_fit,cen_err,cen_data] = set_censor_input(RSIfit,RSIerr,data,parms,'RSI');
  [auto_qcinfo.RSI_censor_info,censor_data] = calc_censor_info(cen_fit,cen_err,cen_data,parms);
  % plot images of censor data
  plot_censor_images(auto_qcinfo.RSI_censor_info,censor_data,parms,fname_censor);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = calc_tSNR_metrics(data,auto_qcinfo,parms)
  % calculate median tSNR within brain mask
  median_brain_tSNR_bvals = zeros(data.nbvals,1);
  for b=1:length(data.unique_bvals)
    bval = data.unique_bvals(b);
    idx = find(data.round_bvals==bval);
    if parms.verbose
      fprintf('%s: calculating tSNR for %d frames at b=%d...\n',...
        mfilename,length(idx),bval);
    end
    vol_tmp = data.vol(:,:,:,idx);
    vol_mean = mean(vol_tmp,4);
    vol_std = std(vol_tmp,[],4);
    vol_snr = vol_mean ./ (vol_std+eps);
    tSNR = median(vol_snr(data.vol_mask>0));
    if parms.verbose
      fprintf('%s: median tSNR = %0.2f\n',mfilename,tSNR);
    end
    median_brain_tSNR_bvals(b) = tSNR;
    % plot tSNR images
    vol_tmp = vol_snr .* data.vol_mask;
    if parms.bids_flag
      if isempty(parms.outstem)
        outstem = sprintf('tsnrB%d',bval);
      else
        outstem = sprintf('%s_desc-b%d_tsnr',parms.outstem,bval);
      end
    else
      if isempty(parms.outstem)
        outstem = sprintf('tSNR_b%d',bval);
      else
        outstem = sprintf('%s_tSNR_b%d',parms.outstem,bval);
      end
    end
    plot_threeviews(vol_tmp,data.M,outstem,parms);
  end
  auto_qcinfo.brain_tSNR.median_brain_tSNR_bvals = median_brain_tSNR_bvals;
return;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = calc_FWHM_metrics(data,auto_qcinfo,parms)
  % measure the RAW data when the caller supplied it. With detrending off the
  % estimate is per-frame, so eddy current and motion correction cannot remove
  % blur from a frame, only add interpolation smoothing - measured at about
  % 1.5x on 10 scans. Callers with no raw data, e.g. post-processing QC, get
  % the volume they passed, which is what they want to measure
  infix = '';
  if isfield(data,'vol_raw') && ~isempty(data.vol_raw)
    data = set_raw_data(data);
    infix = '_raw';
  end
  auto_qcinfo.brain_FWHM = calc_shell_fwhm(data.vol,data.round_bvals,...
    data.unique_bvals,data.M,infix,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = calc_shell_fwhm(vol,round_bvals,unique_bvals,M,infix,parms)
  % spatial smoothness per b-value shell, plus a frame-weighted overall value
  res = [];
  p = [];
  p.outdir = parms.tmpdir;
  p.fname_mask = parms.fname_mask;
  % no detrending: frames within a shell differ by diffusion direction, with
  % no ordering for a polynomial to fit, and AFNI's implicit order
  % (nframes/30) would leave every shell but b=3000 undetrended, putting the
  % shells in two incomparable regimes
  p.detrend_flag = 0;
  p.verbose = 0;
  % NOTE: parms.forceflag, not 0. cmig_3dFWHM raises its runflag only when the
  %   output file is missing or zero length, so hardcoding 0 read an existing
  %   _fwhm.txt back verbatim while the volumes beside it were being recomputed
  p.forceflag = parms.forceflag;
  nbvals = length(unique_bvals);
  brain_FWHM = nan(nbvals,3);
  nframes_bvals = zeros(nbvals,1);
  for b=1:nbvals
    bval = unique_bvals(b);
    idx = find(round_bvals==bval);
    nframes_bvals(b) = length(idx);
    if parms.verbose
      fprintf('%s: calculating FWHM%s for %d frames at b=%d...\n',...
        mfilename,infix,length(idx),bval);
    end
    % NOTE: nifti, since 3dFWHMx reads neither mgz nor mgh; uncompressed
    %   because cleanup_flag deletes it. Cast to single: the dterr volumes are
    %   double, and these are multi-frame, so it halves a large file
    fname_tmp = sprintf('%s/%s%s%s_b%d.nii',...
      parms.tmpdir,parms.tmpstem,parms.outfix,infix,bval);
    if ~exist(fname_tmp,'file') || parms.forceflag
      cmig_write_nifti(single(vol(:,:,:,idx)),fname_tmp,M);
    end
    a = mmil_parms2args(p);
    fwhm = cmig_3dFWHM(fname_tmp,a{:});
    if parms.verbose
      fprintf('%s: FWHM%s = %0.2f, %0.2f, %0.2f\n',...
        mfilename,infix,fwhm(1),fwhm(2),fwhm(3));
    end
    brain_FWHM(b,:) = fwhm;
  end
  res.mean_brain_FWHMx_bvals = brain_FWHM(:,1);
  res.mean_brain_FWHMy_bvals = brain_FWHM(:,2);
  res.mean_brain_FWHMz_bvals = brain_FWHM(:,3);
  res.nframes_bvals = nframes_bvals;
  % overall value as a frame-weighted mean of the per-shell values. Exact,
  % since each per-shell value is itself an average over that shell's frames,
  % so a frame-weighted mean of shell means equals the mean over all frames
  axes_names = {'x','y','z'};
  for k=1:3
    v = brain_FWHM(:,k);
    ind = find(~isnan(v) & nframes_bvals>0);
    if isempty(ind)
      val = NaN;
    else
      val = sum(v(ind).*nframes_bvals(ind))/sum(nframes_bvals(ind));
    end
    res.(sprintf('mean_brain_FWHM%s',axes_names{k})) = val;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [round_bvals,unique_bvals,nbvals] = get_unique_bvals(bvals)
  round_bvals = 50*round(bvals/50);
  unique_bvals = unique(round_bvals);
  nbvals = length(unique_bvals);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol_mask,fname_mask] = create_brainmask(vol,M,bvals,parms);
  vol_mask = []; fname_mask = [];
  % NOTE: nifti so this function does not write two formats.
  %   mmil_synthstrip_brainmask picks its reader from the extension, so the
  %   mask comes back nifti and is read as such below
  fname_b0 = sprintf('%s/%s%s_b0.nii',parms.tmpdir,parms.tmpstem,parms.outfix);
  if ~exist(fname_b0,'file') || parms.forceflag
    idx = find(bvals==0);
    vol_tmp = vol(:,:,:,idx);
    cmig_write_nifti(single(vol_tmp),fname_b0,M);
  end
  p = parms;
  p.outdir = parms.tmpdir;
  p.binarize_flag = 1;
  p.cleanup_flag = 1;
  a = mmil_parms2args(p,parms.brainmask_tags);
  fname_mask = cmig_create_brainmask(fname_b0,a{:});
  vol_mask = cmig_read_nifti(fname_mask);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTfit = calc_DTfit(data,parms)
  nframes = size(data.qmat,1);
  if parms.verbose
    fprintf('%s: running tensor fit calculations with %d frames...\n',...
      mfilename,nframes);
    tic
  end
  DTfit = dti_fit_tensor(data.vol,data.qmat,'bvals',data.bvals,'M',data.M);
  if parms.verbose, toc; end
  if isempty(DTfit)
    fprintf('%s: WARNING: tensor fit calculations failed\n',mfilename);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTmeas = calc_DTmeas(DTfit,parms)
  DTmeas = [];
  if parms.verbose
    fprintf('%s: calculating DT measures...\n',mfilename);
    tic
  end
  % calculate Eigenvectors, FA, etc.
  DTmeas = dti_calc_DTmeas(DTfit,0);
  if parms.verbose, toc; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DTerr = calc_DTerr(DTfit,data,parms)
  DTerr = [];
  if parms.verbose
    fprintf('%s: calculating DT error...\n',mfilename);
    tic
  end

  % calculate synthesized volume
  vol_synth = dti_synth_vol(DTfit);

  DTerr.unique_bvals = data.unique_bvals;
  DTerr.nbvals = data.nbvals;
  DTerr.M = data.M;

  % calculate rms residual error overall and for each b-value
  nslices = DTfit.volsz(3);
  volsz = [DTfit.volsz(1:3)];
  vol_err = zeros(DTfit.volsz);
  vol_rms = zeros(volsz);
  vol_rms_err = zeros(volsz);
  vol_rms_err_norm = zeros(volsz);
  vol_rms_err_rel = zeros(volsz);
  volsz_bvals = [volsz,data.nbvals];
  vol_rms_bvals = zeros(volsz_bvals);
  vol_rms_err_bvals = zeros(volsz_bvals);
  vol_rms_err_norm_bvals = zeros(volsz_bvals);
  vol_rms_err_rel_bvals = zeros(volsz_bvals);
  
  % calculate RMS signal and error across frames
  % one slice at a time to save memory
  for z = 1:nslices
    % calculate RMS error across all frames
    slc = data.vol(:,:,z,:);
    slc_synth = vol_synth(:,:,z,:);
    vol_err(:,:,z,:) = slc_synth-slc;
    vol_rms(:,:,z) = sqrt(mean(slc.^2,4));
    vol_rms_err(:,:,z) = sqrt(mean((slc_synth-slc).^2,4));
    % calculate RMS error across frames for each b-value
    for b=1:data.nbvals
      bval = data.unique_bvals(b);
      idx = find(data.round_bvals==bval);
      slc = data.vol(:,:,z,idx);
      slc_synth = vol_synth(:,:,z,idx);
      vol_rms_bvals(:,:,z,b) = sqrt(mean(slc.^2,4));
      vol_rms_err_bvals(:,:,z,b) = sqrt(mean((slc_synth-slc).^2,4));
    end
  end

  % calculate median rms value of b=0 frames: the common normalizer for the
  %   overall error and for every b-value, so the norm measures share one scale.
  %   A median across all frames would depend on how many frames are at each
  %   b-value, so the same brain under a different shell distribution would
  %   normalize differently
  vol_tmp = vol_rms_bvals(:,:,:,1);
  med_rms_b0 = median(mmil_rowvec(vol_tmp(data.vol_mask>0)));
  % normalize err by the b=0 median rms
  vol_rms_err_norm = vol_rms_err / med_rms_b0;
  % normalize err by rms for each voxel
  vol_rms_err_rel = min(1,vol_rms_err ./ (vol_rms + eps));
  % normalize error for each b-value
  for b=1:data.nbvals
    vol_rms_err_norm_bvals(:,:,:,b) = vol_rms_err_bvals(:,:,:,b) / med_rms_b0;
    vol_rms_err_rel_bvals(:,:,:,b) = min(1,vol_rms_err_bvals(:,:,:,b) ./ (vol_rms_bvals(:,:,:,b) + eps));
  end
      
  % save output volumes
  DTerr.vol = vol_err;
  DTerr.vol_rms_err = vol_rms_err;
  DTerr.vol_rms_err_norm = vol_rms_err_norm;
  DTerr.vol_rms_err_rel = vol_rms_err_rel;
  DTerr.vol_rms_err_norm_bvals = vol_rms_err_norm_bvals;
  DTerr.vol_rms_err_rel_bvals = vol_rms_err_rel_bvals;

  if parms.verbose, toc; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_DTmeas(DTmeas,data,parms)
  % save DT measures for conversion and image plots by caller
  if exist(parms.fname_DTmeas,'file') && ~parms.forceflag, return; end
  M = data.M;
  save(parms.fname_DTmeas,'DTmeas','M');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_DTerr(DTerr,data,parms)
  % save the rms error volume for conversion and image plots by caller
  if exist(parms.fname_DTerr,'file') && ~parms.forceflag, return; end
  M = data.M;
  DTerr = struct('vol_rms_err',DTerr.vol_rms_err);
  save(parms.fname_DTerr,'DTerr','M');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_DTerr = calc_brain_DTerr(DTerr,data,parms)
  brain_DTerr = [];
  if parms.verbose
    fprintf('%s: calculating DTI residual error within brain...\n',mfilename);
  end

  % get brain voxels from brainmask
  ind_brain = find(data.vol_mask>0);

  % calculate per-frame residual error within brain, relative to the data
  %   the fit used; the whole-brain rawqc measures are the mean, median and
  %   standard deviation of that ratio across frames
  nframes = size(data.vol,4);
  brain_rms_vals = zeros(1,nframes);
  brain_rms_err  = zeros(1,nframes);
  for f=1:nframes
    vol_tmp = data.vol(:,:,:,f);
    brain_rms_vals(f) = sqrt(mean(vol_tmp(ind_brain).^2));
    vol_tmp = DTerr.vol(:,:,:,f);
    brain_rms_err(f)  = sqrt(mean(vol_tmp(ind_brain).^2));
  end
  brain_DTerr.brain_rms_err  = brain_rms_err;
  brain_DTerr.brain_rms_vals = brain_rms_vals;
  brain_DTerr.brain_rms_err_ratio = brain_rms_err ./ (brain_rms_vals + eps);
  brain_DTerr.brain_rms_err_ratio_mean   = mean(brain_DTerr.brain_rms_err_ratio);
  brain_DTerr.brain_rms_err_ratio_median = median(brain_DTerr.brain_rms_err_ratio);
  brain_DTerr.brain_rms_err_ratio_std    = std(brain_DTerr.brain_rms_err_ratio);

  % calculate median residual error within brain
  brain_DTerr.median_brain_rms_err_norm = median(DTerr.vol_rms_err_norm(ind_brain));
  brain_DTerr.median_brain_rms_err_rel = median(DTerr.vol_rms_err_rel(ind_brain));
  
  % calculate median residual error within brain for each b-value
  brain_DTerr.median_brain_rms_err_norm_bvals = zeros(1,DTerr.nbvals);
  brain_DTerr.median_brain_rms_err_rel_bvals = zeros(1,DTerr.nbvals);
  for b=1:DTerr.nbvals
    tmp_err = DTerr.vol_rms_err_norm_bvals(:,:,:,b);
    brain_DTerr.median_brain_rms_err_norm_bvals(b) = median(tmp_err(ind_brain));
    tmp_err = DTerr.vol_rms_err_rel_bvals(:,:,:,b);
    brain_DTerr.median_brain_rms_err_rel_bvals(b) = median(tmp_err(ind_brain));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: calculate metrics within wm, gm, and ventricles
function brain_DTmeas = calc_brain_DTmeas(DTmeas,data,parms)
  brain_DTmeas = [];
  if parms.verbose
    fprintf('%s: calculating DTI measures within brain...\n',mfilename);
  end;
  ind_brain = find(data.vol_mask>0);
  for m=1:length(parms.DTI_meas_list)
    meas = parms.DTI_meas_list{m};
    if strcmp(meas,'MD')
      vol = 1000*mean(DTmeas.volE,4);
      vals = vol(ind_brain);
    else
      vals = DTmeas.(['vol' meas])(ind_brain);
    end;
    %% todo: flatten these names; e.g., median_brain_FA, median_brain_MD
    %% todo: remove mean and std?
    brain_DTmeas.(meas) = [];
    brain_DTmeas.(meas).mean = mean(vals);
    brain_DTmeas.(meas).median = median(vals);
    brain_DTmeas.(meas).std = std(vals);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RSIfit = calc_RSIfit(data,parms)
  nframes = size(data.qmat,1);
  if parms.verbose
    fprintf('%s: running RSI fit calculations with %d frames...\n',...
      mfilename,nframes);
    tic
  end
  
  % specify RSI fit parameters
  p = [];
  p.bvals = data.bvals;
  p.M = data.M;
  p.lambda = parms.lambda;
  p.ADC_long = 1e-3; %% todo: allow as input parameter?
  p.ADC_trans_min = 0;
  p.ADC_trans_max = 0.9e-3; %% todo: allow as input parameter?
  p.num_ADC_trans = 2; %% todo: allow as input parameter?
  p.SH_order = 4; %% todo: allow as input parameter?
  p.iso_restricted_flag = 0;
  p.iso_hindered_flag = 0;
  p.iso_free_flag = 1;
  p.ADC_free = 3e-3;
  
  % calculate RSI fit (multiple signal fractions modeled as FODs)
  a = mmil_parms2args(p);
  RSIfit = rsi_fit_MFOD(data.vol,data.qmat,a{:});
  
  if parms.verbose, toc; end
  if isempty(RSIfit)
    fprintf('%s: WARNING: RSI fit calculations failed\n',mfilename);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RSImeas = calc_RSImeas(RSIfit,parms)
  RSImeas = [];
  if parms.verbose
    fprintf('%s: calculating RSI measures...\n',mfilename);
    tic
  end
  % calculate N0, ND, etc.
  RSImeas = rsi_calc_MFmeas(RSIfit,'mask_flag',0);
  if parms.verbose, toc; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RSIerr = calc_RSIerr(RSIfit,data,parms)
  RSIerr = [];
  if parms.verbose
    fprintf('%s: calculating RSI error...\n',mfilename);
    tic
  end

  % calculate synthesized volume
  vol_synth = rsi_synth_vol(RSIfit);

  RSIerr.unique_bvals = data.unique_bvals;
  RSIerr.nbvals = data.nbvals;
  RSIerr.M = data.M;

  % calculate rms residual error overall and for each b-value
  nslices = RSIfit.volsz(3);
  volsz = [RSIfit.volsz(1:3)];
  vol_err = zeros(RSIfit.volsz);
  vol_rms = zeros(volsz);
  vol_rms_err = zeros(volsz);
  vol_rms_err_norm = zeros(volsz);
  vol_rms_err_rel = zeros(volsz);
  volsz_bvals = [RSIfit.volsz(1:3),data.nbvals];
  vol_rms_bvals = zeros(volsz);
  vol_rms_err_bvals = zeros(volsz);
  vol_rms_err_norm_bvals = zeros(volsz);
  vol_rms_err_rel_bvals = zeros(volsz);
  
  % calculate RMS signal and error across frames
  % one slice at a time to save memory
  for z = 1:nslices
    % calculate RMS error across all frames
    slc = data.vol(:,:,z,:);
    slc_synth = vol_synth(:,:,z,:);
    vol_err(:,:,z,:) = slc_synth-slc;
    vol_rms(:,:,z) = sqrt(mean(slc.^2,4));
    vol_rms_err(:,:,z) = sqrt(mean((slc_synth-slc).^2,4));
    % calculate RMS error across frames for each b-value
    for b=1:data.nbvals
      bval = data.unique_bvals(b);
      idx = find(data.round_bvals==bval);
      slc = data.vol(:,:,z,idx);
      slc_synth = vol_synth(:,:,z,idx);
      vol_rms_bvals(:,:,z,b) = sqrt(mean(slc.^2,4));
      vol_rms_err_bvals(:,:,z,b) = sqrt(mean((slc_synth-slc).^2,4));
    end
  end

  % calculate median rms value of b=0 frames: the common normalizer for the
  %   overall error and for every b-value, as in calc_DTerr
  vol_tmp = vol_rms_bvals(:,:,:,1);
  med_rms_b0 = median(mmil_rowvec(vol_tmp(data.vol_mask>0)));
  % normalize err by the b=0 median rms
  vol_rms_err_norm = vol_rms_err / med_rms_b0;
  % normalize err by rms for each voxel
  vol_rms_err_rel = min(1,vol_rms_err ./ (vol_rms + eps));
  % normalize error for each b-value
  for b=1:data.nbvals
    vol_rms_err_norm_bvals(:,:,:,b) = vol_rms_err_bvals(:,:,:,b) / med_rms_b0;
    vol_rms_err_rel_bvals(:,:,:,b) = min(1,vol_rms_err_bvals(:,:,:,b) ./ (vol_rms_bvals(:,:,:,b) + eps));
  end
      
  % save output volumes
  RSIerr.vol = vol_err;
  RSIerr.vol_rms_err_norm = vol_rms_err_norm;
  RSIerr.vol_rms_err_rel = vol_rms_err_rel;
  RSIerr.vol_rms_err_norm_bvals = vol_rms_err_norm_bvals;
  RSIerr.vol_rms_err_rel_bvals = vol_rms_err_rel_bvals;

  if parms.verbose, toc; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_RSIerr = calc_brain_RSIerr(RSIerr,data,parms)
  brain_RSIerr = [];
  if parms.verbose
    fprintf('%s: calculating RSI residual error within brain...\n',mfilename);
  end

  % get brain voxels from brainmask
  ind_brain = find(data.vol_mask>0);

  % calculate median residual error within brain
  brain_RSIerr.median_brain_rms_err_norm = median(RSIerr.vol_rms_err_norm(ind_brain));
  brain_RSIerr.median_brain_rms_err_rel = median(RSIerr.vol_rms_err_rel(ind_brain));
  
  % calculate median residual error within brain for each b-value
  brain_RSIerr.median_brain_rms_err_norm_bvals = zeros(1,RSIerr.nbvals);
  brain_RSIerr.median_brain_rms_err_rel_bvals = zeros(1,RSIerr.nbvals);
  for b=1:RSIerr.nbvals
    tmp_err = RSIerr.vol_rms_err_norm_bvals(:,:,:,b);
    brain_RSIerr.median_brain_rms_err_norm_bvals(b) = median(tmp_err(ind_brain));
    tmp_err = RSIerr.vol_rms_err_rel_bvals(:,:,:,b);
    brain_RSIerr.median_brain_rms_err_rel_bvals(b) = median(tmp_err(ind_brain));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: calculate metrics within wm, gm, and ventricles
function brain_RSImeas = calc_brain_RSImeas(RSImeas,data,parms)
  brain_RSImeas = [];
  if parms.verbose
    fprintf('%s: calculating RSI measures within brain...\n',mfilename);
  end;
  ind_brain = find(data.vol_mask>0);
  for m=1:length(parms.RSI_meas_list)
    meas = parms.RSI_meas_list{m};
    switch meas
      case 'N0'
        vol = RSImeas.volN0(:,:,:,1);
      case 'ND'
        vol = RSImeas.volND(:,:,:,1);
      case 'N0s2'
        vol = RSImeas.volN0(:,:,:,2);
      case 'NDs2'
        vol = RSImeas.volND(:,:,:,2);
      case 'NIf'
        vol = RSImeas.volNIf;
      otherwise
        error('unsupported RSI measure %s',meas);
    end
    vals = vol(ind_brain);
    %% todo: flatten these names? e.g., median_brain_N0, median_brain_ND
    %% todo: remove mean and std?
    brain_RSImeas.(meas) = [];
    brain_RSImeas.(meas).mean = mean(vals);
    brain_RSImeas.(meas).median = median(vals);
    brain_RSImeas.(meas).std = std(vals);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cen_fit,cen_err,cen_data] = set_censor_input(fit,err,data,parms,model)
  % select fit, err, and data for censoring detection
  % use a fit to the raw data when supplied, so detection is slice-aligned
  cen_fit = fit; cen_err = err; cen_data = data;
  if ~isfield(data,'vol_raw') || isempty(data.vol_raw), return; end
  if parms.verbose
    fprintf('%s: detecting %s censoring from raw data...\n',mfilename,model);
  end
  data_raw = set_raw_data(data);
  switch model
    case 'DTI'
      fit_raw = calc_DTfit(data_raw,parms);
    case 'RSI'
      fit_raw = calc_RSIfit(data_raw,parms);
  end
  if isempty(fit_raw)
    fprintf('%s: WARNING: %s fit to raw data failed; detecting censoring from data.vol\n',...
      mfilename,model);
    return;
  end
  fit_raw.volmask = data.vol_mask;
  switch model
    case 'DTI'
      err_raw = calc_DTerr(fit_raw,data_raw,parms);
    case 'RSI'
      err_raw = calc_RSIerr(fit_raw,data_raw,parms);
  end
  cen_fit = fit_raw; cen_err = err_raw; cen_data = data_raw;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_raw = set_raw_data(data)
  % construct data struct for the raw data, for censoring detection
  data_raw = data;
  data_raw.vol = data.vol_raw;
  data_raw.vol_raw = [];
  if isfield(data,'qmat_raw') && ~isempty(data.qmat_raw)
    data_raw.qmat = data.qmat_raw;
  end
  if isfield(data,'bvals_raw') && ~isempty(data.bvals_raw)
    data_raw.bvals = data.bvals_raw;
    [data_raw.round_bvals,data_raw.unique_bvals,data_raw.nbvals] = ...
      get_unique_bvals(data_raw.bvals);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [censor_info,censor_data] = calc_censor_info(fit,err,data,parms)
  censor_info = [];
  args = mmil_parms2args(parms,parms.info_tags);
  [censor_info,censor_data,err_vol_norm] = autoqc_dMRI_censor_info(fit,err,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_err_images(DTerr,data,parms,bids_stem,model)
  if ~exist('bids_stem','var'), bids_stem = ''; end
  if ~exist('model','var') || isempty(model), model = 'DTI'; end
  if parms.bids_flag
    if isempty(bids_stem)
      bids_prefix = sprintf('desc-%s',model);
    else
      bids_prefix = sprintf('%s_desc-%s',bids_stem,model);
    end
    if parms.plot_extra_flag
      vol_tmp = DTerr.vol_rms_err_norm .* data.vol_mask;
      plot_threeviews(vol_tmp,DTerr.M,sprintf('%s_normRMSE',bids_prefix),parms,parms.normmax);
    end
    vol_tmp = DTerr.vol_rms_err_rel .* data.vol_mask;
    plot_threeviews(vol_tmp,DTerr.M,sprintf('%s_relRMSE',bids_prefix),parms,parms.relmax);
    for b=1:DTerr.nbvals
      bval = DTerr.unique_bvals(b);
      if parms.plot_extra_flag
        vol_tmp = DTerr.vol_rms_err_norm_bvals(:,:,:,b) .* data.vol_mask;
        plot_threeviews(vol_tmp,DTerr.M,sprintf('%sb%d_normRMSE',bids_prefix,bval),parms,parms.normmax);
      end
      vol_tmp = DTerr.vol_rms_err_rel_bvals(:,:,:,b) .* data.vol_mask;
      plot_threeviews(vol_tmp,DTerr.M,sprintf('%sb%d_relRMSE',bids_prefix,bval),parms,parms.relmax);
    end
  else
    prefix = bids_stem;
    if isempty(prefix), prefix = upper(model); end
    if parms.plot_extra_flag
      vol_tmp = DTerr.vol_rms_err_norm .* data.vol_mask;
      plot_threeviews(vol_tmp,DTerr.M,sprintf('%s_rms_err_norm',prefix),parms,parms.normmax);
    end
    vol_tmp = DTerr.vol_rms_err_rel .* data.vol_mask;
    plot_threeviews(vol_tmp,DTerr.M,sprintf('%s_rms_err_rel',prefix),parms,parms.relmax);
    for b=1:DTerr.nbvals
      bval = DTerr.unique_bvals(b);
      if parms.plot_extra_flag
        vol_tmp = DTerr.vol_rms_err_norm_bvals(:,:,:,b) .* data.vol_mask;
        plot_threeviews(vol_tmp,DTerr.M,sprintf('%s_rms_err_norm_b%d',prefix,bval),parms,parms.normmax);
      end
      vol_tmp = DTerr.vol_rms_err_rel_bvals(:,:,:,b) .* data.vol_mask;
      plot_threeviews(vol_tmp,DTerr.M,sprintf('%s_rms_err_rel_b%d',prefix,bval),parms,parms.relmax);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_censor_images(censor_info,censor_data,parms,fname_out)
  if ~exist('fname_out','var') || isempty(fname_out)
    fname_out = sprintf('%s/dtiCensor.jpg',parms.outdir);
  end
  if ~exist(fname_out,'file') || parms.forceflag
    figure; clf;
    colormap hot
    subplot(1,2,1);
    imagesc(censor_data.censor_err,[1,censor_info.censor_thresh+2]);
    axis xy;
    title('normalized error','FontSize',parms.fontsize);
    xlabel('frame','FontSize',parms.fontsize);
    ylabel('slice','FontSize',parms.fontsize);
    set(gca,'FontSize',parms.fontsize);
    subplot(1,2,2);
    imagesc(censor_data.censor_mat);
    axis xy;
    title('censored frame-slices','FontSize',parms.fontsize);
    xlabel('frame','FontSize',parms.fontsize);
    ylabel('slice','FontSize',parms.fontsize);
    set(gca,'FontSize',parms.fontsize);
    mmil_save_fig(gcf,fname_out,parms.censor_fig_size,parms.tif_dpi);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function plot_threeviews(vol,M,outstem,parms,valrange)
  if ~exist('valrange','var'), valrange = []; end
  tparms = [];
  tparms.outdir = parms.outdir;
  if parms.bids_flag
    tparms.outstem = outstem;
    fname_out = sprintf('%s/%s.jpg',parms.outdir,outstem);
  else
    tparms.outstem = sprintf('%s_tv',outstem);
    fname_out = sprintf('%s/%s.jpg',parms.outdir,tparms.outstem);
  end
  if ~exist(fname_out,'file') || parms.forceflag
    fname = sprintf('%s/%s.mgh',parms.tmpdir,outstem);
    if ~exist(fname,'file') || parms.forceflag
      fs_save_mgh(vol,fname,M);
    end
    if isempty(valrange)
      tparms.plim = parms.threeviews_plim;
    else
      tparms.plim_flag = 0;
      tparms.clim = [0,valrange];
    end
    tparms.fig_size = parms.threeviews_fig_size;
    tparms.tif_dpi = parms.tif_dpi;
    tparms.planestrings = parms.planestrings;
    tparms.slice_fract = parms.threeviews_slice_fract;
    tparms.colorbar_flag = parms.colorbar_flag;
    tparms.forceflag = parms.forceflag;
    args = mmil_parms2args(tparms);
    mmil_subplots(fname,args{:});
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if exist(parms.tmpdir,'dir') && ~strcmp(parms.tmpdir,parms.indir) &&...
                                  ~strcmp(parms.tmpdir,parms.outdir)
    fprintf('%s: deleting %s...\n',mfilename,parms.tmpdir);
    rmdir(parms.tmpdir,'s');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
