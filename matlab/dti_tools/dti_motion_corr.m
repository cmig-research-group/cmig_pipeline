function [fname_out,fname_qmat,censor_info] = dti_correct_motion(fname,varargin)
%function [fname_out,fname_qmat] = dti_correct_motion(fname,[options])
%
% Purpose: Correct head motion for DTI scans using
%            rigid-body volume registration between data and tensor fit
%
% Usage:
%  [fname_out,fname_qmat] = dti_correct_motion(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Parameters:
%  'fname_out': output file name
%    if empty, will add '_mc' to file stem of input file name
%    {default = []}
%  'fname_qmat': output file name containing adjusted qmat
%     and motion estimates
%    if empty, will add '_mc_qmat' to file stem of input file name
%    {default = []}
%  'interpm' - interpolation method
%     0 = nearest neighbor, 1 = linear, 2 = cubic
%     3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%    { default: 2 }
%  'maskoutput': apply brain mask to set background to zero
%    {default = 0}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': overwrite existing output files
%     {default: 0}
%
% Optional Parameters for motion estimation:
%  'intra_flag': perform within-scan motion correction
%     {default = 1}
%  'qmat': matrix of diffusion direction vectors
%    if empty, will do between-scan estimation only
%    {default = []}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    if empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%    if empty, will generate one from fname_ref using cmig_create_brainmask
%    {default = []}
%  'fname_reg': output file name for between-scan registration matrix
%    {default = []}
%  'censor_niter': number of iterations of censoring
%    {default = 1}
%  'censor_min_ndirs': minimum number of diffusion directions (not including
%    b=0 images) required for tensor fit after censoring
%    will not do censoring if it means reducing number of directions below min
%    {default = 6}
%  'censor_thresh': error threshold for censoring bad frames
%    normalized to median error for each slice
%    higher values mean less censoring
%    {default = 3.2}
%  'censor_mat': matrix of slices by frame to exclude from fit
%    ignored if censor_niter = 0
%    {default = []}
%  'nonlin_flag': [0|1] use nonlinear optimization
%     with initial parameters from linear fit
%     {default = 0}
%  'b0_thresh': threshold used for considering a b-value to be 0
%    {default = 10}
%  'bvals': vector of b values (one for each diffusion direction)
%    If single value supplied, will use same for all
%    {default = [1000]}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%  'afni_flag': [0|1] use AFNI's 3dvolreg for within-scan registration
%     otherwise use rbreg_vol2vol_icorr_mask
%     {default = 0}
%  'censor_reg_flag': [0|1] exclude censored slices from the per-frame
%     registration instead of replacing them with synthesized values
%     before registration; replacement then happens after alignment:
%     the tensor is refit on the aligned data, each frame's synthesized
%     image is resampled back to its acquisition space with the inverse
%     motion estimate, outlier slices are detected there (slice-wise
%     residual vs the frame's actual data, normalized by median slice
%     error, same censor_thresh), and those slices are replaced before
%     the frame is realigned
%     ignored if censor_niter = 0
%     {default = 0}
%  'tmpdir': directory for temporary files used by AFNI registration
%     if empty, will use directory containing fname
%     {default = []}
%
% Output:
%   fname_out: motion correct mgh/mgz file containing 4D diffusion volume
%   fname_qmat: mat file containing:
%    'qmat': diffusion direction matrix (ndirs x 3) adjusted for rotation
%      from motion correction
%    'M_motion': cell array containing registration matrices M_ref_to_orig
%      for each frame
%   censor_info: struct with censor_mat and censor_err from the
%     post-alignment acquisition-space detection
%     empty unless censor_reg_flag = 1 and motion correction ran
%
% Created:  02/11/13 by Don Hagler
% Prev Mod: 03/23/16 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Prev Mod: 08/12/2026 by Don Hagler
% Last Mod: 08/26/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_qmat = []; censor_info = [];
if ~mmil_check_nargs(nargin,2), return; end;

parms = check_input(fname,varargin);

if output_exists(parms)
  fname_out = parms.fname_out;
  fname_qmat = parms.fname_qmat;
  % recover the post-alignment censoring from its intermediate file
  if parms.censor_reg_flag
    fname_censor = censor_fname(parms.fname_qmat);
    if exist(fname_censor,'file')
      censor_info = load(fname_censor);
    end;
  end;
  return;
end;

[parms,fname_qmat] = estimate_motion(parms);

[fname_out,censor_info] = correct_motion(parms);

% save the post-alignment censoring as an intermediate next to the
%   motion file, so skipped reruns can recover it without recomputing
if ~isempty(censor_info)
  fname_censor = censor_fname(parms.fname_qmat);
  censor_mat = censor_info.censor_mat;
  censor_err = censor_info.censor_err;
  save(fname_censor,'censor_mat','censor_err');
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'fname_out',[],[],...
    'fname_qmat',[],[],...
    'maskoutput',false,[false true],...
    'interpm',2,[0:5],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'intra_flag',true,[false true],...
    'qmat',[],[],...
    'fname_ref',fname,[],...
    'fname_mask',[],[],...
    'fname_reg',[],[],...
    'censor_niter',1,[0,100],...
    'censor_min_ndirs',6,[],...
    'censor_thresh',3.2,[],...
    'censor_mat',[],[],...
    'nonlin_flag',false,[false true],...
    'b0_thresh',10,[0,500],...
    'bvals',1000,[],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'afni_flag',false,[false true],...
    'censor_reg_flag',false,[false true],...
    'tmpdir',[],[],...
  ... % hidden
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
    'smf',1e-5,[1e-100,1e-1],...
    'motion_radius',50,[],... % for calculating distance from angle
  ...
    'estimate_tags',{'intra_flag','qmat',...
                     'fname_ref','fname_mask','fname_reg','censor_niter',...
                     'censor_min_ndirs','censor_thresh','censor_mat',...
                     'nonlin_flag','b0_thresh',...
                     'bvals','min_trans','min_rot','afni_flag',...
                     'censor_reg_flag','tmpdir',...
                     'verbose','forceflag','mstep','scales','smf'},[],...
    'censor_tags',{'volmask','censor_min_ndirs','nb0'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(fname,'file')
    error('file %s not found',fname);
  end;
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end;

  [fpath,fstem,fext] = fileparts(parms.fname_in);
  if ~isempty(parms.fname_out)
    parms.outdir = fileparts(parms.fname_out);
  else
    parms.outdir = fpath;
  end;
  if isempty(parms.fname_out)
    parms.fname_out = [parms.outdir '/' fstem '_mc' fext];
  end;
  if isempty(parms.fname_qmat)
    parms.fname_qmat = [parms.outdir '/' fstem '_mc_qmat.mat'];
  end;

  % get info about input volume
  [parms.M,parms.volsz] = ...
    mmil_load_mgh_info(parms.fname_in,parms.forceflag);
  parms.nframes = parms.volsz(4);

  % get info about input volume
  [parms.M_ref,parms.volsz_ref] = ...
    mmil_load_mgh_info(parms.fname_ref,parms.forceflag);

  % create brainmask if needed
  if parms.maskoutput
    if ~isempty(parms.fname_mask)
      [fpath,fstem_ref,fext] = fileparts(parms.fname_ref);
      parms.fname_mask = [parms.outdir '/' fstem_ref '_mask' fext];
    end;
    if ~exist(parms.fname_mask,'file') || parms.forceflag
      cmig_create_brainmask(parms.fname_ref,'fname_mask',parms.fname_mask,...
        'binarize_flag',1);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end;
  if ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if ~exist(parms.fname_qmat,'file')
    exist_flag = 0;
  end;
  if parms.censor_reg_flag && parms.censor_niter>0 && ...
     parms.intra_flag && ~isempty(parms.qmat) && ...
     ~exist(censor_fname(parms.fname_qmat),'file')
    exist_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = censor_fname(fname_qmat)
  % name of the censoring intermediate paired with a motion file
  fname = regexprep(fname_qmat,'_qmat\.mat$','_censor.mat');
  if strcmp(fname,fname_qmat)
    fname = regexprep(fname_qmat,'\.mat$','_censor.mat');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname_qmat] = estimate_motion(parms)
  if ~exist(parms.fname_qmat,'file') || parms.forceflag
    args = mmil_parms2args(parms,parms.estimate_tags);
    [M_motion,motion_tseries,qmat,DTfit] = ...
      dti_estimate_motion(parms.fname_in,args{:});
    if ~isempty(parms.qmat)
      bvals = parms.bvals;
      orig_qmat = parms.qmat;
      % calculate absolute difference in head position or rotation
      [mean_motion,mean_trans,mean_rot] = ...
        mmil_mean_motion(motion_tseries,parms.motion_radius);
      save(parms.fname_qmat,'M_motion','motion_tseries',...
        'mean_motion','mean_trans','mean_rot',...
        'qmat','orig_qmat','bvals','DTfit');
    end;
  else
    load(parms.fname_qmat);
  end;
  parms.M_motion = M_motion;
  parms.DTfit = DTfit;
  parms.qmat_adj = qmat;
  fname_qmat = parms.fname_qmat;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,censor_info] = correct_motion(parms)
  censor_info = [];
  if ~exist(parms.fname_out,'file') || parms.forceflag
    [vol_orig,M] = fs_load_mgh(parms.fname_in);
    % censor_reg_flag: censored slices were excluded from registration,
    %   so replace them AFTER alignment (see replace_censored_aligned);
    %   otherwise replace before the resample (legacy behavior)
    censor_post = parms.censor_reg_flag && parms.censor_niter>0 && ...
                  ~isempty(parms.DTfit) && ~isempty(parms.qmat);
    % replace censored slices in vol with those synthesized from tensor
    if parms.censor_niter>0 && ~isempty(parms.DTfit) && ~censor_post
      vol_synth = dti_synth_vol(parms.DTfit);
      parms.volmask = parms.DTfit.volmask;
      parms.nb0 = length(parms.DTfit.i_b0);
      args = mmil_parms2args(parms,parms.censor_tags);
      vol_orig = dti_censor_vol(vol_orig,vol_synth,...
                                parms.DTfit.censor_mat,args{:});
    end;
    % initialize output
    vol = zeros([parms.volsz_ref(1:3),parms.nframes],'single');
    if parms.maskoutput
      volmask = fs_load_mgh(parms.fname_mask);
    end;
    if parms.verbose
      if parms.nframes>1
        plurstr = 's';
      else
        plurstr = [];
      end;      
      fprintf('%s: resampling %d frame%s...\n',...
        mfilename,parms.nframes,plurstr);
    end;
    for f = 1:parms.nframes
      vol_tmp = squeeze(vol_orig(:,:,:,f));
      vol_tmp = mmil_resample_vol(vol_tmp,M,...
        'M_ref',parms.M_ref,'nvox_ref',parms.volsz_ref,...
        'M_reg',parms.M_motion{f},...
        'interpm',parms.interpm);
      if parms.maskoutput
        vol_tmp(~volmask) = 0;
      end;
      vol(:,:,:,f) = vol_tmp;
    end;
    % replace censored slices after alignment
    if censor_post
      [vol,censor_info] = replace_censored_aligned(parms,vol,vol_orig,M);
    end;
    fs_save_mgh(vol,parms.fname_out,parms.M_ref);
  end;
  fname_out = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,censor_info] = replace_censored_aligned(parms,vol,vol_orig,M)
  % refit tensor on the aligned data with the motion-adjusted qmat,
  % then detect outlier slices IN EACH FRAME'S ACQUISITION SPACE:
  % resample the synthesized frame back with the inverse motion estimate
  % and compare against the frame's actual (unaligned) data, slice by
  % slice, so detection and replacement index the same physical slices
  if parms.verbose
    fprintf('%s: replacing censored slices after alignment...\n',mfilename);
  end;
  DTrefit = dti_fit_tensor(double(vol),parms.qmat_adj,...
    'bvals',parms.bvals,'M',parms.M_ref,...
    'censor_niter',parms.censor_niter,...
    'censor_mat',parms.DTfit.censor_mat,...
    'censor_min_ndirs',parms.censor_min_ndirs,...
    'censor_thresh',parms.censor_thresh,...
    'nonlin_flag',parms.nonlin_flag,...
    'b0_thresh',parms.b0_thresh);
  i_b0 = DTrefit.i_b0;
  i_bx = setdiff(1:parms.nframes,i_b0);
  % synthesize volume using tensor fit to aligned data
  vol_synth = dti_synth_vol(DTrefit);
  volmask = 1.0*(DTrefit.volmask>0);
  % pass 1: slice-wise residuals in acquisition space
  %   (normalization mirrors dti_fit_tensor fit_tensor_censor)
  nz = parms.volsz(3);
  mean_err_mat = zeros(nz,parms.nframes);
  for f=1:parms.nframes
    M_inv = inv(parms.M_motion{f});
    synth_f = mmil_resample_vol(squeeze(vol_synth(:,:,:,f)),parms.M_ref,...
      'M_ref',M,'nvox_ref',parms.volsz,...
      'M_reg',M_inv,'interpm',parms.interpm);
    mask_f = mmil_resample_vol(volmask,parms.M_ref,...
      'M_ref',M,'nvox_ref',parms.volsz,...
      'M_reg',M_inv,'interpm',0);
    vol_f = double(squeeze(vol_orig(:,:,:,f)));
    for z=1:nz
      im_mask = mask_f(:,:,z)>0;
      if nnz(im_mask)>100
        im_err = vol_f(:,:,z) - synth_f(:,:,z);
        mean_err_mat(z,f) = sqrt(mean(im_err(im_mask).^2));
      end;
    end;
  end;
  % normalize by median error within slice across diffusion-weighted frames
  med_err = median(mean_err_mat(:,i_bx),2);
  censor_err = mean_err_mat ./ max(eps,med_err*ones(1,parms.nframes));
  censor_err(:,i_b0) = 0;
  censor_mat = 1.0*(censor_err > parms.censor_thresh);
  % first b=0 and b>0 frames anchor the registration: never censored
  censor_mat(:,i_b0(1)) = 0;
  censor_mat(:,i_bx(1)) = 0;
  % skip censoring for slices with too few remaining directions
  %   (as in dti_censor_vol)
  nleft = parms.nframes - sum(censor_mat,2);
  censor_mat(nleft < parms.censor_min_ndirs + length(i_b0),:) = 0;
  censor_info = [];
  censor_info.censor_mat = censor_mat;
  censor_info.censor_err = censor_err;
  if parms.maskoutput
    volmask_out = fs_load_mgh(parms.fname_mask);
  end;
  frames = find(any(censor_mat,1));
  if parms.verbose
    fprintf('%s: patching %d frames with censored slices...\n',...
      mfilename,length(frames));
  end;
  % pass 2: patch each censored frame in its acquisition space, realign
  for f=frames
    M_inv = inv(parms.M_motion{f});
    synth_f = mmil_resample_vol(squeeze(vol_synth(:,:,:,f)),parms.M_ref,...
      'M_ref',M,'nvox_ref',parms.volsz,...
      'M_reg',M_inv,'interpm',parms.interpm);
    mask_f = mmil_resample_vol(volmask,parms.M_ref,...
      'M_ref',M,'nvox_ref',parms.volsz,...
      'M_reg',M_inv,'interpm',0);
    % replace censored slices within brain mask
    vol_f = double(squeeze(vol_orig(:,:,:,f)));
    for z=mmil_rowvec(find(censor_mat(:,f)))
      im_orig = vol_f(:,:,z);
      im_synth = synth_f(:,:,z);
      im_mask = mask_f(:,:,z);
      im_orig(im_mask>0) = im_synth(im_mask>0);
      vol_f(:,:,z) = im_orig;
    end;
    % apply the motion estimate to the patched frame
    vol_tmp = mmil_resample_vol(vol_f,M,...
      'M_ref',parms.M_ref,'nvox_ref',parms.volsz_ref,...
      'M_reg',parms.M_motion{f},...
      'interpm',parms.interpm);
    if parms.maskoutput
      vol_tmp(~volmask_out) = 0;
    end;
    vol(:,:,:,f) = vol_tmp;
  end;
return;
