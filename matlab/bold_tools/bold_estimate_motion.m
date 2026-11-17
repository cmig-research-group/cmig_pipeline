function [M_motion,motion_tseries,start_info] = bold_estimate_motion(fname,varargin)
%function [M_motion,motion_tseries,start_info] = bold_estimate_motion(fname,[options])
%
% Purpose: Estimate head motion for BOLD scans using
%            rigid-body volume registration
%
% Between-scan registration (fname_ref differs from fname):
%   the first registration runs from two starts, the header and the
%   brain-mask centroids aligned (cmig_create_brainmask of the reference
%   and of the input mean); if the two results are within 3 mm (mean
%   displacement over the reference brain mask), the header result is kept;
%   otherwise the result with the lower final cost, with a WARNING
%
% Usage:
%  [M_motion,motion_tseries,start_info] = bold_estimate_motion(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Input:
%  'intra_flag': estimate within-scan motion
%     {default = 1}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    If empty, will use fname as reference
%    {default = []}
%  'ref_frame': reference frame number (1-based)
%    if ref_frame = 0, will use mean across frames, excluding first skipTRs
%    {default = 1}
%  'skipTRs': number of TRs at beginning of scan to be ignored in calculating mean
%    only used if ref_frame = 0
%    {default = 0}
%  'fname_mask': name of mask volume
%     If empty, will generate one from fname_ref
%     {default = []}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%  'verbose': [0|1] display status messages
%    {default = 1}
%
% Output:
%   M_motion: cell array containing 4x4 rigid-body transformation matrices
%             for each frame
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%         also called: dL, dP, dS, pitch, yaw, roll
%                      to match output of mmil_load_motion_1D
%     with translations in mm and rotations in degrees clockwise
%   start_info: struct with fields
%     start: start of the first registration kept ('header' or 'centroid')
%     distance: mean displacement (mm) over the reference brain mask between
%       the header and centroid results; NaN if only the header start ran
%       (within-scan, or an empty brain mask)
%
% Created:  08/27/12 by Don Hagler
% Prev Mod: 02/02/21 by Don Hagler
% Prev Mod: 05/26/21 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Last Mod: 11/12/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_motion = []; motion_tseries = [];
start_info = struct('start','header','distance',NaN);
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters, check files exist
parms = check_input(fname,varargin);

% load input volume, reference, and mask
parms = load_data(parms);

% estimate head motion
[M_motion,start_info] = estimate_motion(parms);

% derive motion parameter time series from registration matrices
motion_tseries = set_motion_tseries(parms,M_motion);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'intra_flag',true,[false true],...
    'fname_ref',fname,[],...
    'ref_frame',1,[0 Inf],...
    'skipTRs',0,[0,Inf],...
    'fname_mask',[],[],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'verbose',true,[false true],...
  ...
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
    'start_maxdiff',3,[0,Inf],... % mm
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end;
  if ~isempty(parms.fname_mask) && ~exist(parms.fname_mask,'file'), 
    error('file %s not found',parms.fname_mask);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = load_data(parms)
  % load input file
  [parms.vol,parms.M,mrp,parms.volsz] = fs_load_mgh(parms.fname_in);
  parms.nframes = size(parms.vol,4);

  % check ref_frame
  [parms.M_ref,parms.volsz_ref] = fs_read_header(parms.fname_ref);
  if parms.ref_frame~=1
    if parms.ref_frame==0 % indicates mean as reference
      % ensure skipTRs < num frames - 1
      parms.skipTRs = min(parms.skipTRs,parms.volsz_ref(4)-1);
    else % single frame reference, other than first
      % ensure ref_frame is less than num frames
      parms.ref_frame = min(parms.ref_frame,parms.volsz_ref(4));
    end;
  end;
  
  % load reference volume
  if parms.ref_frame>0
    if ~strcmp(parms.fname_in,parms.fname_ref)
      [vol_ref,M_ref] =...
        fs_load_mgh(parms.fname_ref,[],parms.ref_frame);
    else
      vol_ref = parms.vol(:,:,:,parms.ref_frame);
      M_ref = parms.M;
    end;
  else % calculate mean if ref_frame = 0
    if parms.verbose
      fprintf('%s: calculating mean of reference volume %s...\n',...
        mfilename,parms.fname_ref);
      tic;
    end;
    if ~strcmp(parms.fname_in,parms.fname_ref)
      [vol_ref,M_ref] =...
        fs_load_mgh(parms.fname_ref);
    else
      vol_ref = parms.vol;
      M_ref = parms.M;
    end;
    if parms.skipTRs>0
      vol_ref = vol_ref(:,:,:,parms.skipTRs+1:end);
    end;
    vol_ref = mean(vol_ref,4);
    if parms.verbose, toc; end;
  end;
  parms.vol_ref_ctx = ctx_mgh2ctx(vol_ref,M_ref);
    
  % load mask volume or create from reference volume
  if ~isempty(parms.fname_mask)
    parms.vol_mask_ctx = ctx_load_mgh(parms.fname_mask);
  else
    fname_tmp_vol = [mmil_tempfname('bold_estimate_motion') '.mgz'];
    fname_tmp_mask = [mmil_tempfname('bold_estimate_motion') '.mgz'];
    ctx_save_mgh(parms.vol_ref_ctx,fname_tmp_vol);
    cmig_create_brainmask(fname_tmp_vol,'fname_mask',fname_tmp_mask,...
      'binarize_flag',1);
    parms.vol_mask_ctx = ctx_load_mgh(fname_tmp_mask);
    delete(fname_tmp_vol);
    delete(fname_tmp_mask);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_motion,start_info] = estimate_motion(parms)
  M_motion = cell(parms.nframes,1);
  start_info = struct('start','header','distance',NaN);

  % between-scan registration (reference is another scan)
  between_flag = ~strcmp(parms.fname_in,parms.fname_ref);

  % start from the header
  M_adj = parms.M;
  M_init = eye(4);

  % estimate motion for each frame
  for f=1:parms.nframes
    %% todo: only if intra_flag=true
    if f==1 || parms.intra_flag
      %% todo: reduce number of scales for faster reg on each frame (if intra_flag=1)
      if parms.ref_frame==0 && ~parms.intra_flag
        % register mean image to ref
        vol = mean(parms.vol(:,:,:,parms.skipTRs+1:end),4);
      else
        % register frame to ref
        vol = squeeze(parms.vol(:,:,:,f));
      end;
      vol_ctx = ctx_mgh2ctx(vol,M_adj);
      [M_reg,min_cost] = rbreg_vol2vol_icorr_mask(parms.vol_ref_ctx,vol_ctx,...
        parms.vol_mask_ctx,0,parms.mstep,0,parms.scales,...
        parms.min_trans,parms.min_rot*pi/180);
      if f==1 && between_flag
        % also start from the brain-mask centroids, keep the lower cost
        [M_adj,M_init,M_reg,start_info] = check_centroid_start(parms,vol,M_reg,min_cost);
      end;
    end;
    M_motion{f} = M_init*M_reg;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_adj,M_init,M_reg,start_info] = check_centroid_start(parms,vol,M_reg_hdr,cost_hdr)
  % header start (default)
  M_adj = parms.M;
  M_init = eye(4);
  M_reg = M_reg_hdr;
  start_info = struct('start','header','distance',NaN);

  % brain-mask centroid of the reference
  [vol_mask_ref,M_mask_ref] = ctx_ctx2mgh(parms.vol_mask_ctx);
  vol_mask_ref = vol_mask_ref>0;
  ref_cen = calc_centroid(vol_mask_ref,M_mask_ref);

  % brain-mask centroid of the input mean
  vol_mean = mean(parms.vol(:,:,:,min(parms.skipTRs,parms.nframes-1)+1:end),4);
  fname_tmp_vol = [mmil_tempfname('bold_estimate_motion') '.mgz'];
  fname_tmp_mask = [mmil_tempfname('bold_estimate_motion') '.mgz'];
  fs_save_mgh(vol_mean,fname_tmp_vol,parms.M);
  cmig_create_brainmask(fname_tmp_vol,'fname_mask',fname_tmp_mask,...
    'binarize_flag',1);
  vol_mask_in = fs_load_mgh(fname_tmp_mask)>0;
  delete(fname_tmp_vol);
  delete(fname_tmp_mask);
  in_cen = calc_centroid(vol_mask_in,parms.M);
  if any(isnan(in_cen)) || any(isnan(ref_cen))
    fprintf('%s: WARNING: empty brain mask, using header start only for %s\n',...
      mfilename,parms.fname_in);
    return;
  end;

  % shift input header by centroid offset and register
  offset = ref_cen - in_cen;
  M_adj_cen = parms.M;
  M_adj_cen(1:3,4) = M_adj_cen(1:3,4) + offset;
  M_init_cen = M_RAS_TO_LPH*parms.M*inv(M_RAS_TO_LPH*M_adj_cen);
  vol_ctx = ctx_mgh2ctx(vol,M_adj_cen);
  [M_reg_cen,cost_cen] = rbreg_vol2vol_icorr_mask(parms.vol_ref_ctx,vol_ctx,...
    parms.vol_mask_ctx,0,parms.mstep,0,parms.scales,...
    parms.min_trans,parms.min_rot*pi/180);

  % mean displacement of reference mask voxels between the two results
  [i,j,k] = ind2sub(size(vol_mask_ref),find(vol_mask_ref));
  lph = M_RAS_TO_LPH*M_mask_ref*[i';j';k';ones(1,length(i))];
  D = (M_init_cen*M_reg_cen)*inv(M_reg_hdr) - eye(4);
  dmean = mean(sqrt(sum((D(1:3,:)*lph).^2,1)));

  start_info.distance = dmean;

  % header result if the two agree, else the lower final cost
  if dmean <= parms.start_maxdiff
    if parms.verbose
      fprintf(['%s: kept start: header (centroid agrees) (cost header %0.4f, ',...
        'centroid %0.4f; %0.1f mm apart; centroid offset %0.1f %0.1f %0.1f mm)\n'],...
        mfilename,cost_hdr,cost_cen,dmean,offset);
    end;
    return;
  end;
  if cost_cen < cost_hdr
    start_info.start = 'centroid';
    M_adj = M_adj_cen;
    M_init = M_init_cen;
    M_reg = M_reg_cen;
  end;
  fprintf(['%s: WARNING: header and brain-mask centroid starts end %0.1f mm apart ',...
    '(mean over reference brain mask) for %s with %s as reference; ',...
    'kept %s start (cost header %0.4f, centroid %0.4f)\n'],...
    mfilename,dmean,parms.fname_in,parms.fname_ref,start_info.start,cost_hdr,cost_cen);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = calc_centroid(vol_mask,M)
  % mean voxel index mapped through vox2ras (RAS)
  [i,j,k] = ind2sub(size(vol_mask),find(vol_mask));
  c = M(1:3,:)*[mean(i);mean(j);mean(k);1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_tseries = set_motion_tseries(parms,M_motion)
  motion_tseries = zeros(parms.nframes,6);
  for f=1:parms.nframes
    tmp_motion = mmil_M_mat2vec(inv(M_motion{f}));
    motion_tseries(f,:) = tmp_motion;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

