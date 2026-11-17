function [fname_out,fname_mask,M_atl_to_subj,reg_failed,corr_atl,fvox_miss,corr_atl_fail,fvox_miss_fail] = cmig_reg_T1w_to_atlas(fname_in,varargin)
%function [fname_out,fname_mask,M_atl_to_subj,reg_failed,corr_atl,fvox_miss,corr_atl_fail,fvox_miss_fail] = cmig_reg_T1w_to_atlas(fname_in,[options])
%
% Required Input:
%   fname_in: full path file name of T1-weighted image volume (mgz format)
%
% Optional Parameters (key,'value'):
%  'outdir': output directory
%    if not specified, will use path of fname_in
%    {default = []}
%  'atlas_flag': [0|1] whether to resample output in atlas space
%              (rigid-body registration only)
%    {default = 1}
%  'native_flag': [0|1] whether to keep native resolution
%     otherwise, resample to 1mm, 256^3, LIA
%    {default = 0}
%  'atlas_mask_flag': [0|1] whether to derive brain mask from atlas registration
%    otherwise use mask derived in native space
%    if atlas_mask_flag=1 and atlas_flag=0
%      will perform registration to atlas to derive mask but not resample to atlas
%    {default = 1}
%
% Output:
%   fname_out: full path file name of output volume
%   fname_mask: full path file name of output mask volume
%   M_atl_to_subj: 4x4 registration matrix for atlas to subject space, as
%     applied (for the fallback, the offset that centers the volume)
%   reg_failed: [0|1] whether registration to atlas failed
%   corr_atl: correlation between atlas and subject brain, for the
%     registration applied (for the fallback, of the centered volume)
%   fvox_miss: fraction of brain voxels missing, for the registration applied
%   corr_atl_fail: corr_atl of the failed registration (NaN if not failed)
%   fvox_miss_fail: fvox_miss of the failed registration (NaN if not failed)
%
% Created:  08/07/2025 by Don Hagler
% Prev Mod: 08/13/2025 by Don Hagler
% Prev Mod: 08/14/2025 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Last Mod: 10/30/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Based on MMIL_Register_and_Resample_MRI_Volumes
% Created:  03/18/08 by Don Hagler
% Last Mod: 07/04/25 by Don Hagler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_mask = []; M_atl_to_subj = []; reg_failed = []; corr_atl = []; fvox_miss = [];
corr_atl_fail = NaN; fvox_miss_fail = NaN;
if ~mmil_check_nargs(nargin,1), return; end

% check input parameters
parms = check_input(fname_in,varargin);

% create initial brain mask
fname_mask_init = init_brain_mask(fname_in,parms);

% register volumes
[M_atl_to_subj,fname_atlas_res,fname_atlas_mask_res,fname_mask,reg_failed,corr_atl,fvox_miss,parms] = register_volumes(fname_in,fname_mask_init,parms);

% resample volumes
[fname_out,fname_mask,M_atl_to_subj] = resample_volumes(fname_in,fname_mask,fname_atlas_res,fname_atlas_mask_res,M_atl_to_subj,parms);

% registration attempted but not applied: check the volume as resampled
if ~parms.atlas_flag && ~isnan(reg_failed)
  if reg_failed
    corr_atl_fail = corr_atl; fvox_miss_fail = fvox_miss;
  end
  [~,corr_atl,fvox_miss] = check_atlas_reg(M_atl_to_subj,fname_in,fname_mask_init,fname_out,parms,'_applied');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  parms = mmil_args2parms(options, { ...
    'outdir',[],[],...
    'atlas_flag',true,[false true],...
    'native_flag',false,[false true],...
    'atlas_mask_flag',true,[false true],...
  ...
    'image_range',[0,255],[],...
    'forceflag',false,[false true],...
    'cleanupflag',true,[false true],...
    'hatl_offset_RAS',[0 0 0],[],...
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas_ABCD',[],... % ABCD atlas
  ...%    'atlasname','T1_Atlas/T1_atlas',[],...
  ...
    'reg_atl_niters',2,[1,10],...
    'min_corr_atl',0.3,[0.01,0.99],...
    'max_fvox_miss',0.001,[0.0001,0.99],...
  });
  
  if ~exist(fname_in,'file'), error('file %s not found',fname_in); end
  [parms.indir,parms.fstem,parms.ext] = fileparts(fname_in);
  if isempty(parms.indir), parms.indir = pwd; end
  if isempty(parms.outdir), parms.outdir = parms.indir; end
  if isempty(parms.atlasdir), parms.atlasdir = cmig_getenv('CMIG_PIPELINE_ATLASDIR'); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_mask_init = init_brain_mask(fname_in,parms)
  % create brain mask for reference volume in native space
  fname_mask_init = sprintf('%s/%s_brainmask%s',parms.outdir,parms.fstem,parms.ext);
  if ~exist(fname_mask_init,'file') || parms.forceflag
    fprintf('%s: creating initial brainmask...\n',mfilename);
    tic;
    mmil_brainmask(fname_in,'fname_mask',fname_mask_init,'binarize_flag',1);
    toc;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_atl_to_subj,fname_atlas_res,fname_atlas_mask_res,fname_mask,reg_failed,corr_atl,fvox_miss,parms] = register_volumes(fname_in,fname_mask_init,parms);
  M_atl_to_subj = []; fname_mask = []; reg_failed = NaN; fname_atlas_res = []; fname_atlas_mask_res = [];
  if parms.atlas_flag || parms.atlas_mask_flag
    % registration to atlas
    [M_atl_to_subj,fname_atlas_mask,fname_atlas_res,fname_atlas_mask_res] = reg_to_atlas(fname_in,parms);
    % check registration to atlas
    [reg_failed,corr_atl,fvox_miss] = check_atlas_reg(M_atl_to_subj,fname_in,fname_mask_init,fname_atlas_res,parms);
    if reg_failed
      fprintf('%s: WARNING: registration to atlas failed... setting atlas_flag = 0\n',mfilename);
      parms.atlas_flag = 0;
      % use initial mask instead of atlas reg derived mask
      fname_mask = fname_mask_init;
      M_atl_to_subj = eye(4);
    else
      fname_mask = fname_atlas_mask;
    end
  else
    fname_mask = fname_mask_init;
    M_atl_to_subj = eye(4);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M_atl_to_subj,fname_atlas_mask,fname_atlas_res,fname_atlas_mask_res] = reg_to_atlas(fname_in,parms)
  M_atl_to_subj = []; fname_mask = []; fname_atlas_res = []; fname_atlas_mask_res = [];
  % set output file names
  fname_regmat = sprintf('%s/%s_reg2atl.mat',parms.outdir,parms.fstem);
  fname_atlas_mask = sprintf('%s/%s_atlas_mask%s',parms.outdir,parms.fstem,parms.ext);
  fname_atlas_res = sprintf('%s/%s_atlas%s',parms.outdir,parms.fstem,parms.ext);
  fname_atlas_mask_res = sprintf('%s/%s_atlas_mask_res%s',parms.outdir,parms.fstem,parms.ext); 
  if ~exist(fname_regmat,'file')
    fprintf('%s: registering to atlas %s...\n',mfilename,fname_in);
  end
  % register to atlas and create brain mask
  p = [];
  p.fname_regmat = fname_regmat;
  p.fname_mask = fname_atlas_mask;
  p.fname_res = fname_atlas_res; 
  p.fname_mask_res = fname_atlas_mask_res; 
  p.hatl_offset_RAS = parms.hatl_offset_RAS;
  p.atlasdir = parms.atlasdir;
  p.atlasname = parms.atlasname;
  p.outdir = parms.outdir;
  p.outstem = parms.fstem;
  p.reg_atl_niters = parms.reg_atl_niters;
  p.forceflag = parms.forceflag;
  a = mmil_parms2args(p);
  M_atl_to_subj = mmil_reg_atlas(fname_in,a{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reg_failed,corr_atl,fvox_miss] = check_atlas_reg(M_atl_to_subj,fname_in,fname_mask_init,fname_res,parms,infix)
  reg_failed = 0; corr_atl = NaN; fvox_miss = NaN;
  if ~exist('infix','var'), infix = ''; end
  [fpath,fstem] = fileparts(fname_in);
  fname_check_reg = sprintf('%s/%s_reg2atl%s_check.mat',parms.outdir,fstem,infix);
  if ~exist(fname_check_reg,'file') || parms.forceflag
    fprintf('%s: checking atlas registration...\n',mfilename);
    fname_atl = sprintf('%s/%s.mgz',parms.atlasdir,parms.atlasname);
    if ~exist(fname_atl,'file')
      error('%s: %s file missing',mfilename,parms.fname_atl);
    end
    fname_atl_mask = sprintf('%s/%s_mask.mgz',parms.atlasdir,parms.atlasname);
    if ~exist(fname_atl_mask,'file')
       error('%s: %s file missing',mfilename,fname_atl_mask);
    end
    [vol_atl,M_atl,~,volsz_atl] = fs_load_mgh(fname_atl);
    vol_atl_mask = fs_load_mgh(fname_atl_mask);
    [vol_res,M_res] = fs_load_mgh(fname_res);
    % resample vol_res image to atlas space (dims and vox2ras may differ)
    vol_res_atl = mmil_resample_vol(vol_res,M_res,...
      'M_ref',M_atl,'nvox_ref',volsz_atl(1:3),...
      'M_reg',eye(4));
    idx_mask = find(vol_atl_mask>0);
    vals_res = vol_res_atl(idx_mask);
    vals_atl = vol_atl(idx_mask);
    % calculate correlation between resampled image and atlas
    %   in double precision: in single, the result depends on the number of
    %   MATLAB threads (the order of the sums changes the last digits)
    corr_atl = corr(double(vals_res),double(vals_atl));
    fprintf('%s: correlation of resampled brain with atlas = %0.3f\n',mfilename,corr_atl);
    % resample vol_res to native space (to check for brain cutoff)
    [vol_in,M_in,~,volsz_in] = fs_load_mgh(fname_in);
    vol_in_mask = fs_load_mgh(fname_mask_init);
    vol_res_in = mmil_resample_vol(vol_res,M_res,...
      'M_ref',M_in,'nvox_ref',volsz_in(1:3),...
      'M_reg',inv(M_atl_to_subj));
    idx_mask = find(vol_in_mask>0);
    vals_res = vol_res_in(idx_mask);
    vals_in = vol_in(idx_mask);
    % calculate fraction of missing brain voxels
    nvox_miss = abs(nnz(vals_res==0) - nnz(vals_in==0));
    nvox = length(idx_mask);
    fvox_miss = nvox_miss / nvox;
    fprintf('%s: fraction of resampled brain voxels missing = %0.3f\n',mfilename,fvox_miss);
    % decide whether reg to atlas failed
    if corr_atl < parms.min_corr_atl ||...
       fvox_miss > parms.max_fvox_miss
      reg_failed = 1;
    end
    save(fname_check_reg,'corr_atl','fvox_miss','reg_failed');
  else
    output = mmil_load_matfile(fname_check_reg,{'corr_atl','fvox_miss','reg_failed'});
    corr_atl = output.corr_atl;
    fvox_miss = output.fvox_miss;
    reg_failed = output.reg_failed;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,fname_mask_out,M_atl_to_subj] = resample_volumes(fname_in,fname_mask,fname_atlas_res,fname_atlas_mask_res,M_atl_to_subj,parms)
  fname_out = sprintf('%s/%s_res%s',parms.outdir,parms.fstem,parms.ext); 
  fname_mask_out = sprintf('%s/%s_res_mask%s',parms.outdir,parms.fstem,parms.ext);
 
  % update M_atl_to_subj
  if parms.atlas_flag
    if parms.native_flag
      [M_in,volsz_in] = mmil_load_mgh_info(fname_in,parms.forceflag);
      [M_res,volsz_res] = mmil_load_mgh_info(fname_atlas_res,parms.forceflag);
      M_subj_to_atl = center_volume(M_res,volsz_res,M_in,volsz_in);
      M_atl_to_subj = M_atl_to_subj * M_subj_to_atl;
    end
  else % atlas_flag = 0
    if parms.native_flag
      M_atl_to_subj = eye(4);
    else
      [M_in,volsz_in] = mmil_load_mgh_info(fname_in,parms.forceflag);
      % output volume parameters    
      M_out =  [-1     0     0   129;...
                 0     0     1  -129;...
                 0    -1     0   129;...
                 0     0     0     1];
      volsz_out = [256 256 256];
      % adjust offsets so brain does not get clipped
      M_atl_to_subj = center_volume(M_in,volsz_in,M_out,volsz_out);
    end
  end

  % save final outputs  
  if ~exist(fname_out,'file') || ~exist(fname_mask_out,'file') || parms.forceflag
    if parms.atlas_flag
      if ~parms.native_flag
        if ~isempty(parms.image_range)
          % load fname_atlas_res, clip values, save as fname_out
          [vol_res,M_res,mr_parms] = fs_load_mgh(fname_atlas_res);
          vol_res(vol_res<parms.image_range(1)) = parms.image_range(1);
          vol_res(vol_res>parms.image_range(2)) = parms.image_range(2);
          fs_save_mgh(vol_res,fname_out,M_res,mr_parms);
        else
          % copy fname_atlas_res to fname_out
          copyfile(fname_atlas_res,fname_out);
        end
        % copy fname_atlas_mask_res to fname_mask_out
        copyfile(fname_atlas_mask_res,fname_mask_out);
      else
        fprintf('%s: adjusting offsets to center volume...\n',mfilename);
        % load fname_in and fname_mask
        [vol_in,M_in,mr_parms,volsz_in] = fs_load_mgh(fname_in);
        vol_mask_in = fs_load_mgh(fname_mask);
        % resample to native ref space
        vol_res = mmil_resample_vol(vol_in,M_in,...
          'M_ref',M_in,'nvox_ref',volsz_in(1:3),'M_reg',M_atl_to_subj);
        % clip values
        if ~isempty(parms.image_range)
          vol_res(vol_res<parms.image_range(1)) = parms.image_range(1);
          vol_res(vol_res>parms.image_range(2)) = parms.image_range(2);
        end
        fs_save_mgh(vol_res,fname_out,M_in);
        % resamaple vol_mask to native ref space
        vol_mask_res = mmil_resample_vol(vol_mask_in,M_in,...
          'M_ref',M_in,'nvox_ref',volsz_in(1:3),'M_reg',M_atl_to_subj);
        vol_mask_res = 1.0*(vol_mask_res>0.5);
        fs_save_mgh(vol_mask_res,fname_mask_out,M_in);
      end
    else % atlas_flag = 0
      if parms.native_flag
        %% todo: is copying necessary? or should just change the output file names?
        % copy fname_in to fname_out
        copyfile(fname_in,fname_out);
        % copy fname_mask to fname_mask_out
        copyfile(fname_mask,fname_mask_out);
      else
        fprintf('%s: adjusting offsets to center volume...\n',mfilename);
        % load reference volume and mask
        [vol,M,mr_parms,volsz] = fs_load_mgh(fname_in);    
        vol_mask = fs_load_mgh(fname_mask);    
        % resample to 256^3, 1mm^3, coronal slices
        fprintf('%s: resampling %s to LIA...\n',mfilename,fname_in);
        vol = mmil_resample_vol(vol,M,...
          'M_ref',M_out,'nvox_ref',volsz_out(1:3),'M_reg',M_atl_to_subj);
        % clip values
        if ~isempty(parms.image_range)
          vol(vol<parms.image_range(1)) = parms.image_range(1);
          vol(vol>parms.image_range(2)) = parms.image_range(2);
        end
        fs_save_mgh(vol,fname_out,M_out,mr_parms);
        % resample mask
        fprintf('%s: resampling %s to LIA...\n',mfilename,fname_mask);
        vol_mask = mmil_resample_vol(vol_mask,M,...
          'M_ref',M_out,'nvox_ref',volsz_out(1:3),'M_reg',M_atl_to_subj);
        vol_mask = 1.0*(vol_mask>0.5);
        fs_save_mgh(vol_mask,fname_mask_out,M_out,mr_parms);        
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_reg = center_volume(M_in,volsz_in,M_out,volsz_out)
  M_reg = [];
  % transform vox2ras matrix to vox2lph
  M_in_lph = M_RAS_TO_LPH * M_in;
  M_out_lph = M_RAS_TO_LPH * M_out;
  % calculate offset needed to center volume
  C_in = M_in_lph(1:3,:)*[volsz_in(1:3)/2 1]';
  C_out = M_out_lph(1:3,:)*[volsz_out(1:3)/2 1]';
  T = C_in - C_out;
  % set M_reg using offset
  M_reg = Mtrans(T(1),T(2),T(3));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


