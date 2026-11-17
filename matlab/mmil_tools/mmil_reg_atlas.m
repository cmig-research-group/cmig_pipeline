function M_atl_to_subj = mmil_reg_atlas(fname,varargin)
%function M_atl_to_subj = mmil_reg_atlas(fname,varargin)
%
% Required parameters:
%  fname: input file name
%   T1w volume as mgh/mgz file expected
%
% Optional parameters:
%  'outdir': output directory
%    if empty, set to pwd
%    {default = []}
%  'outstem': output file stem
%    if empty, set to file stem of input file name
%    {default = []}
%  'fname_regmat': output file name for registration mat file
%    if empty, based on outdir and outstem
%    {default = []}
%  'fname_mask': output file name for brain mask
%    if empty, based on outdir and outstem
%    {default = []}
%  'fname_res': output file name for resampled volume
%    if empty, based on outdir and outstem
%    {default = []}
%  'fname_mask_res': output file name for atlas-registered mask
%    if empty, based on outdir and outstem
%    {default = []}
%  'atlasdir': full path of atlas directory
%    {default =  cmig_getenv('CMIG_PIPELINE_ATLASDIR')}}
%  'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%    {default = 'T1_Atlas/T1_atlas'}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%    {default = '.mgz'}
%  'forceflag': [0|1] whether to overwrite existing output
%    {default = 0}
%
% Created:  07/24/20 by Don Hagler
% Prev Mod: 09/30/20 by Don Hagler
% Prev Mod: 10/19/20 by Don Hagler
% Prev Mod: 10/01/24 by Don Hagler
% Prev Mod: 04/04/25 by Anders Dale
% Prev Mod: 04/11/25 by Don Hagler
% Prev Mod: 05/14/25 by Don Hagler
% Prev Mod: 08/13/25 by Don Hagler
% Last Mod: 10/26/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: based on a portion of
%   MMIL_Register_and_Resample_MRI_Volumes
% Created:  03/18/08 by Don Hagler
% Prev Mod: 07/14/17 by Don Hagler
% Prev Mod: 08/06/17 by Don Hagler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_atl_to_subj = [];
if ~mmil_check_nargs(nargin,1), return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin, { ...
  'outdir',[],[],...
  'outstem',[],[],...
  'atlasdir',[],[],...
  'atlasname','T1_Atlas/T1_atlas',[],...
  'preregflag',true,[false true],...
  'cleanupflag',true,[false true],...
  'forceflag',false,[false true],...
...
  'fname_regmat',[],[],...
  'fname_mask',[],[],...
  'fname_res',[],[],...
  'fname_mask_res',[],[],...
...
  'reg_atl_niters',2,[1,10],...
  'hatl_offset_RAS',[0 0 0],[],...
  'mask_smooth',15,[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

% check fname_atl
if isempty(parms.atlasdir)
  parms.atlasdir = cmig_getenv('CMIG_PIPELINE_ATLASDIR');
end;
parms.fname_atl = sprintf('%s/%s.mgz',parms.atlasdir,parms.atlasname);
if ~exist(parms.fname_atl,'file')
  error('%s: %s file missing',mfilename,parms.fname_atl);
end;
parms.fname_atl_mask = sprintf('%s/%s_mask.mgz',parms.atlasdir,parms.atlasname);
if ~exist(parms.fname_atl_mask,'file')
   error('%s: %s file missing',mfilename,parms.fname_atl_mask);
end;
parms.fname_atl_maskBroad = sprintf('%s/%s_maskBroad.mgz',parms.atlasdir,parms.atlasname);
if ~exist(parms.fname_atl_maskBroad,'file')
  error('%s: %s file missing',mfilename,parms.fname_atl_maskBroad);
end;

% construct offset matrix
M_hatl_offset = mmil_construct_M('trans',parms.hatl_offset_RAS);

if isempty(parms.outdir), parms.outdir = pwd; end;
mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% register to atlas and create brain mask
[fpath,fstem] = fileparts(fname);
if isempty(parms.outstem), parms.outstem = fstem; end;

% set output file names if not supplied
if isempty(parms.fname_regmat)
  parms.fname_regmat = sprintf('%s/%s_reg2atl.mat',parms.outdir,parms.outstem);
end;
if isempty(parms.fname_mask)
  parms.fname_mask = sprintf('%s/%s_mask%s',parms.outdir,parms.outstem,parms.ext);
end;
if isempty(parms.fname_res)
  parms.fname_res = sprintf('%s/%s_res%s',parms.outdir,parms.outstem,parms.ext); 
end;
if isempty(parms.fname_mask_res)
  parms.fname_mask_res = sprintf('%s/%s_res_mask%s',parms.outdir,parms.outstem,parms.ext);
end;

if ~exist(parms.fname_res,'file') || ~exist(parms.fname_mask_res,'file') || parms.forceflag
  if ~exist(parms.fname_mask,'file') || ~exist(parms.fname_regmat,'file') ||...
      ~exist(parms.fname_mask_res,'file') || parms.forceflag
    % pre-registration to atlas
    if parms.preregflag
      % initial registration to atlas
      fprintf('%s: pre-registration to atlas of %s...\n',...
        mfilename,fname);
      fname_prereg = sprintf('%s/%s_prereg%s',...
        parms.outdir,parms.outstem,parms.ext); 
      M_atl_to_subj_prereg = initial_rigidreg(fname,fname_prereg,parms);
      % registration to atlas to create brain mask
      fprintf('%s: nonlinear registration to atlas of %s...\n',...
        mfilename,fname_prereg);
      fname_mask_prereg = sprintf('%s/%s_prereg_mask%s',...
        parms.outdir,parms.outstem,parms.ext);
    else
      M_atl_to_subj_prereg = eye(4);
      fname_prereg = fname;
      fname_mask_prereg = sprintf('%s/%s_mask%s',...
        parms.outdir,parms.outstem,parms.ext);
    end
    [volmask_prereg,M_atl_to_subj_postreg] = mmil_dct_brainmask([],...
      'fname_in',fname_prereg,...
      'fname_mask',fname_mask_prereg,...
      'atlasdir',parms.atlasdir,'atlasname',parms.atlasname,...
      'smooth',parms.mask_smooth,'forceflag',1);
    [volmask_prereg,M_prereg] = ctx_ctx2mgh(volmask_prereg);
    % combine M_atl_to_subj_prereg and M_atl_to_subj_postreg
    M_atl_to_subj = M_atl_to_subj_prereg * M_atl_to_subj_postreg;
    % offset output location
    M_atl_to_subj = M_atl_to_subj*M_hatl_offset;  
    % resample ref to res
    fprintf('%s: resampling volume...\n',mfilename);
    [M_res,volsz_res] = mmil_load_mgh_info(fname_prereg,parms.forceflag);
    [vol_ref,M_ref,mr_parms,volsz_ref] = fs_load_mgh(fname);
    vol_ref_res = mmil_resample_vol(vol_ref,M_ref,...
      'M_ref',M_res,'nvox_ref',volsz_res(1:3),...
      'M_reg',M_atl_to_subj);
    fs_save_mgh(vol_ref_res,parms.fname_res,M_res);
    % resample volmask_prereg to volmask_res
    fprintf('%s: resampling mask volumes...\n',mfilename);
    volmask_res = mmil_resample_vol(volmask_prereg,M_prereg,...
      'M_ref',M_res,'nvox_ref',volsz_res(1:3),...
      'M_reg',M_atl_to_subj_postreg);
    volmask_res = 1.0*(volmask_res>0.5);
    fs_save_mgh(volmask_res,parms.fname_mask_res,M_res);
    % resample mask back to reference volume
    volmask_ref = mmil_resample_vol(volmask_res,M_res,...
      'M_ref',M_ref,'nvox_ref',volsz_ref(1:3),...
      'M_reg',inv(M_atl_to_subj));
    volmask_ref = 1.0*(volmask_ref>0.5);
    fs_save_mgh(volmask_ref,parms.fname_mask,M_ref);
    clear volmask volmask_res volmask_ref;    
    fname_mask = parms.fname_mask;
    save(parms.fname_regmat,'M_atl_to_subj','fname','fname_mask');
  end
else
  load(parms.fname_regmat);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_atl_to_subj_init = initial_rigidreg(fname,fname_out,parms)
  % load volumes
  vol_t1_ctx = ctx_load_mgh(fname);
  volm_ctx = ctx_load_mgh(parms.fname_atl);
  mask_ctx = ctx_load_mgh(parms.fname_atl_mask);
  
  % smooth mask
  mask_dilated = real(mmil_smooth3d(bwdist(mask_ctx.imgs>0.5)<40,20,20,20));
  mask_dilated_ctx = mask_ctx; 
  mask_dilated_ctx.imgs = mask_dilated;

  % register volumes
  tic
  fprintf('%s: registering volumes...\n',mfilename);
  [dummy regStruct] = register_volumes_affine(vol_t1_ctx, volm_ctx, mask_dilated, 'icorr', parms.reg_atl_niters); % note: does not work well with a tight mask
  toc

  % save registered and resampled volume
  M_atl_to_subj_init = regStruct.M_atl_to_vol_rb;
  M_LIA_to_orig = mmil_resample_to_LIA(fname,fname_out,'M_reg',M_atl_to_subj_init);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
