function M_volA_to_volB = epi_register_scans(fnameA,fnameB,varargin)
%function M_volA_to_volB = epi_register_scans(fnameA,fnameB,[options])
%
% Purpose: register one EPI scan to another
%
% Required Parameters:
%   fnameA: full path name of mgh/mgz file containing
%     reference 4D EPI volume
%   fnameB: full path name of mgh/mgz file containing 4D EPI volume
%     to be registered to reference
%
% Optional Parameters:
%   'fname_mask': input/output file name for brain mask of fnameA
%     If not empty and does not exist, will create one and save as fname_mask
%     If empty, will create one temporarily
%     {default=[]}
%   'fname_reg': output file name for registration matrix
%     {default=[]}
%   'forceflag': overwrite existing output files
%     {default: 0}
%
% Created:  02/27/10 by Don Hagler
% Prev Mod: 02/11/13 by Don Hagler
% Prev Mod: 07/16/20 by Don Hagler
% Prev Mod: 04/10/25 by Don Hagler
% Prev Mod: 04/23/25 by Don Hagler
% Last Mod: 06/25/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_volA_to_volB = [];
if ~mmil_check_nargs(nargin, 1), return; end
parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'fname_reg',[],[],...
  'forceflag',false,[false true],...
...
  'tmpdir',[],[],...
  'cleanupflag',true,[false true],...
... % mask dilation
  'mask_smooth1',10,[],...
  'mask_thresh1',0.05,[],...
  'mask_smooth2',10,[],...
  'mask_thresh2',0.5,[],...
... % parameters for rbreg_vol2vol_icorr_mask
  'mstep',4,[],...
  'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  'min_trans',0.05,[0,1],... % mm
  'min_rot',0.05,[0,1],... % degrees
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input files
if ~exist(fnameA,'file')
  error('file %s not found',fnameA);
end
if ~exist(fnameB,'file')
  error('file %s not found',fnameB);
end

% check for output file
if ~isempty(parms.fname_reg) && exist(parms.fname_reg,'file') &&...
   ~parms.forceflag
  load(parms.fname_reg);
end

if isempty(parms.tmpdir)
  if ~isempty(parms.fname_reg)
    [tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_reg);
  else
    [tmp_path,tmp_fstem,tmp_ext] = fileparts(fnameA);
  end    
  if isempty(tmp_path), tmp_path = pwd; end
  parms.tmpdir = [tmp_path '/tmp_EPI_reg'];
elseif exist(parms.tmpdir,'dir')
  parms.cleanupflag = false;
end

if isempty(parms.fname_mask)
  parms.fname_mask = [parms.tmpdir '/mask.mgz'];
end

mmil_mkdir(parms.tmpdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create brain mask from fnameA (if does not already exist)
cmig_create_brainmask(fnameA,...
  'fname_mask',parms.fname_mask,...
  'binarize_flag',1,...
  'verbose',1,...
  'forceflag',parms.forceflag);

% load volumes
volA = ctx_load_mgh(fnameA);
volB = ctx_load_mgh(fnameB);
volmask = ctx_load_mgh(parms.fname_mask);

% dilate mask
if parms.mask_smooth1>0
  volmask = mmil_dilate_mask(volmask,...
    'smooth1',parms.mask_smooth1,'thresh1',parms.mask_thresh1,...
    'smooth2',parms.mask_smooth2,'thresh2',parms.mask_thresh2);
end

% register volA and volB
M_volA_to_volB = rbreg_vol2vol_icorr_mask(volA,volB,volmask,...
  0,parms.mstep,0,parms.scales,...
  parms.min_trans,parms.min_rot*pi/180);

% save output file
if ~isempty(parms.fname_reg)
  save(parms.fname_reg,'M_volA_to_volB');
end

% remove tmp dir
if parms.cleanupflag
  cmd = sprintf('rm -r %s',parms.tmpdir);
  [status,result]=unix(cmd);
  if status
    fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n',mfilename,result);
  end
end


