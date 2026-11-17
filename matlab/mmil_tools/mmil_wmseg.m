function [fname_wm,fname_bm,fname_gm] = mmil_wmseg(fname_in,varargin)
%function [fname_wm,fname_bm,fname_gm] = mmil_wmseg(fname_in,varargin)
%
% Purpose: create white matter segmentation volume from T1-weighted image
%
% Usage: mmil_wmseg(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_wm': output file name for white matter mask
%    {default = 'wm.mgz'}
%  'fname_bm': output file name for brain mask
%    {default = 'bm.mgz'}
%  'fname_gm': output file name for gray matter mask
%    {default = 'gm.mgz'}
%  'erode_wm_flag': [0|1] erode white matter mask
%    {default = 0}
%  'dilate_bm_flag': [0|1] dilate brain mask
%    {default = 0}
%  'verbose': [0|1] display progress messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  07/13/17 by Don Hagler
% Prev Mod: 08/07/17 by Don Hagler
% Prev Mod: 09/10/20 by Don Hagler
% Last Mod: 09/16/20 by Don Hagler
%

% Based on corSeg_amd
% Last Mod:  07/07/17 by Anders Dale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_wm = []; fname_bm = []; fname_gm = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);

% check output file
if exist(parms.fname_wm,'file') &&...
   exist(parms.fname_bm,'file') && ~parms.forceflag, return; end;

% load data
vol_brain = load_data(parms);

% bias correction (N3) for head
vol_brain = bias_corr(parms,vol_brain,'head');

% perform skull stripping
[vol_brain,seg_brain] = skull_strip(parms,vol_brain);

% bias correction (N3) for brain
vol_brain = bias_corr(parms,vol_brain,'brain');

% atlas registration and brain matter segmentation
[vol_wm,vol_gm,seg_brain] = brain_segment(parms,vol_brain,seg_brain);

% optionally erode white matter mask
if parms.erode_wm_flag
  vol_wm = erode_wm(parms,vol_wm);
end;

% optionally dilate brain mask
if parms.dilate_bm_flag
  vol_brain = dilate_bm(parms,vol_brain);
end;

% save output
ctx_save_mgh(vol_wm,parms.fname_wm);
ctx_save_mgh(vol_brain,parms.fname_bm);
ctx_save_mgh(vol_gm,parms.fname_gm);

fname_wm = parms.fname_wm;
fname_bm = parms.fname_bm;
fname_gm = parms.fname_gm;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_wm','wm.mgz',[],...
    'fname_bm','bm.mgz',[],...
    'fname_gm','gm.mgz',[],...
    'erode_wm_flag',false,[false true],...
    'dilate_bm_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % for wm segmentation
    'wm_codes',[1,9,17,19,21],[],... % 1 = cortical wm, 9 = cerebellar wm, 17 = ventral DC, 19 = brainstem, 21 = wm hypointensity
  ... % for gm segmentation
    'gm_codes',[2,11,12,14,15],[],...
  ... % for N3 bias correction
    'estop',0.001,[],...
    'maxiter',100,[],...
  ... % for erosion of white matter mask
    'nvoxels',1,[1:100],...
    'sigma',1,[1e-2,10],...
    'thresh_init',1e-5,[0,Inf],...
    'thresh',0.99,[0.1,0.999],...
  ... % for dilation of brain mask
    'smooth1',30,[],...
    'thresh1',0.05,[],...
    'smooth2',30,[],...
    'thresh2',0.05,[],...
    'smooth3',0,[],...
    'thresh3',0,[],...  
  ...
    'dilate_tags',{'thresh0','smooth1','thresh1',...
                   'smooth2','thresh2','smooth3','thresh3',...
                   'erode_flag'},[],... % 'fname_in','fname_out','forceflag'
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_brain = load_data(parms)
  if parms.verbose
    fprintf('%s: loading data from %s...\n',mfilename,parms.fname_in);
  end;
  vol_brain = ctx_load_mgh(parms.fname_in);
  vol_brain.minI = 0;
  vol_brain.maxI = max(vol_brain.imgs(:));
  vol_brain.slicegap = 0;
  vol_brain.sf = 1.0;
  if vol_brain.maxI > 4096
    vol_brain.imgs = vol_brain.imgs*4096/(vol_brain.maxI-vol_brain.minI)+vol_brain.minI;
    vol_brain.maxI=4096;
    vol_brain.minI=0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_brain = bias_corr(parms,vol_brain,label)
  if parms.verbose
    fprintf('%s: bias correction for %s...\n',mfilename,label);
  end;
  vol_brain = vol_correct_bias_field_n3(vol_brain,parms.estop,parms.maxiter);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol_brain,seg_brain] = skull_strip(parms,vol_brain)
  if parms.verbose
    fprintf('%s: creating brain mask...\n',mfilename);
  end;
  %% todo: use mmil_dct_brainmask?
  sampling =[4 4 4];
  nK = [5 5 5];
  boutput = false;
  sf = 1;
  [bm_gen, seg_brain.M_hatl_to_head_rb, regStruct] = dctMorph_SkullStrip(vol_brain, sampling, nK, boutput, sf, false);
  seg_brain.brainmesh = SkullStripAMDwrapper(vol_brain, regStruct);
  seg_brain.MI = regStruct.min_cost_rb;
  vol_brain.imgs = vol_brain.imgs *regStruct.sfm;
  vol_brain.maxI = 600;
  vol_brain.minI= 0;
  clear regStruct volm;
  vol_brain = getmaskvol(vol_brain, seg_brain.brainmesh, eye(4,4));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol_wm,vol_gm,seg_brain] = brain_segment(parms,vol_brain,seg_brain)
  if parms.verbose
    fprintf('%s: registration to atlas...\n',mfilename);
  end;
  [regStruct, M_atl_to_vol_rb, volm] = brain_reg2atl(vol_brain,...
    seg_brain.M_hatl_to_head_rb,false);
  seg_brain.sfm = regStruct.sfm;
  vxlmap = getDCTVxlMapping(vol_brain, regStruct, volm);
  clear volm;
  if parms.verbose
    fprintf('%s: brain segmentation...\n',mfilename);
  end;
  vol_brain.imgs = vol_brain.imgs * seg_brain.sfm;
  vol_brain.maxI = 600;
  vol_brain.minI = 0;
  load('ctxatl_1_4_rbm.mat');
  seg_brain.topLabel = ctxSeg(vol_brain,ctxatl,vxlmap,3,5);
  seg_brain.Left_cm3 = seg_brain.topLabel.Left_cm3;
  seg_brain.Right_cm3 = seg_brain.topLabel.Right_cm3;
  clear ctxatl vxlmap;
  seg_brain.topLabel.imgs = double(seg_brain.topLabel.imgs);
  seg_brain.vol2 = vol_brain;
  seg_brain.brainmask = seg_brain.vol2;
  seg_brain.icv = length(find((seg_brain.topLabel.imgs > 0) &...
                              (seg_brain.topLabel.imgs~=22)))/1000;
  seg_brain.M_atl_to_vols = eye(4,4);
  seg_brain.anatomical_atlas.NumClasses = 23;
  seg_brain.feature_atlas.NumClasses = 23;
  seg_brain.numSpectra = 1;
  seg_brain.pvFactor = -1;
  seg_brain.pvList = zeros(0,5);
  vol_wm = vol_brain;
  vol_wm.imgs = ismember(seg_brain.topLabel.imgs,parms.wm_codes);
  vol_gm = vol_brain;
  vol_gm.imgs = ismember(seg_brain.topLabel.imgs,parms.gm_codes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_wm = erode_wm(parms,vol_wm)
  if parms.verbose
    fprintf('%s: eroding white matter mask...\n',mfilename);
  end;
  vol = 1.0 * (vol_wm.imgs>1-parms.thresh_init);
  for j=1:parms.nvoxels
    vol = mmil_smooth3d(vol,parms.sigma,parms.sigma,parms.sigma);
    vol = 1.0*(vol>=parms.thresh);
  end;
  vol_wm.imgs = vol;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_bm = dilate_bm(parms,vol_bm)
  if parms.verbose
    fprintf('%s: dilating brain mask...\n',mfilename);
  end;
  args = mmil_parms2args(parms,parms.dilate_tags);
  vol_bm = mmil_dilate_mask(vol_bm,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

