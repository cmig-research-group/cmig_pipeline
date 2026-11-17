function fname_mask = mmil_brainmask(fname_in,varargin)
%function fname_mask = mmil_brainmask(fname_in,[options])
%
% Purpose: create brain mask volume from T1-weighted image
%
% Usage: mmil_brainmask(fname_in,'key1',value1,...)
%
% Required Parameters:
%  fname_in: input T1-weighted file name
%
% Optional Parameters:
%  'fname_mask': output file name
%    if not supplied, will generate a name based on file stem of fname and outdir
%    {default = []}
%  'outdir': output directory
%    if not supplied, will use path of fname_mask
%    will use pwd if both outdir and fname_mask not supplied
%    {default = []}
%  'n3_flag': [0|1] perform N3 intensity correction
%    {default = 1}
%  'binarize_flag': [0|1] output binary mask
%    otherwise return masked brain image
%    {default = 0}
%  'verbose': [0|1] display progress messages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  07/04/25 by Don Hagler
% Last Mod: 07/04/25 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on mmil_wmseg
%   Created:  07/13/17 by Don Hagler
%   Last Mod: 09/16/20 by Don Hagler
% Based on corSeg_amd
%   Last Mod: 07/07/17 by Anders Dale
% Based on mmil_fsl_brainmask
%   Created:  02/20/16 by Don Hagler
%   Last Mod: 01/31/22 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mask = []; fname_brain = [];
if ~mmil_check_nargs(nargin,1), return; end

% check input parameters
parms = check_input(fname_in,varargin);

% quit if output already exists
if exist(parms.fname_mask,'file') && ~parms.forceflag
  return;
end

% load data
vol_brain = load_data(parms);

% bias correction (N3) for head
if parms.n3_flag
  vol_brain = bias_corr(parms,vol_brain,'head');
end

% perform skull stripping
[vol_brain,seg_brain] = skull_strip(parms,vol_brain);

% bias correction (N3) for brain
if parms.n3_flag
  vol_brain = bias_corr(parms,vol_brain,'brain');
end

% create vol_mask from vol_brain
if parms.binarize_flag
  vol_brain = get_mask(vol_brain);
end

% save output
ctx_save_mgh(vol_brain,parms.fname_mask);

fname_mask = parms.fname_mask;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_mask','mask.mgz',[],...
    'outdir',[],[],...
    'n3_flag',true,[false true],...
    'binarize_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % for N3 bias correction
    'estop',0.001,[],...
    'maxiter',100,[],...
  ...
    'dilate_tags',{'thresh0','smooth1','thresh1',...
                   'smooth2','thresh2','smooth3','thresh3',...
                   'erode_flag'},[],... % 'fname_in','fname_out','forceflag'
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end

  if isempty(parms.outdir)
    if isempty(parms.fname_mask)
      parms.outdir = pwd;
    else
      parms.outdir = fileparts(parms.fname_mask);
      if isempty(parms.outdir), parms.outdir = pwd; end
    end
  end

  [fpath,fstem,fext] = fileparts(parms.fname_in);
  suffix_brain = '_brainmask';
  if parms.binarize_flag
    suffix = '_brainmask_mask';
  else
    suffix = suffix_brain;
  end

  % set output file name
  if isempty(parms.fname_mask)
    parms.fname_mask = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,fext);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_brain = load_data(parms)
  if parms.verbose
    fprintf('%s: loading data from %s...\n',mfilename,parms.fname_in);
  end
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
  end
  vol_brain = vol_correct_bias_field_n3(vol_brain,parms.estop,parms.maxiter);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol_brain,seg_brain] = skull_strip(parms,vol_brain)
  if parms.verbose
    fprintf('%s: creating brain mask...\n',mfilename);
  end
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

function vol_brain = get_mask(vol_brain)
  vol_brain.imgs = 1.0*(vol_brain.imgs>0);
  vol_brain.maxI = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


