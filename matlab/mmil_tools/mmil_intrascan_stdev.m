function results = mmil_intrascan_stdev(fnamelist,varargin)
%function results = mmil_intrascan_stdev(fnamelist,[options])
%
% Required input:
%   fnamelist: cell array of file names (mgh/mgz format)
%
% Optional input:
%   'outdir': output directory
%     {default = [pwd '/output]}
%   'outstem': output file stem
%     {default = 'data'}
%   'fname_brain': full path of input brain mask volume file
%     if not supplied, will create from mean of input images
%     {default = []}
%   'skipTRs': number of initial repetitions to omit from analysis
%     {default = 0}
%   'detrend': [0|1|2] whether and how to detrend input timeseries
%     0: no detrend (not recommended)
%     1: linear detrend
%     2: quadratic detrend
%     {default = 2}
%   'plotflag': [0|1] plot images of inter- and intrascan variability
%     {default = 0}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  08/14/2017 by Don Hagler
% Prev Mod: 11/06/2018 by Feng Xue
% Prev Mod: 08/07/2020 by Don Hagler
% Prev Mod: 03/14/2025 by Don Hagler
% Prev Mod: 04/23/2025 by Don Hagler
% Prev Mod: 06/12/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Last Mod: 10/20/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'outdir',[pwd '/output'],[],...
  'outstem','data',[],...
  'fname_brain',[],[],...
  'skipTRs',0,[0 Inf],...
  'detrend',2,[0:2],...
  'plotflag',false,[false true],...
  'verbose',true,[false true],...
  'forceflag',false,[false true],...
  ...
  'fig_size',[],[],...
  'tif_dpi',300,[10,10000],...
  'planestrings',{'HOR','SAG','COR'},[],...
  'slice_fract',[0.5,0.6,0.5],[0,1],...
  ...
  'plot_tags',{'outdir','outstem','border_flag','fig_size','tif_dpi',...
    'frame','planestrings','slice_fract','clim','visible_flag','forceflag',...
    'plim_flag','plim'},[],...
});

if ~iscell(fnamelist), fnamelist = {fnamelist}; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfiles = length(fnamelist);
keepflags = ones(nfiles,1);
for f=1:nfiles
  fname_in = fnamelist{f};
  % check that file exists
  if ~exist(fname_in,'file')
    if parms.verbose
      fprintf('%s: WARNING: skiping scan with file not found\n%s\n',mfilename,fname_in);
    end;
    keepflags(f) = 0;
    continue;
  end;
  % check how many frames
  [M,volsz,mr_parms] = mmil_load_mgh_info(fname_in,parms.forceflag,parms.outdir);
  if volsz(4)<=parms.skipTRs+3;
    if parms.verbose
      fprintf('%s: WARNING: skipping scan with too few frames (%d)\n%s\n',...
        mfilename,volsz(4),fname_in);
    end;
    keepflags(f) = 0;
    continue;
  end;
end;
fnamelist = fnamelist(find(keepflags));
nfiles = length(fnamelist);

if ~nfiles
  fprintf('%s: ERROR: all files in fnamelist are missing or have too few frames\n',mfilename);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

% calculate intrascan mean, stdev, and  tSNR for each volume timeseries
fnamelist_mean = {};
fnamelist_stdev = {};
fnamelist_tSNR = {};
for i=1:nfiles
  fname_in = fnamelist{i};
  [~,fstem,~] = fileparts(fname_in);
  fname_mean = sprintf('%s/%s_mean.mgz',parms.outdir,fstem);
  fname_stdev = sprintf('%s/%s_stdev.mgz',parms.outdir,fstem);
  % calculate mean and standard deviation
  vol_mean = []; vol_stdev = [];
  if ~exist(fname_mean,'file') || ~exist(fname_stdev,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating mean and stdev of %s...\n',mfilename,fname_in);
    end;
    % load data
    [vol,M] = fs_load_mgh(fname_in);
    % remove initial TRs and detrend
    p = [];
    p.skipTRs = parms.skipTRs;
    p.detrend = parms.detrend;
    p.norm_flag = false;
    p.mean_flag = true;
    a = mmil_parms2args(p);
    vol = mmil_detrend_vol(vol,a{:});
    % calculate mean and stdev
    vol_mean = mean(vol,4);
    vol_stdev = std(vol,[],4);
    fs_save_mgh(vol_mean,fname_mean,M);
    fs_save_mgh(vol_stdev,fname_stdev,M);
  end;
  % calculate tSNR
  fname_tSNR = sprintf('%s/%s_tSNR.mgz',parms.outdir,fstem);
  if ~exist(fname_tSNR,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating tSNR of %s...\n',mfilename,fname_in);
    end;
    % load mean and stdev
    if isempty(vol_mean), [vol_mean,M] = fs_load_mgh(fname_mean); end;
    if isempty(vol_stdev), [vol_stdev,M] = fs_load_mgh(fname_stdev); end;
    % save tSNR
    vol_tSNR = vol_mean ./ vol_stdev;
    fs_save_mgh(vol_tSNR,fname_tSNR,M);
  end;
  fnamelist_mean{i} = fname_mean;
  fnamelist_stdev{i} = fname_stdev;
  fnamelist_tSNR{i} = fname_tSNR;
  clear vol;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interscan mean (of intrascan means)
fname_inter_mean = sprintf('%s/%s_inter_mean.mgz',parms.outdir,parms.outstem);
% interscan stdev (of intrascan means)
fname_inter_stdev = sprintf('%s/%s_inter_stdev.mgz',parms.outdir,parms.outstem);
vol_inter_mean = []; vol_inter_stdev = [];
if ~exist(fname_inter_mean,'file') || ~exist(fname_inter_stdev,'file') || parms.forceflag
  vol_mean = [];
  for i=1:nfiles
    fname_mean = fnamelist_mean{i};
    if parms.verbose
      fprintf('%s: loading %s...\n',mfilename,fname_mean);
    end;
    [vol,M] = fs_load_mgh(fname_mean);
    vol_mean = cat(4,vol_mean,vol);
  end;
  if parms.verbose
    fprintf('%s: calculating interscan mean of intrascan means...\n',mfilename);
  end;
  vol_inter_mean = mean(vol_mean,4);
  fs_save_mgh(vol_inter_mean,fname_inter_mean,M);
  if parms.verbose
    fprintf('%s: calculating interscan stdev of intrascan means...\n',mfilename);
  end;
  vol_inter_stdev = std(vol_mean,[],4);
  fs_save_mgh(vol_inter_stdev,fname_inter_stdev,M);
end;

% interscan tSNR
fname_inter_tSNR = sprintf('%s/%s_inter_tSNR.mgz',parms.outdir,parms.outstem);
if ~exist(fname_inter_tSNR,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: calculating interscan tSNR of %s...\n',mfilename,fname_in);
  end;
  % load mean and stdev
  if isempty(vol_inter_mean), [vol_inter_mean,M] = fs_load_mgh(fname_inter_mean); end;
  if isempty(vol_inter_stdev), [vol_inter_stdev,M] = fs_load_mgh(fname_inter_stdev); end;
  % save tSNR
  vol_inter_tSNR = vol_inter_mean ./ vol_inter_stdev;
  fs_save_mgh(vol_inter_tSNR,fname_inter_tSNR,M);
end;

% calculate interscan average of intrascan standard deviation
fname_intra_stdev = sprintf('%s/%s_intra_stdev.mgz',parms.outdir,parms.outstem);
if ~exist(fname_intra_stdev,'file') || parms.forceflag
  vol_stdev = [];
  for i=1:nfiles
    fname_stdev = fnamelist_stdev{i};
    if parms.verbose
      fprintf('%s: loading %s...\n',mfilename,fname_stdev);
    end;
    [vol,M] = fs_load_mgh(fname_stdev);
    vol_stdev = cat(4,vol_stdev,vol);
  end;
  if parms.verbose
    fprintf('%s: calculating mean of intrascan stdev...\n',mfilename);
  end;
  vol_intra_stdev = mean(vol_stdev,4);
  fs_save_mgh(vol_intra_stdev,fname_intra_stdev,M);
end;

% calculate interscan average of intrascan tSNR
fname_intra_tSNR = sprintf('%s/%s_intra_tSNR.mgz',parms.outdir,parms.outstem);
if ~exist(fname_intra_tSNR,'file') || parms.forceflag
  vol_tSNR = [];
  for i=1:nfiles
    fname_tSNR = fnamelist_tSNR{i};
    if parms.verbose
      fprintf('%s: loading %s...\n',mfilename,fname_tSNR);
    end;
    [vol,M] = fs_load_mgh(fname_tSNR);
    vol_tSNR = cat(4,vol_tSNR,vol);
  end;
  if parms.verbose
    fprintf('%s: calculating mean of intrascan tSNR...\n',mfilename);
  end;
  vol_intra_tSNR = mean(vol_tSNR,4);
  fs_save_mgh(vol_intra_tSNR,fname_intra_tSNR,M);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check fname_brain
if ~isempty(parms.fname_brain)
  fname_mask = parms.fname_brain;
else
  % create brain mask from mean image
  [~,fstem,fext] = fileparts(fname_inter_mean);
  fname_mask = sprintf('%s/%s_mask%s',parms.outdir,fstem,fext);
  cmig_create_brainmask(fname_inter_mean,...
    'fname_mask',fname_mask,...
    'outdir',parms.outdir,'binarize_flag',1,...
    'verbose',parms.verbose,'forceflag',parms.forceflag);
end;

% summarize results in a mat file
% NOTE: a companion _results.json was written here until 08/2026, by jsonlab's
%   savejson. Nothing read it: callers use the returned struct, and this mat
%   file is what the cache check below loads
fname_mat = sprintf('%s/%s_results.mat',parms.outdir,parms.outstem);
calc_flag = ~exist(fname_mat,'file') || parms.forceflag;
if ~calc_flag
  load(fname_mat);
  % recalculate a results file written before the median fields were added
  calc_flag = ~isfield(results,'inter_tSNR_median');
end;
if calc_flag
  results = [];
  % load mask
  if parms.verbose
    fprintf('%s: loading %s...\n',mfilename,fname_mask);
  end;
  vol_mask = fs_load_mgh(fname_mask);
  ind_mask = find(vol_mask>0);

  % calculate average interscan mean within brain mask
  if parms.verbose
    fprintf('%s: calculating average interscan mean within brain mask...\n',mfilename);
  end;
  results.inter_mean = calc_brain_stats(fname_inter_mean,ind_mask);
   
  % calculate average interscan stdev within brain mask
  if parms.verbose
    fprintf('%s: calculating average interscan stdev within brain mask...\n',mfilename);
  end;
  results.inter_stdev = calc_brain_stats(fname_inter_stdev,ind_mask);

  % calculate average intrascan stdev value within brain mask
  if parms.verbose
    fprintf('%s: calculating average intrascan stdev within brain mask...\n',mfilename);
  end;
  results.intra_stdev = calc_brain_stats(fname_intra_stdev,ind_mask);

  % calculate average interscan tSNR within brain mask
  if parms.verbose
    fprintf('%s: calculating average interscan tSNR within brain mask...\n',mfilename);
  end;
  % finite voxels only: a voxel whose intrascan means are equal across scans
  %   has interscan stdev 0 and tSNR Inf (or NaN for 0/0)
  [results.inter_tSNR,results.inter_tSNR_median] = ...
    calc_brain_stats(fname_inter_tSNR,ind_mask,1);
  
  % calculate average intrascan tSNR within brain mask
  if parms.verbose
    fprintf('%s: calculating average intrascan tSNR within brain mask...\n',mfilename);
  end;
  [results.intra_tSNR,results.intra_tSNR_median] = ...
    calc_brain_stats(fname_intra_tSNR,ind_mask);
  
  % calculate coefficient of variation
  results.inter_cov = results.inter_stdev / results.inter_mean;
  results.intra_cov = results.intra_stdev / results.inter_mean;

  % calculate ratio between interscan and intrascan stdev
  results.inter_intra_stdev_ratio = results.inter_stdev / results.intra_stdev;

  % calculate ratio between interscan and intrascan tSNR
  results.inter_intra_tSNR_ratio = results.inter_tSNR / results.intra_tSNR;
  results.inter_intra_tSNR_median_ratio = ...
    results.inter_tSNR_median / results.intra_tSNR_median;

  save(fname_mat,'results');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot images for each output file
if parms.plotflag

  % plot interscan mean
  outstem = [parms.outstem '_inter_mean'];
  plot_images(fname_inter_mean,outstem,'interscan mean',parms);
  
  % plot mask
  outstem = [parms.outstem '_mask'];
  if ~isempty(parms.fname_brain)
    desc = 'brain mask supplied';
  else
    desc = 'mask derived from interscan mean';
  end;
  plot_images(fname_mask,outstem,desc,parms);

  % plot interscan stdev (calculated across average scans)
  outstem = [parms.outstem '_inter_stdev'];
  plot_images(fname_inter_stdev,outstem,'interscan stdev',parms);

  % plot intrascan stdev (averaged across scans)
  outstem = [parms.outstem '_intra_stdev'];
  plot_images(fname_intra_stdev,outstem,'intrascan stdev',parms);

  % plot interscan tSNR (calculated across average scans)
  outstem = [parms.outstem '_inter_tSNR'];
  plot_images(fname_inter_tSNR,outstem,'interscan tSNR',parms);

  % plot intrascan tSNR (averaged across scans)
  outstem = [parms.outstem '_intra_tSNR'];
  plot_images(fname_intra_tSNR,outstem,'intrascan tSNR',parms);
  
  % plot intrascan mean, stdev, and tSNR for individual scans
  for i=1:nfiles
    % plot intrascan mean
    fname_in = fnamelist_mean{i};
    [~,outstem,~] = fileparts(fname_in);
    desc = sprintf('intrascan mean from %s',outstem);
    plot_images(fname_in,outstem,desc,parms);
    
    % plot intrascan stdev
    fname_in = fnamelist_stdev{i};
    [~,outstem,~] = fileparts(fname_in);
    desc = sprintf('intrascan stdev from %s',outstem);
    plot_images(fname_in,outstem,desc,parms);

    % plot intrascan tSNR
    fname_in = fnamelist_tSNR{i};
    [~,outstem,~] = fileparts(fname_in);
    desc = sprintf('intrascan tSNR from %s',outstem);
    plot_images(fname_in,outstem,desc,parms);
  end;  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanval,medval,stdval] = calc_brain_stats(fname,ind_mask,finite_flag)
  if ~exist('finite_flag','var') || isempty(finite_flag), finite_flag = 0; end;
  [vol,~] = fs_load_mgh(fname,[],1);
  vals = vol(ind_mask);
  % optionally exclude Inf and NaN values
  if finite_flag, vals = vals(isfinite(vals)); end;
  meanval = mean(vals);
  if nargout>1, medval = median(vals); end;
  if nargout>2, stdval = std(vals); end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_images(fname,outstem,desc,parms)
  parms.outstem = outstem;
  fname_out = sprintf('%s/%s.jpg',parms.outdir,outstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: plotting %s...\n',mfilename,desc);
    end;
    args = mmil_parms2args(parms,parms.plot_tags);
    mmil_subplots(fname,args{:});
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
