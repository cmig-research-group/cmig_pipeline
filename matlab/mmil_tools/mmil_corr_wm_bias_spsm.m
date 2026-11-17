function [vol_corr,vol_ratio,vol_pred,vol_wm_masked] = ...
  mmil_corr_wm_bias_spsm(vol,vol_wm,vol_bm,vol_gm,varargin)
%function [vol_corr,vol_ratio,vol_pred,vol_wm_masked] = ...
%  mmil_corr_wm_bias_spsm(vol,vol_wm,vol_bm,vol_gm,[options])
%
% Purpose: estimate B1 bias field using white matter segmentation mask
%   and sparse smoothing in volume
%
% Required Input:
%   vol: brain image volume
%   vol_wm: white matter mask volume
%   vol_bm: brain mask volume
%   vol_gm: gray matter mask volume
%
% Optional Input ('key',value pairs):
%   'vol_target': target volume for white matter intensity
%     may be a scalar value
%     if empty, will use median white matter intensity value
%     {default = []}
%   'ratio_range': range of values to constrain ratio
%     relative to median ratio within brain mask
%     if empty, no clipping
%     {default = []}
%   'image_range': range of values to clip output image
%     if empty, no clipping
%     {default = []}
%   'outlier_thresh': exclude outliers from wm mask
%     difference from local max/min relative to mean wm intensity
%     {default = 0.1}
%   'outlier_edge_flag': [0|1] only exclude outlier voxels from
%     edge of white matter mask
%     {default = 1}
%   'lambda1': smoothing factor inside mask
%     {default = 20}
%   'lambda2': smoothing factor outside mask
%     {default = 1000}
%   'postsmooth_flag': [0|1] whether to apply post-smoothing
%       to the bias field ratio
%     {default = 1}
%
% Output:
%   vol_corr: corrected volume
%   vol_ratio: bias field ratio volume
%   vol_pre: predicted intensity volume
%   vol_wm_masked: white matter mask excluding outliers
%
% Created:  02/23/15 by Don Hagler
% Prev Mod: 11/28/15 by Don Hagler
% Prev Mod: 11/06/17 by Anders Dale
% Prev Mod: 11/06/17 by Don Hagler
% Last Mod: 10/07/20 by Don Hagler
%
% Based on code by Anders Dale
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_corr = []; vol_ratio = []; vol_pred = []; vol_wm_masked = [];
if ~mmil_check_nargs(nargin,3), return; end;
parms = mmil_args2parms(varargin,{...
  'vol_target',[],[],...
  'ratio_range',[],[],...
  'image_range',[],[],...
  'outlier_thresh',0.1,[0,1],...
  'outlier_edge_flag',true,[false true],...
  'lambda1',20,[1e-10,1e10],...
  'lambda2',1e3,[1e-10,1e10],...
  'postsmooth_flag',true,[false true],...
});

vol_wm = 1.0*(vol_wm>0);
vol_bm = 1.0*(vol_bm>0);
vol_gm = 1.0*(vol_gm>0);

if isempty(parms.vol_target)
  parms.vol_target = median(vol(vol_wm>0));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[i1max i1min] = maxmin(find(sum(sum(vol_bm,3),2)));
[i2max i2min] = maxmin(find(sum(sum(vol_bm,3),1)));
[i3max i3min] = maxmin(find(sum(sum(vol_bm,2),1)));

% compare wm and gm values to determine contrast type
meanval_wm = mean(vol(find(vol_wm>0)));
meanval_gm = mean(vol(find(vol_gm>0)));
% if wm brighter than gm, then T1w images
T1flag = (meanval_wm-meanval_gm > 0);
% if wm darker than gm, then T2w images

if parms.outlier_edge_flag
  % set nvoxels depending on number of voxels in vol
  nvoxels = max(1,ceil(size(vol_wm,1)/100));
  % erode white matter mask
  vol_wm_eroded = erode_mask(vol_wm,nvoxels);
end;

% remove voxels from mask if intensity
%   is much less than local maximnum within mask
%   (or much more than local minimum if T2w)
vol_val_wm = vol.*vol_wm;       
vol_val_locmax = vol_val_wm;
vol_exclude = false;
niter = 2;
for iter = 1:niter
  vol_tmp = cat(4,...
    circshift(vol_val_locmax,[-1,0,0]),circshift(vol_val_locmax,[1,0,0]),...
    circshift(vol_val_locmax,[0,-1,0]),circshift(vol_val_locmax,[0,1,0]),...
    circshift(vol_val_locmax,[0,0,-1]),circshift(vol_val_locmax,[0,0,1]));
  if T1flag
    vol_val_locmax = max(vol_tmp,[],4);
    vol_wm_diff = vol_val_wm-vol_val_locmax;
    vol_wm_masked = vol_wm.*(vol_wm_diff>-parms.outlier_thresh*meanval_wm);
  else
    vol_tmp(vol_tmp==0) = nan;
    vol_val_locmax = min(vol_tmp,[],4);
    vol_wm_diff = vol_val_wm-vol_val_locmax;
    vol_wm_masked = vol_wm.*(vol_wm_diff<parms.outlier_thresh*meanval_wm);
  end
  
  % restore values for voxels in eroded wm mask
  if parms.outlier_edge_flag
    vol_wm_masked(vol_wm_eroded>0) = 1;
  end;
  
  % exclude voxels with poor fit last iteration
  vol_wm_masked = (~vol_exclude).*vol_wm_masked;
  
  % find bias field ratio
  tmp_mask = vol_wm_masked(i1min:i1max,i2min:i2max,i3min:i3max) .* (vol(i1min:i1max,i2min:i2max,i3min:i3max)>0.3*meanval_wm);
  vol_tmp_in = vol(i1min:i1max,i2min:i2max,i3min:i3max);
  vol_tmp_in(~tmp_mask) = NaN;
  tmp_dist = bwdist(tmp_mask);
  vol_tmp1 = smoothn(vol_tmp_in,parms.lambda1,'robust');
  vol_tmp2 = smoothn(vol_tmp_in,parms.lambda2,'robust');
  vol_w = logit((tmp_dist-2)*1.0,1); vol_w(tmp_dist==0) = 0;
  vol_tmp = max(0.1*meanval_wm,(1-vol_w).*vol_tmp1+vol_w.*vol_tmp2);
  vol_pred = nan(size(vol));
  vol_pred(i1min:i1max,i2min:i2max,i3min:i3max) = vol_tmp;
  vol_ratio = parms.vol_target./vol_pred;
  vol_corr = vol .* vol_ratio;
  vol_exclude = vol_wm_masked.*(abs(vol_corr-parms.vol_target)>0.2*meanval_wm);
end

% restrict range relative to median
if ~isempty(parms.ratio_range)
  vol_ratio(vol_ratio<med_ratio*parms.ratio_range(1)) = ...
    med_ratio*parms.ratio_range(1);
  vol_ratio(vol_ratio>med_ratio*parms.ratio_range(2)) = ...
    med_ratio*parms.ratio_range(2);
end;

% smooth bias field ratio
if parms.postsmooth_flag
  vol_ratio_sm = vol_ratio; vol_ratio_sm(~vol_wm_masked) = NaN;
  vol_ratio_sm = smoothn(vol_ratio_sm,0.1);
  vol_ratio = vol_ratio_sm;
end;

% use simple constant extrapolation outside bounding box
dims = size(vol);
vol_ratio_padded = vol_ratio;
vol_ratio_padded(1:i1min-1,:,:) = repmat(vol_ratio_padded(i1min,:,:),[i1min-1 1 1]);
vol_ratio_padded(i1max+1:dims(1),:,:) = repmat(vol_ratio_padded(i1max,:,:),[dims(1)-i1max 1 1]);
vol_ratio_padded(:,1:i2min-1,:) = repmat(vol_ratio_padded(:,i2min,:),[1 i2min-1 1]);
vol_ratio_padded(:,i2max+1:dims(2),:) = repmat(vol_ratio_padded(:,i2max,:),[1 dims(2)-i2max 1]);
vol_ratio_padded(:,:,1:i3min-1) = repmat(vol_ratio_padded(:,:,i3min),[1 1 i3min-1]);
vol_ratio_padded(:,:,i3max+1:dims(3)) = repmat(vol_ratio_padded(:,:,i3max),[1 1 dims(3)-i3max]);
vol_ratio = vol_ratio_padded;

% apply correction
vol_corr = vol.*vol_ratio;
vol_corr(isnan(vol_corr)) = 0;
vol_pred = parms.vol_target.*vol_ratio;
vol_pred(isnan(vol_pred)) = 0;

% clip values to input range
if ~isempty(parms.image_range)
  vol_corr(vol_corr<parms.image_range(1)) = parms.image_range(1);
  vol_corr(vol_corr>parms.image_range(2)) = parms.image_range(2);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol = erode_mask(vol,nvoxels,sigma,thresh)
  if ~exist('nvoxels','var') || isempty(nvoxels), nvoxels=1; end;
  if ~exist('sigma','var') || isempty(sigma), sigma=0.3; end;
  if ~exist('thresh','var') || isempty(thresh), thresh=0.999; end;
  for j=1:nvoxels
    vol = mmil_smooth3d(vol,sigma,sigma,sigma);
    vol = 1.0*(vol>=thresh);
  end;
return;

