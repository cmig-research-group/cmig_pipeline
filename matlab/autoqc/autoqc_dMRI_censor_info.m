function [censor_info,censor_data,err_vol_norm] = autoqc_dMRI_censor_info(fit,err,varargin)
%function [censor_info,censor_data,err_vol_norm] = autoqc_dMRI_censor_info(fit,err,[options])
%
% Created:  10/04/22 by Don Hagler
% Prev Mod: 10/11/22 by Don Hagler
% Prev Mod: 11/01/22 by Don Hagler
% Last Mod: 01/30/23 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

censor_info = []; censor_data = []; err_vol_norm = [];

if ~mmil_check_nargs(nargin,2), return; end

parms = mmil_args2parms(varargin,{...
  'norm_bvals_flag',true,[false true],...
  'norm_slices_flag',true,[false true],...
  'winsorize_flag',true,[false true],...
  'clip_neg_flag',true,[false true],...
  'erode_nvoxels',0,[0,10],...
  'smooth_nvoxels',0,[0,100],...
  'censor_thresh',[],[],...
  'censor_prctile',99,[],...
  'min_fract',0.05,[],...
  'sigma',1,[1e-2,10],...
  'thresh',0.99,[0.1,0.999],...
});  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% erode brain mask
if parms.erode_nvoxels>0
  volmask = 1.0*(fit.volmask>0);
  for j=1:parms.erode_nvoxels
    volmask = mmil_smooth3d(volmask,parms.sigma,parms.sigma,parms.sigma);
    volmask = 1.0*(volmask>=parms.thresh);
  end
else
  volmask = fit.volmask;
end

% set volume dimensions
nx = fit.volsz(1); ny = fit.volsz(2); nz = fit.volsz(3); nf = fit.volsz(4);

% get unique b-values
bvals = fit.bvals;
round_bvals = 50*round(bvals/50);
unique_bvals = unique(round_bvals);
nbvals = length(unique_bvals);

% initial error volume
err_vol = err.vol;

% apply mask to error volume
err_vol = bsxfun(@times,err_vol,volmask);

% Winsorize error (clip values at +- 3 stdev) within each 3D volume
if parms.winsorize_flag
  for f=1:nf
    for z = 1:nz
      mask_vec = mmil_rowvec(squeeze(volmask(:,:,z)));
      if nnz(mask_vec)/length(mask_vec) < parms.min_fract
        continue;
      end
      img_vec = mmil_rowvec(squeeze(err_vol(:,:,z,f)));
      idx_mask = find(mask_vec>0);
      vals = img_vec(idx_mask);
      mean_val = mean(vals);
      std_val = std(vals);
      lo_val = mean_val - 3*std_val;
      hi_val = mean_val + 3*std_val;
      img_vec(img_vec < lo_val) = lo_val;
      img_vec(img_vec > hi_val) = hi_val;
      err_vol(:,:,z,f) = reshape(img_vec,[nx,ny]);
    end
  end
end

% apply 2D smoothing within slice
if parms.smooth_nvoxels>0
  for f=1:nf
    for z = 1:nz
      img_vec = mmil_rowvec(squeeze(err_vol(:,:,z,f)));
      img_vec = mmil_smooth2d(img_vec,parms.smooth_nvoxels,parms.smooth_nvoxels);
      err_vol(:,:,z,f) = reshape(img_vec,[nx,ny]);
    end
  end
end

% set negative err values to 0
if parms.clip_neg_flag
  err_vol(err_vol<0) = 0;
end

% calculate rms of residual error for each frame-slice
rms_err_mat = zeros(nz,nf);
for z = 1:nz
  mask_vec = mmil_rowvec(squeeze(volmask(:,:,z)));
  if nnz(mask_vec)/length(mask_vec) < parms.min_fract
    continue;
  end
  err_mat = reshape(err_vol(:,:,z,:),[nx*ny,nf])';
  err_mat = err_mat(:,mask_vec>0);
  rms_err_mat(z,:) = rms(err_mat,2);
end

% normalize by b-value
rms_err_mat_norm_bval = rms_err_mat;
err_vol_norm_bval = err_vol;
if parms.norm_bvals_flag
  % calculate median error within b-values
  med_err_bval = zeros(1,nbvals);
  for b=1:nbvals
    i_b = find(round_bvals == unique_bvals(b));
    tmp_err = median(rms_err_mat(:,i_b),2);
    med_err_bval(b) = median(tmp_err(tmp_err>0));
    % normalize by median error for each b-value
    rms_err_mat_norm_bval(:,i_b) = rms_err_mat(:,i_b) / med_err_bval(b);
    err_vol_norm_bval(:,:,:,i_b) = err_vol(:,:,:,i_b) / med_err_bval(b);
  end
elseif ~parms.norm_slices_flag
  % calculate median error across non-zero b-values
  i_b = find(round_bvals ~= 0);
  tmp_err = median(rms_err_mat(:,i_b),2);
  med_err_bval = median(tmp_err(tmp_err>0));
  rms_err_mat_norm_bval(:,i_b) = rms_err_mat(:,i_b) / med_err_bval;
  err_vol_norm_bval(:,:,:,i_b) = err_vol(:,:,:,i_b) / med_err_bval;
  rms_err_mat_norm_bval(:,fit.i_b0) = 0;
  err_vol_norm_bval(:,:,:,fit.i_b0) = 0;
end

% normalize by slice
if parms.norm_slices_flag
  % normalize by median error within slice across diffusion weighted frames
  med_err_slice = median(rms_err_mat(:,fit.i_bx),2);
  rms_err_mat_norm_slice = rms_err_mat ./ max(eps,med_err_slice*ones(1,nf));
  rms_err_mat_norm_slice(:,fit.i_b0) = 0;

  % normalize by median error within slice across diffusion weighted frames
  med_err_bval_slice = median(rms_err_mat_norm_bval,2);
  rms_err_mat_norm_bval_slice = rms_err_mat_norm_bval ./ max(eps,med_err_bval_slice*ones(1,nf));
  if ~parms.norm_bvals_flag
    rms_err_mat_norm_bval_slice(:,fit.i_b0) = 0;
  end

  % normalize err vol by median error by slice
  err_vol_norm_bval_slice = err_vol_norm_bval;
  tmp_err = med_err_bval_slice;
  tmp_err(tmp_err==0) = median(med_err_bval_slice);
  tmp_vol = repmat(reshape(tmp_err,[1,1,nz]),[nx,ny,1,nf]);
  err_vol_norm_bval_slice = err_vol_norm_bval ./ tmp_vol;
  if ~parms.norm_bvals_flag
    err_vol_norm_bval_slice(:,:,:,fit.i_b0) = 0;
  end

  % set final censor error
  censor_err = rms_err_mat_norm_bval_slice;
  err_vol_norm = err_vol_norm_bval_slice;
else
  % set final censor error
  censor_err = rms_err_mat_norm_bval;
  err_vol_norm = err_vol_norm_bval;
end

% set censor_thresh
if isempty(parms.censor_thresh)
  parms.censor_thresh = prctile(censor_err(censor_err>0),parms.censor_prctile);
  fprintf('%s: %0.0f%%tile censor threshold = %0.2f\n',...
    mfilename,parms.censor_prctile,parms.censor_thresh);
else
  fprintf('%s: censor threshold = %0.2f\n',...
    mfilename,parms.censor_thresh);
end

% set final censor error and censor frames
censor_mat = (censor_err>parms.censor_thresh);

% save b-values and unique b-values
censor_info.bvals = mmil_rowvec(bvals);
censor_info.unique_bvals = unique_bvals;
censor_info.nframes = nf;
censor_info.nslices = nz;

% calculate censoring metrics
censor_info.censor_thresh = parms.censor_thresh;
censor_info.nbad_frame_slices = sum(sum(censor_mat));
censor_info.nbad_frames = sum(sum(censor_mat,1) > 0);
censor_info.nbad_slices = sum(sum(censor_mat,2) > 0);
censor_info.nbad_frames_per_slice = sum(censor_mat,2);
censor_info.max_nbad_frames_per_slice = ...
  max(censor_info.nbad_frames_per_slice);
censor_info.nbad_slices_per_frame = sum(censor_mat,1);
censor_info.max_nbad_slices_per_frame = ...
  max(censor_info.nbad_slices_per_frame);

% calculate censoring metrics per shell
censor_info.nframes_per_bval = zeros(1,nbvals);
censor_info.nbad_frame_slices_per_bval = zeros(1,nbvals);
censor_info.nbad_frames_per_bval = zeros(1,nbvals);
censor_info.nbad_slices_per_bval = zeros(1,nbvals);
for b=1:nbvals
  i_b = find(round_bvals == unique_bvals(b));
  censor_info.nframes_per_bval(b) = length(i_b);
  censor_info.nbad_frame_slices_per_bval(b) = sum(sum(censor_mat(:,i_b)));
  censor_info.nbad_frames_per_bval(b) = sum(sum(censor_mat(:,i_b),1) > 0);
  censor_info.nbad_slices_per_bval(b) = sum(sum(censor_mat(:,i_b),2) > 0);
end

% save censoring data
censor_data.censor_err = censor_err;
censor_data.censor_mat = censor_mat;

censor_data.err_raw = rms_err_mat;
if parms.norm_slices_flag
  censor_data.err_norm_slice = rms_err_mat_norm_slice;
end
if parms.norm_bvals_flag
  censor_data.err_norm_bval = rms_err_mat_norm_bval;
end
if parms.norm_bvals_flag && parms.norm_slices_flag
  censor_data.err_norm_bval_slice = rms_err_mat_norm_bval_slice;
end


