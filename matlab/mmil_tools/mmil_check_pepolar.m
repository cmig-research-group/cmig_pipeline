function [mi_diff,mi_for,mi_rev]=mmil_check_pepolar(fname_main,fname_for,fname_rev,varargin)
%function [mi_diff,mi_for,mi_rev]=mmil_check_pepolar(fname_main,fname_for,fname_rev,varargin)
%
% Purpose: determine the polarity conformance of main scan vs forward or
%          reverse scans
%         
% Usage:  [mi_diff,mi_for,mi_rev] = mmil_check_pepolar(fname_main,fname_for,fname_rev,varargin)
%
% Required Inputs:
%   fname_main : 3D volume for main scan
%   fname_for  : 3D volume for forward sca
%   fname_rev  : 3D volume for reverse scan
%
% Optional Input:
%   'mask_flag': [0|1] restrict mutual information calculations to voxels
%     within the region defined by the mask
%     {default = 1}
%   'fname_mask': file name of mask file matching fname_main
%     If not suppplied, will generate a mask from fname_main
%       based on the bottom quartile of intensities
%     {default = []}
%   'reg_flag': [0|1] register images before comparing images
%     {default = 1}
%   'fname_plot': full path of output image file (use jpg or tif extension)
%     If supplied, save montage of slices from three input volumes
%     {default = []}
%   'forceflag': [0|1} overwrite existing files
%     {default = 0}
%
% Output:
%   mi_diff: normalized difference in mutual information (MI)
%            (forward to main) vs (reverse to main)
%     if >> 0  => main scan is in better alignmnet with forward scan  
%     if << 0  => main scan is in better alignment with reverse scan 
%     if ~= 0  => ambgious polarity ("forward" and "reverse" may be same)
%   mi_for: MI between forward and main scan
%   mi_rev: MI between reverse and main scan
%
% Created:  02/12/2020 by Payam Saisan
% Prev Mod: 02/23/2020 by Don Hagler
% Last Mod: 03/04/2020 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mi_diff = nan; mi_for = nan; mi_rev = nan;
if ~mmil_check_nargs(nargin,3), return; end;

parms = mmil_args2parms(varargin,{...
  'mask_flag',true,[false true],...
  'fname_mask',[],[],...
  'reg_flag',true,[false true],...
  'fname_plot',[],[],...
  'forceflag',false,[false true],...
  ...
  'plane','SAG',{'SAG','COR','HOR'},...
  'slice_fract',0.5,[0,1],...
  'fig_size',[5,5],[],...
  'img_dpi',300,[],...
  'visible_flag',false,[false true],...
  ...
  'mstep',4,[],...
  'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  'min_trans',0.05,[0,1],... % mm
  'min_rot',0.05,[0,1],... % degrees
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the three volumes (first frame only)
[vol_for,M_for,~,volsz_for] = fs_load_mgh(fname_for,[],1);
[vol_rev,M_rev,~,volsz_rev] = fs_load_mgh(fname_rev,[],1);
[vol_main,M_main,~,volsz_main] = fs_load_mgh(fname_main,[],1);

% load or create mask volume
if ~parms.mask_flag
  idx=1:numel(vol_main);
else
  if ~isempty(parms.fname_mask)
    % load vol mask
    [vol_mask,M_mask,~,volsz_mask] = fs_load_mgh(parms.fname_mask);
    % check whether volsz matches mask and main
    if any(volsz_mask(1:3)~=volsz_main(1:3))
      % resample mask to main voxels
      [vol_mask,M_mask] = mmil_resample_vol(vol_mask,M_mask,...
        'M_ref',M_main,'nvox_ref',volsz_main);
    end;
    idx = find(vol_mask(:)>0.5);
  else
    % create vol mask
    [mm,nn] = hist(vol_main(:),4);
    idx = find(vol_main>nn(1));
  end
end;
vol_mask = zeros(size(vol_main));
vol_mask(idx) = 1;

if parms.reg_flag

  % convert to ctx structs
  volmain = ctx_mgh2ctx(vol_main,M_main);
  volfor = ctx_mgh2ctx(vol_for,M_for);
  volrev = ctx_mgh2ctx(vol_rev,M_rev);
  volmask = ctx_mgh2ctx(vol_mask,M_main);
  
  % register for to main
  M_main_to_for = rbreg_vol2vol_icorr_mask(volmain,volfor,...
        volmask,0,parms.mstep,0,parms.scales,...
        parms.min_trans,parms.min_rot*pi/180);
  % resample for to main
  [vol_for,M_for] = mmil_resample_vol(vol_for,M_for,...
    'M_ref',M_main,'nvox_ref',volsz_main,'M_reg',M_main_to_for);

  % register rev to main
  M_main_to_rev = rbreg_vol2vol_icorr_mask(volmain,volrev,...
        volmask,0,parms.mstep,0,parms.scales,...
        parms.min_trans,parms.min_rot*pi/180);
  % resample rev to main
  [vol_rev,M_rev] = mmil_resample_vol(vol_rev,M_rev,...
    'M_ref',M_main,'nvox_ref',volsz_main,'M_reg',M_main_to_rev);
    
else

  % check whether volsz matches for and main
  if any(volsz_for(1:3)~=volsz_main(1:3))
    % resample "for" to main resolution
    [vol_for,M_for] = mmil_resample_vol(vol_for,M_for,...
      'M_ref',M_main,'nvox_ref',volsz_main);
  end;

  % check whether volsz matches rev and main
  if any(volsz_rev(1:3)~=volsz_main(1:3))
    % resample "rev" to main voxels
    [vol_rev,M_rev] = mmil_resample_vol(vol_rev,M_rev,...
      'M_ref',M_main,'nvox_ref',volsz_main);
  end;
  
end;

% calculate mutual information metrics
mi_for = mmil_information(vol_for(idx)',vol_main(idx)');
mi_rev = mmil_information(vol_rev(idx)',vol_main(idx)');  
mi_diff = (mi_for-mi_rev)*2/(mi_rev+mi_for); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot images for later review
if ~isempty(parms.fname_plot) &&...
  (~exist(parms.fname_plot) || parms.forceflag)
  
  % reorient to standard orientation
  vol_main_PRS = fs_reorient(vol_main,M_main,'PRS');
  vol_mask_PRS = fs_reorient(vol_mask,M_main,'PRS');
  vol_for_PRS = fs_reorient(vol_for,M_for,'PRS');
  vol_rev_PRS = fs_reorient(vol_rev,M_rev,'PRS');

  % extract slices  
  s_main = extract_slice(vol_main_PRS,parms.plane,parms.slice_fract);
  s_mask = extract_slice(vol_mask_PRS,parms.plane,parms.slice_fract);
  s_for = extract_slice(vol_for_PRS,parms.plane,parms.slice_fract);
  s_rev = extract_slice(vol_rev_PRS,parms.plane,parms.slice_fract);

  figure; clf;
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  
  subplot(2,2,1); % MAIN
  imagesc(s_main);
  colormap(bone);
  title('Main Scan');
  set(gca,'XTickLabel','','YTickLabel','');
  xlabel(sprintf('norm MI diff = %0.2f',mi_diff));
  
  subplot(2,2,2); % FOR
  imagesc(s_for);
  colormap(bone);
  title('Forward Scan');
  set(gca,'XTickLabel','','YTickLabel','');
  xlabel(sprintf('MI(Main,Fwd) = %0.2f',mi_for));

  subplot(2,2,3); % MASK
  imagesc(s_mask);
  colormap(bone);
  title('Mask');
  set(gca,'XTickLabel','','YTickLabel','');

  subplot(2,2,4); % REV
  imagesc(s_rev);
  colormap(bone);
  title('Reverse Scan');
  set(gca,'XTickLabel','','YTickLabel','');
  xlabel(sprintf('MI(Main,Rev) = %0.2f',mi_rev));

  mmil_save_fig(gcf,parms.fname_plot,...
    parms.fig_size,parms.img_dpi,0,parms.visible_flag);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img = extract_slice(vol,plane,slice_fract)
  dim = size(vol);
  s = choose_slice(dim,plane,slice_fract);
  switch plane
    case 'COR'
      img = rot90(squeeze(vol(s,:,:)));
    case 'SAG'
      img = rot90(squeeze(vol(:,s,:)));
    case 'HOR'
      img = squeeze(vol(:,:,s));
  end
return;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slice_num = choose_slice(dim,plane,slice_fract)
  switch plane
    case 'COR'
      d = 1;
    case 'SAG'
      d = 2;
    case 'HOR'
      d = 3;
  end
  slice_num = round(dim(d)*slice_fract);
  slice_num = min(max(slice_num,1),dim(d));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
