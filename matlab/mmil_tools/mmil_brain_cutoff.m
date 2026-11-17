function [cutoff,dorsal,ventral]=mmil_brain_cutoff(RegInfo,fname_brain,varargin)
%function [cutoff,dorsal,ventral]=mmil_brain_cutoff(RegInfo,fname_brain,[options])
%
% Purpose: estimate the fraction of brain excluded
%          from imaging FOV
%
% Usage: [cutoff,dorsal,ventral] = mmil_brain_cutoff(RegInfo,fname_brain,'key1', value1,...);
%
% Required Input:
%   RegInfo: struct array with fields:
%     M_T1: vox2ras matrix for T1
%     volsz_T1: number of voxels for each dimension for T1
%     M_T2: vox2ras matrix for T2
%     volsz_T2: number of voxels for each dimension for T2
%     M_T1_to_T2: registration matrix between T1 and T2
%   fname_brain: full path of input brain mask volume file
%     This file should be an accurate mask of the brain
%     in the space of the T1
%         
% Optional Input:
%   'fspath': Freesurfer path
%     will use aseg.mgz to create brain mask
%     if fname_brain does not exist
%   'fname_T2': full path of the T2 (e.g. EPI) image (mgh or mgz)
%     if supplied, its nonzero voxels (first frame, holes filled) mark the
%     imaged region,
%     so zero padding of a resampled volume counts as outside the FOV;
%     if not, the whole T2 grid does
%     {default = []}
%   'fname_plot': full path of output image file (*.png)
%     If supplied, resulting image showing intersection of brain mask
%     with frames will be saved in fname_out
%   'forceflag': [0|1} overwrite existing files
%     {default = 0}
%
% Output:
%   cutoff: total cutoff score (sum of dorsal and ventral)
%   dorsal: cutoff score for dorsal edge (superior)
%   ventral: cutoff score for ventral edge (inferior)
%
% Created:  02/09/20 by Don Hagler
% Prev Mod: 02/13/20 by Don Hagler
% Last Mod: 10/01/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cutoff = nan; dorsal = nan; ventral = nan;
if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'fspath',[],[],...
  'fname_T2',[],[],...
  'fname_plot',[],[],...
  'forceflag',false,[false true],...
  ...
  'snum',120,[],...
  'fig_size',[5,5],[],...
  'img_dpi',300,[],...
  'visible_flag',false,[false true],...
});

if ~exist(fname_brain,'file')
  if isempty(parms.fspath)
    error('fname_brain %s does not exist and fspath was not supplied',fname_brain);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create aseg brain mask
if ~exist(fname_brain,'file')|| parms.forceflag
%  fprintf('%s: creating brain mask from FreeSurfer aseg...\n',...
%    mfilename);
%  tic
  p = [];
  p.fname_mask = fname_brain;
  p.brain_flag = 1; % set to 0 for binary mask
  p.forceflag = parms.forceflag;
  a = mmil_parms2args(p);
  mmil_aseg_brainmask(parms.fspath,a{:});
%  toc
end;
%fprintf('%s: loading brain mask from %s...\n',...
%  mfilename,fname_brain);
%tic
[vol_brain,M_brain] = fs_load_mgh(fname_brain);
%toc

% resample mask of all slices from T2 to T1
%fprintf('%s: resampling T2 acquisition volume to T1 space...\n',...
%  mfilename);
%tic
if ~isempty(parms.fname_T2)
  % the imaged region: nonzero voxels of the first frame, holes filled so
  %   isolated zero voxels inside the image count as imaged
  vol_T2 = 1.0*imfill(fs_load_mgh(parms.fname_T2,[],1)~=0,'holes');
else
  vol_T2 = ones(RegInfo.volsz_T2(1:3));
end;
p = [];
p.M_ref = RegInfo.M_T1;
p.nvox_ref = RegInfo.volsz_T1(1:3);
p.M_reg = RegInfo.M_T1_to_T2;
a = mmil_parms2args(p);
[vol_T2,M_T2] = mmil_resample_vol(vol_T2,RegInfo.M_T2,a{:});
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate total cutoff
n_brain = nnz(vol_brain>0);
n_brain_T2 = nnz(vol_brain>0 & vol_T2>0);
cutoff = (n_brain - n_brain_T2) / n_brain;

% split brain mask into dorsal and ventral halves
avg_slice = squeeze(mean(vol_brain,1));
avg_row = squeeze(mean(avg_slice,2));
mask_vec = avg_row>0;
center = round(mmil_wtd_mean([1:length(avg_row)]',avg_row));
vol_dorsal = vol_brain;
vol_dorsal(:,center:end,:) = 0;
vol_ventral = vol_brain;
vol_ventral(:,1:center,:) = 0;

% calculate dorsal cutoff
n_dorsal = nnz(vol_dorsal>0);
n_dorsal_T2 = nnz(vol_dorsal>0 & vol_T2>0);
dorsal = (n_dorsal - n_dorsal_T2) / n_dorsal;

% calculate ventral cutoff
n_ventral = nnz(vol_ventral>0);
n_ventral_T2 = nnz(vol_ventral>0 & vol_T2>0);
ventral = (n_ventral - n_ventral_T2) / n_ventral;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(parms.fname_plot) &&...
  (~exist(parms.fname_plot) || parms.forceflag)
  s_T2 = squeeze(vol_T2(parms.snum,:,:));
  s_brain = squeeze(vol_brain(parms.snum,:,:));
  s_rgb = zeros([size(s_brain),3]);
  s_rgb(:,:,1) = s_brain/128 + (~s_T2)/4;
  s_rgb(:,:,2) = s_brain/128 + s_T2/4;
  s_rgb(:,:,3) = s_brain/128 + s_T2/7;
  figure; clf;
  imagesc(s_rgb);
  title(sprintf('cutoff = %0.2f%%',cutoff*100));
  ylabel(sprintf('dorsal cutoff = %0.1f%%',dorsal*100));
  xlabel(sprintf('ventral cutoff = %0.1f%%',ventral*100));
  mmil_save_fig(gcf,parms.fname_plot,...
    parms.fig_size,parms.img_dpi,0,parms.visible_flag);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

