function [regT1_corr,regT1_norm,regT1_warp,results,errcode] = mmil_check_regT1(fname_in,fname_T1,varargin)
%function [regT1_corr,regT1_norm,regT1_warp,results,errcode] = mmil_check_regT1(fname_in,fname_T1,[options])
%
% Purpose: compute the voxel-wise correlation between two images
%          optionally within brain mask
%
% Usage: regT1_corr = mmil_check_regT1(fname_in,fname_T1, 'key1', value1,...);
%
% Required Input:
%   fname_in   : full path of input image volume file (mgz)
%   fname_T1   : full path of T1-weighted volume file (mgz)
%
% Optional Input:
%   'fname_mask': full path of mask volume file (mgz)
%     if not supplied, will include all voxels
%     {default = []}
%   'fname_plot': full path of output image file (jpg or tif)
%     {default = []}
%   'visible_flag': [0|1] show figure on screen
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   regT1_corr: correlation between input images within mask
%   regT1_norm: combined measure of rigid (trans+rot) transformation
%   regT1_warp: measure of quasi-affine warp (xy strech) 
%   errcode: returned as 1 if there was an error computing metric, 0 otherwise
%
% Created:  03/02/20 by Don Hagler
% Prev Mod: 05/01/20 by Payam Saisan
% Last Mod: 05/01/20 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0; regT1_corr = nan; regT1_norm = nan; regT1_warp = nan;
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
    'fname_mask',[],[],...
    'fname_plot',[],[],...
    'forceflag',false,[false true],...
    ...
    'title','Check Registration to T1',[],...
    'fig_size',[7,5],[],...
    'img_dpi',300,[],...
    'visible_flag',false,[false true],...
    ...
    'vox2ras_maxdiff',1e-3,[],...
    ...
    'tags',{'fname_plot','forceflag',...
            'title','fig_size','img_dpi','visible_flag'},[],...
    });

% load input image (first frame only)
[vol_in,M_in,mrp,volsz_in] = fs_load_mgh(fname_in,[],1);

% load T1-weighted image
[vol_T1,M_T1,mrp,volsz_T1] = fs_load_mgh(fname_T1,[],1);

% compare image dimension and M matrics
errcode = check_vols(fname_in,volsz_in,M_in,...
                     fname_T1,volsz_T1,M_T1,parms.vox2ras_maxdiff);
if errcode, return; end;

% load mask
if ~isempty(parms.fname_mask)
  [vol_mask,M_mask,mrp,volsz_mask] = fs_load_mgh(parms.fname_mask,[],1);
  vol_mask = 1.0*(vol_mask>1e-5);
else
  vol_mask = ones(volsz_T1(1:3));
  M_mask = M_T1;
  volsz_mask = volsz_T1;
end;

% compare image dimension and M matrics
errcode = check_vols(parms.fname_mask,volsz_mask,M_mask,...
                     fname_T1,volsz_T1,M_T1,parms.vox2ras_maxdiff);
if errcode, return; end;

% compute corr between images within mask
regT1_corr = abs(corr(vol_in(vol_mask>0),vol_T1(vol_mask>0)));

% estimate quasi-affine transformation between vol_in and vol_T1
%   using bounding box around masked region as a window of interest
a = mmil_parms2args(parms,parms.tags);
[results,errcode] = mmil_check_affine_trans(vol_in,vol_T1,vol_mask,a{:});

regT1_norm = results.TotRegT1Norm;
regT1_warp = results.WarpNorm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compare image dimension and M matrics
function errcode = check_vols(fname_A,volsz_A,M_A,fname_B,volsz_B,M_B,vox2ras_maxdiff)
  errcode = 0;
  % compare image dimensions
  if any(volsz_A(1:3)~=volsz_B(1:3))
    fprintf('%s: WARNING: volume dimensions do not match for %s and %s\n',...
        mfilename,fname_A,fname_B);
    errcode = 1;
    return;
  end;
  % compare M matrices
  if any(abs(M_A(:)-M_B(:))>vox2ras_maxdiff)
    fprintf('%s: WARNING: vox2ras matrices do not match for %s and %s\n',...
        mfilename,fname_A,fname_B);
    errcode = 1;
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
