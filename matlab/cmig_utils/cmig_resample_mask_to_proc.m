function fname_mask_proc = cmig_resample_mask_to_proc(fname_mask, fname_ref, outdir, verbose, forceflag)
%function fname_mask_proc = cmig_resample_mask_to_proc(fname_mask, fname_ref, outdir, [verbose], [forceflag])
%
% Purpose: resample a mask volume to match the space and resolution of a reference volume
%
% Required Input:
%   fname_mask: mask file name to resample
%   fname_ref: reference volume defining target space and resolution
%   outdir: output directory
%
% Optional Input:
%   verbose: [0|1] display status messages
%     {default = 0}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_mask_proc: output resampled mask ({outdir}/aseg_corpus_callosum_proc.mgh)
%
% Created:  05/21/2026 by Don Hagler
% Last Mod: 05/21/2026 by Don Hagler
%

fname_mask_proc = [];
if ~mmil_check_nargs(nargin,3), return; end
if ~exist('verbose','var')  || isempty(verbose),  verbose  = false; end
if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end

fname_mask_proc = sprintf('%s/aseg_corpus_callosum_proc.mgh',outdir);
if ~exist(fname_mask_proc,'file') || forceflag
  if verbose
    fprintf('%s: resampling CC mask to proc space...\n',mfilename);
    tic;
  end
  [vol_mask,M_mask] = fs_load_mgh(fname_mask);
  [M_ref,nvox_ref] = mmil_load_mgh_info(fname_ref,forceflag,outdir);
  [vol_res,M_res] = mmil_resample_vol(vol_mask,M_mask,...
    'M_ref',M_ref,'nvox_ref',nvox_ref(1:3),'interpm',1,'bclamp',1);
  vol_res(vol_res<0.5) = 0;
  fs_save_mgh(vol_res,fname_mask_proc,M_res);
  if verbose && ~isempty(fname_mask_proc), toc; end
end

return;
