function fname_cc = cmig_create_aseg_ccmask(fspath, outdir, verbose, forceflag)
%function fname_cc = cmig_create_aseg_ccmask(fspath, outdir, [verbose], [forceflag])
%
% Purpose: create a dilated corpus callosum mask from FreeSurfer aseg.mgz
%
% Required Input:
%   fspath: FreeSurfer subject directory
%   outdir: output directory
%
% Optional Input:
%   verbose: [0|1] display status messages
%     {default = 0}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_cc: output CC mask file name ({outdir}/aseg_corpus_callosum.mgh)
%
% Created:  05/21/2026 by Don Hagler
% Last Mod: 05/21/2026 by Don Hagler
%

fname_cc = [];
if ~mmil_check_nargs(nargin,2), return; end
if ~exist('verbose','var')  || isempty(verbose),  verbose  = false; end
if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end

cc_roi_vals = [251 252 253 254 255];
fname_cc = sprintf('%s/aseg_corpus_callosum.mgh',outdir);
if ~exist(fname_cc,'file') || forceflag
  if verbose
    fprintf('%s: creating corpus callosum mask from FreeSurfer aseg...\n',mfilename);
    tic;
  end
  fname_aseg = sprintf('%s/mri/aseg.mgz',fspath);
  [vol_aseg,M] = fs_load_mgh(fname_aseg);
  i_mask = find(ismember(vol_aseg,cc_roi_vals));
  vol_mask = zeros(size(vol_aseg));
  vol_mask(i_mask) = 1;
  vol_mask_dil = mmil_dilate_mask(ctx_mgh2ctx(vol_mask,M),...
    'smooth1',45,'thresh1',0.02,'smooth2',0,'smooth3',0);
  ctx_save_mgh(vol_mask_dil,fname_cc);
  if verbose && ~isempty(fname_cc), toc; end
end

return;
