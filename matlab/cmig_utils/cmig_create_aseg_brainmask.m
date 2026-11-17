function fname_brainmask = cmig_create_aseg_brainmask(fspath, outdir, verbose, forceflag)
%function fname_brainmask = cmig_create_aseg_brainmask(fspath, outdir, [verbose], [forceflag])
%
% Purpose: create a brain mask from FreeSurfer aseg.mgz
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
%   fname_brainmask: output brain mask file name ({outdir}/aseg_brainmask.mgz)
%
% Created:  05/21/2026 by Don Hagler
% Last Mod: 08/20/2026 by Don Hagler
%

fname_brainmask = [];
if ~mmil_check_nargs(nargin,2), return; end
if ~exist('verbose','var')  || isempty(verbose),  verbose  = false; end
if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end

fname_brainmask = sprintf('%s/aseg_brainmask.mgz',outdir);
if ~exist(fname_brainmask,'file') || forceflag
  if verbose
    fprintf('%s: creating brain mask from FreeSurfer aseg...\n',mfilename);
    tic;
  end
  p = [];
  p.fname_mask = fname_brainmask;
  p.brain_flag = 0;
  p.forceflag = forceflag;
  a = mmil_parms2args(p);
  mmil_aseg_brainmask(fspath,a{:});
  if verbose && ~isempty(fname_brainmask), toc; end
end

return;
