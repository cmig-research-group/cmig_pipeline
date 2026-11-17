function fname_mgz = cmig_import_to_mgz(fname, outdir, fstem, forceflag)
%function fname_mgz = cmig_import_to_mgz(fname, outdir, [fstem], [forceflag])
%
% Purpose: convert a NIfTI (or other FreeSurfer-compatible) file to MGZ format
%   skips conversion if output already exists; returns empty if conversion fails
%
% Required Input:
%   fname: input file name
%   outdir: output directory
%
% Optional Input:
%   fstem: output file stem; if empty, derived from fname
%     {default = []}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_mgz: output MGZ file name; empty if input is empty or conversion failed
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Last Mod: 07/13/2026 by Don Hagler
%

fname_mgz = [];
if ~mmil_check_nargs(nargin,2), return; end
if ~exist('fstem','var') || isempty(fstem), fstem = []; end
if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end
if isempty(fname), return; end
if isempty(fstem)
  [~,fstem] = cmig_get_fileparts(fname);
end
fname_mgz = sprintf('%s/%s.mgz',outdir,fstem);
if ~exist(fname_mgz,'file') || forceflag
  fs_mri_convert(fname,fname_mgz);
end
if ~exist(fname_mgz,'file'), fname_mgz = []; end

return;
