function [cost0,errcode,errmsg] = cmig_check_regT1_segreg(fname_mov,fspath,varargin)
%function [cost0,errcode,errmsg] = cmig_check_regT1_segreg(fname_mov,fspath,[options])
%
% Purpose: boundary-based registration cost of a volume already resampled
%   to T1w space, from mri_segreg with no optimization (--nmax 0) along the
%   FreeSurfer white surfaces
%
% Required Input:
%   fname_mov: volume resampled to T1w space (mgh/mgz format)
%   fspath: FreeSurfer subject directory of the T1w
%
% Optional Input ('key',value):
%   'contrast': ['bold'|'t2'] mri_segreg contrast option
%     'bold' for EPI (gray matter brighter than white matter)
%     't2' for T2-weighted volumes
%     {default = 'bold'}
%   'tmpdir': working directory for the registration and cost files
%     {default = pwd}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   cost0: mri_segreg cost of the existing registration (NaN on error)
%   errcode: 0 if successful, 1 if error
%   errmsg: error message string (empty if no error)
%
% Notes:
%   the registration is tkregister2 --regheader between fname_mov and the
%     subject's orig.mgz; the identity when both share scanner coordinates
%
% Created:  10/05/2026 by Don Hagler
% Last Mod: 10/05/2026 by Don Hagler
%

cost0 = NaN; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,2), errcode = 1; return; end
parms = mmil_args2parms(varargin,{...
  'contrast','bold',{'bold','t2'},...
  'tmpdir',pwd,[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
});

[subjdir,subj] = fileparts(fspath);
fnames_need = {fname_mov,sprintf('%s/mri/orig.mgz',fspath)};
for h={'lh','rh'}
  fnames_need{end+1} = sprintf('%s/surf/%s.white',fspath,h{1});
  fnames_need{end+1} = sprintf('%s/label/%s.cortex.label',fspath,h{1});
end
for f=1:length(fnames_need)
  if ~exist(fnames_need{f},'file')
    errmsg = sprintf('file %s not found',fnames_need{f});
    errcode = 1;
    return;
  end
end

[~,fstem] = cmig_get_fileparts(fname_mov);
outdir = sprintf('%s/%s_segreg',parms.tmpdir,fstem);
mmil_mkdir(outdir);
fname_reg = sprintf('%s/reg.dat',outdir);
fname_cost = sprintf('%s/cost0.init',outdir);

if ~exist(fname_cost,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: mri_segreg cost of %s...\n',mfilename,fname_mov);
  end
  cmd = sprintf(['setenv SUBJECTS_DIR %s; ',...
    'tkregister2 --mov %s --targ %s/mri/orig.mgz --regheader --reg %s --s %s --noedit'],...
    subjdir,fname_mov,fspath,fname_reg,subj);
  [s,r] = mmil_unix(cmd);
  if s || ~exist(fname_reg,'file')
    errmsg = sprintf('tkregister2 failed:\n%s',r);
    errcode = 1;
    return;
  end
  cmd = sprintf(['setenv SUBJECTS_DIR %s; ',...
    'mri_segreg --mov %s --init-reg %s --s %s --out-reg %s/cost0.dat --%s ',...
    '--interp trilinear --wm-proj-abs 2 --c0 0 --dof 6 --surf white ',...
    '--nmax 0 --initcost %s'],...
    subjdir,fname_mov,fname_reg,subj,outdir,parms.contrast,fname_cost);
  [s,r] = mmil_unix(cmd);
  if s || ~exist(fname_cost,'file')
    errmsg = sprintf('mri_segreg failed:\n%s',r);
    errcode = 1;
    return;
  end
end

% cost0.init: cost, white matter mean, gray matter mean, percent contrast
vals = load(fname_cost);
if isempty(vals)
  errmsg = sprintf('%s is empty',fname_cost);
  errcode = 1;
  return;
end
cost0 = vals(1);

return;

