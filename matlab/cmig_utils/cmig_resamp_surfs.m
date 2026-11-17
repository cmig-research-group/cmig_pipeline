function errcode = cmig_resamp_surfs(fspath,outdir,forceflag)
%function errcode = cmig_resamp_surfs(fspath,outdir,[forceflag])
%
% Required parameters:
%   fspath: full path of FreeSurfer output directory
%   outdir: output directory
%
% Optional parameters:
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  10/31/2025 by Don Hagler
% Prev Mod: 11/04/2025 by Don Hagler
% Last Mod: 07/13/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,2), return; end
if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end

hemilist = {'lh','rh'};
surflist = {'white','pial'};

if ~exist(fspath,'dir')
  fprintf('%s: WARNING: fspath %s not found\n',mfilename,fspath);
  errcode = 1;
  return;
end

% get location of fsverage
fsavg = sprintf('%s/subjects/fsaverage',getenv('FREESURFER_HOME'));
if ~exist(fsavg,'dir')
  fprintf('%s: WARNING: fsaverage %s not found\n',mfilename,fsavg);
  errcode = 1;
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(outdir);

for h=1:length(hemilist)
  hemi = hemilist{h};
  fname_subj_reg = sprintf('%s/surf/%s.sphere.reg',fspath,hemi);
  fname_atlas_reg = sprintf('%s/surf/%s.sphere.reg',fsavg,hemi);
  for s=1:length(surflist)
    surfname = surflist{s};
    fname_in = sprintf('%s/surf/%s.%s',fspath,hemi,surfname);
    if ~exist(fname_in,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
      errcode = 1;
      return;
    end
    fname_out = sprintf('%s/%s.%s',outdir,hemi,surfname);
    if ~exist(fname_out,'file') || forceflag
      fprintf('%s: resampling %s to fsaverage...\n',mfilename,fname_in);
      cmd = sprintf('mris_resample --atlas_reg %s --subject_reg %s --subject_surf %s --out %s',...
        fname_atlas_reg,fname_subj_reg,fname_in,fname_out);
      [s,r] = mmil_unix(cmd);
      if s
        fprintf('%s: WARNING: cmd %s failed:\n%s\n',mfilename,cmd,r);
        errcode = 1;
        return;
      end
    end
  end
end


