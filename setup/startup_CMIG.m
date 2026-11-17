function startup_CMIG
%function startup_CMIG
%
% Purpose: add surface_morphometry_toolbox, gradient_distortion_fields and
%   cmig_pipeline's matlab/ to the MATLAB path, in that order (cmig_pipeline
%   highest), from SURFACE_MORPHOMETRY_TOOLBOX_DIR,
%   GRADIENT_DISTORTION_FIELDS_DIR and CMIG_PIPELINE_DIR (set by SetUpCMIG.csh
%   or SetUpMMPS.csh); .git folders are left out, an unset variable is skipped
%
% Usage: called by startup_MMPS.m; without MMPS, from startup.m:
%   run([getenv('CMIG_PIPELINE_DIR') '/setup/startup_CMIG.m']);
%
% Created:  10/27/2026 by Don Hagler
% Last Mod: 10/27/2026 by Don Hagler
%

dirlist = {getenv('SURFACE_MORPHOMETRY_TOOLBOX_DIR'),...
           getenv('GRADIENT_DISTORTION_FIELDS_DIR')};
cmig_dir = getenv('CMIG_PIPELINE_DIR');
if ~isempty(cmig_dir), dirlist{end+1} = [cmig_dir '/matlab']; end;
for i=1:length(dirlist)
  if isempty(dirlist{i}) || ~exist(dirlist{i},'dir'), continue; end;
  subdirs = regexp(genpath(dirlist{i}),pathsep,'split');
  subdirs = subdirs(~cellfun(@isempty,subdirs) & ...
                    cellfun(@isempty,regexp(subdirs,'/\.git(/|$)')));
  addpath(subdirs{:});
end;

return;
