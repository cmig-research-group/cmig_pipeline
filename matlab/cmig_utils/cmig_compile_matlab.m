function cmig_compile_matlab(varargin)
%function cmig_compile_matlab([options])
%
% Purpose: compile a cmig_pipeline MATLAB function as a standalone
%   application, with only the toolbox, the gradient distortion fields and
%   the cmig_pipeline tree on the MATLAB path (no MMPS)
%
% Optional Parameters ('key',value pairs):
%   'funcname': name of the function to compile
%     must be in pipeline_dir/matlab or matlab_funcdir, or given as full
%     path file name
%     {default = 'cmig_run_bids_module'}
%   'pipeline_dir': cmig_pipeline tree (matlab/ under it)
%     {default = [repos '/cmig_pipeline_dev']}
%   'toolbox_dir': surface_morphometry_toolbox tree
%     {default = [repos '/surface_morphometry_toolbox']}
%   'gradient_dir': gradient_distortion_fields tree
%     {default = [repos '/gradient_distortion_fields']}
%   repos: the environment variable CMIG_REPOS_DIR, or else the folder
%     holding the cmig_pipeline tree this file is in
%   'toolbox_data_dir': directory with the toolbox data files, in the
%     toolbox's folder layout (MorphometryToolbox/, SurfaceGeometry/); the
%     data are not in the toolbox repository
%     {default = []}
%   'matlab_funcdir': directory of MATLAB functions taking precedence over
%     the others, e.g. working copies
%     {default = []}
%   'outdir': output directory
%     {default = [pwd '/output']}
%
% The path, lowest first: toolbox, gradient_dir, pipeline_dir/matlab,
%   matlab_funcdir. Added as files the dependency analysis does not find:
%   calmi_hatl.m, named in a string by mmil_dct_morph; the toolbox's
%   icosahedron.off; and the toolbox data files the code loads by name
%   (from toolbox_data_dir, required).
%
% NOTE: requires Matlab R2023a or higher
%
% Created:  10/26/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Last Mod: 11/16/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on cmig_compile_mmps
%   Created:  11/20/2025 by Don Hagler
%   Last Mod: 11/20/2025 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

repos = getenv('CMIG_REPOS_DIR');
if isempty(repos)
  % the folder holding this tree (<repos>/<tree>/matlab/cmig_utils/)
  repos = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
parms = mmil_args2parms(varargin,{...
  'funcname','cmig_run_bids_module',[],...
  'pipeline_dir',[repos '/cmig_pipeline_dev'],[],...
  'toolbox_dir',[repos '/surface_morphometry_toolbox'],[],...
  'gradient_dir',[repos '/gradient_distortion_fields'],[],...
  'toolbox_data_dir',[],[],...
  'matlab_funcdir',[],[],...
  'outdir',[pwd '/output'],[],...
});

if isempty(parms.toolbox_data_dir), error('toolbox_data_dir not given'); end
for d = {parms.pipeline_dir,parms.toolbox_dir,parms.gradient_dir,parms.toolbox_data_dir}
  if ~exist(d{1},'dir'), error('directory %s not found',d{1}); end
end
if ~isempty(parms.matlab_funcdir) && ~exist(parms.matlab_funcdir,'dir')
  error('matlab_funcdir %s not found',parms.matlab_funcdir);
end

% only these trees on the path, so a dependency missing from them fails
%   the compile instead of being taken from MMPS
restoredefaultpath;
add_tree(parms.toolbox_dir);
addpath(parms.gradient_dir);
add_tree([parms.pipeline_dir '/matlab']);
if ~isempty(parms.matlab_funcdir)
  add_tree(parms.matlab_funcdir);
end

if exist(parms.funcname,'file') == 2 && ~isempty(regexp(parms.funcname,'^/','once'))
  fname_func = parms.funcname;
else
  fname_func = which(parms.funcname);
end
if isempty(fname_func) || ~exist(fname_func,'file')
  error('file for funcname %s not found',parms.funcname);
end

% the toolbox data files the code loads by name (load('...'), load ...)
data_files = {...
  'MorphometryToolbox/2_Viewers/tkmcolors.mat',...
  'MorphometryToolbox/4_Segmentation/atlases/batl_all_rbm.mat',...
  'MorphometryToolbox/4_Segmentation/atlases/bm_hatl_afm.mat',...
  'MorphometryToolbox/4_Segmentation/atlases/ctxatl_1_4_rbm.mat',...
  'MorphometryToolbox/4_Segmentation/atlases/hatl_all_rbm.mat',...
  'MorphometryToolbox/7_SkullStripping/SkullStripAMD/SkullStripAMDparams.mat'};
extra_files = [{...
  which('calmi_hatl'),...
  sprintf('%s/SurfaceGeometry/1_IO/offFiles/icosahedron.off',parms.toolbox_dir)},...
  strcat(parms.toolbox_data_dir,'/',data_files)];
for i = 1:numel(extra_files)
  if isempty(extra_files{i}) || ~exist(extra_files{i},'file')
    error('additional file %s not found',extra_files{i});
  end
end

fprintf('%s: creating output dir %s...\n',mfilename,parms.outdir);
mmil_mkdir(parms.outdir);

fprintf('%s: compiling matlab code...\n',mfilename);
orig_state = warning;
warning('off','all');
compiler.build.standaloneApplication(fname_func,...
  'AdditionalFiles',extra_files,'OutputDir',parms.outdir);
warning(orig_state);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function add_tree(topdir)
  % genpath of topdir without .git folders
  dirlist = regexp(genpath(topdir),pathsep,'split');
  dirlist = dirlist(~cellfun(@isempty,dirlist) & ...
                    cellfun(@isempty,regexp(dirlist,'/\.git(/|$)')));
  addpath(dirlist{:});
return;
