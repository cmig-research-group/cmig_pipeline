function [qcinfo,errcode,errmsg] = cmig_detect_brain_cutoff(qcinfo, fnames, fnames_info, fnamestem, varargin)
%function [qcinfo,errcode,errmsg] = cmig_detect_brain_cutoff(qcinfo, fnames, fnames_info, fnamestem, [options])
%
% Purpose: detect brain cutoff across a set of EPI scans
%   uses FreeSurfer aseg brain mask and EPI-to-T1 registration
%
% Required Input:
%   qcinfo: session-level QC info struct to populate
%   fnames: cell array of EPI file names
%   fnames_info: cell array of proc JSON sidecar file names (may be empty)
%   fnamestem: string identifier for this modality (e.g. 'func', 'DTI')
%
% Optional Input ('key',value):
%   'tmpdir': temporary/working directory
%     {default = ''}
%   'outdir': output directory for plots
%     {default = ''}
%   'outstem': output file stem for plot filenames
%     {default = ''}
%   'fspath': FreeSurfer subject directory
%     {default = ''}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   qcinfo: updated QC info struct with brain cutoff metrics
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 10/01/2026 by Don Hagler
% Last Mod: 11/02/2026 by Don Hagler
%

errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,4), errcode = 1; return; end
parms = mmil_args2parms(varargin,{...
  'tmpdir','',[],...
  'outdir','',[],...
  'outstem','',[],...
  'fname_plot','',[],...
  'fspath','',[],...
  'verbose',false,[],...
  'forceflag',false,[],...
});

checkname = sprintf('%s brain cutoff',fnamestem);

if isempty(parms.fspath)
  fprintf('%s: WARNING: skipping %s: no fspath\n',mfilename,checkname);
  return;
end
fname_t1 = sprintf('%s/mri/nu.mgz',parms.fspath);
if ~exist(fname_t1,'file')
  fprintf('%s: WARNING: skipping %s: %s not found\n',mfilename,checkname,fname_t1);
  return;
end
fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
if ~exist(fname_aseg,'file')
  fprintf('%s: WARNING: skipping %s: %s not found\n',mfilename,checkname,fname_aseg);
  return;
end

% get registration_matrix_T1 (T1->EPI) from first available proc JSON
M_T1_to_EPI = [];
for f=1:length(fnames)
  if ~isempty(fnames_info) && iscell(fnames_info) && length(fnames_info)>=f
    fname_info_f = fnames_info{f};
  else
    fname_info_f = [];
  end
  if isempty(fname_info_f) || ~exist(fname_info_f,'file')
    fname_info_f = regexprep(fnames{f},'\.nii.*','.json');
  end
  if ~exist(fname_info_f,'file'), continue; end
  try
    json = jsondecode(fileread(fname_info_f));
    M_reg = cmig_get_reg_T1w(json);
    if ~isempty(M_reg)
      M_T1_to_EPI = M_reg;
      break;
    end
  catch
  end
end
if isempty(M_T1_to_EPI)
  fprintf('%s: WARNING: skipping %s: registration_matrix_T1 not found in proc JSON\n',mfilename,checkname);
  return;
end

fname_cache = sprintf('%s/%s_brain_cutoff.mat',parms.tmpdir,fnamestem);
if ~exist(fname_cache,'file') || parms.forceflag
  % get T1 (FreeSurfer nu.mgz) vox2ras and dimensions
  [M_T1,volsz_T1] = mmil_load_mgh_info(fname_t1,parms.forceflag,parms.tmpdir);

  % get EPI vox2ras and dimensions from first proc NIfTI
  bc_tmpdir = sprintf('%s/%s_brain_cutoff',parms.tmpdir,fnamestem);
  mmil_mkdir(bc_tmpdir);
  fname_epi_mgz = cmig_import_to_mgz(fnames{1},bc_tmpdir);
  if isempty(fname_epi_mgz)
    errmsg = sprintf('failed to import EPI file: %s',fnames{1});
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end
  [M_T2,volsz_T2] = mmil_load_mgh_info(fname_epi_mgz,parms.forceflag,bc_tmpdir);

  % build RegInfo
  RegInfo = [];
  RegInfo.M_T1       = M_T1;
  RegInfo.volsz_T1   = volsz_T1(1:3);
  RegInfo.M_T2       = M_T2;
  RegInfo.volsz_T2   = volsz_T2(1:3);
  RegInfo.M_T1_to_T2 = M_T1_to_EPI;

  % create aseg brain mask in T1 space
  % brain mask carrying the T1 intensities (brain_flag 1), so the plot shows
  %   the brain; the cutoff counts its nonzero voxels as a binary mask would
  fname_brain = sprintf('%s/aseg_brain.mgz',parms.tmpdir);
  if ~exist(fname_brain,'file') || parms.forceflag
    mmil_aseg_brainmask(parms.fspath,'fname_mask',fname_brain,'brain_flag',1,...
      'forceflag',parms.forceflag);
  end

  if ~isempty(parms.fname_plot)
    fname_plot = parms.fname_plot;
  else
    fname_plot = sprintf('%s/%s_%s_brain_cutoff.jpg',parms.outdir,parms.outstem,fnamestem);
  end
  if parms.verbose
    fprintf('%s: detecting %s brain cutoff...\n',mfilename,fnamestem);
    tic;
  end
  % the EPI's nonzero voxels mark the imaged region, so the zero padding of
  %   a resampled volume counts as outside it
  [cutoff,dorsal,ventral] = mmil_brain_cutoff(RegInfo,fname_brain,...
    'fname_T2',fname_epi_mgz,'fname_plot',fname_plot,...
    'fspath',parms.fspath,'forceflag',parms.forceflag);
  if parms.verbose, toc; end
  save(fname_cache,'cutoff','dorsal','ventral');
else
  output = mmil_load_matfile(fname_cache,{'cutoff','dorsal','ventral'});
  cutoff  = output.cutoff;
  dorsal  = output.dorsal;
  ventral = output.ventral;
end

qcinfo.brain_cutoff         = cutoff;
qcinfo.brain_cutoff_dorsal  = dorsal;
qcinfo.brain_cutoff_ventral = ventral;

return;
