function cmig_extract_bids_func_tseries(fnames_data,fnames_info,fspath,outdir,varargin)
%function cmig_extract_bids_func_tseries(fnames_data,fnames_info,fspath,outdir,[options])
%
% Purpose: extract fMRI time series for each run in a session
%
% Required Input:
%   fnames_data: cell array of processed fMRI data file names (.nii.gz)
%   fnames_info: cell array of JSON sidecar file names (one per run)
%   fspath: FreeSurfer subject directory
%   outdir: output directory
%
% Optional Input ('key',value):
%   'fs_version': FreeSurfer version number
%     {default = []}
%   'extract_ico_order': ico order(s) for surface resampling
%     if empty, use cmig_extract_fmri_tseries default
%     {default = []}
%   'extract_ico_smooth': number of smoothing steps for each ico order,
%     passed to fs_surf2surf as smooth_in
%     if empty, use cmig_extract_fmri_tseries default
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  05/20/2026 by Don Hagler
% Prev Mod: 05/20/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 09/07/2026 by Don Hagler
% Last Mod: 10/21/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,4), return; end
parms = mmil_args2parms(varargin,{...
  'fs_version',[],[],...
  'extract_ico_order',[],[],...
  'extract_ico_smooth',[],[],...
  'forceflag',false,[],...
});

run_flag = check_extract_output(fnames_data,outdir);
if ~run_flag, return; end

for s=1:length(fnames_data)
  fname_data = fnames_data{s};
  fname_info = fnames_info{s};
  if isempty(fname_data)
    fprintf('%s: no usable data for this scan\n',mfilename);
    continue;
  end
  % load info file
  try
    scan_info = jsondecode(fileread(fname_info));
  catch me
    error('failed to load json file %s:\n%s',fname_info,me.message);
  end
  % get registration matrix
  M_reg = cmig_get_reg_T1w(scan_info);
  if isempty(M_reg)
    error('RegistrationMatrixT1w missing from %s',fname_info);
  end
  % create per-scan directory in outdir
  [~,fstem] = cmig_get_fileparts(fname_data);
  tmpdir = sprintf('%s/%s',outdir,fstem);
  mmil_mkdir(tmpdir)
  % extract time series
  p = [];
  p.fs_version = parms.fs_version;
  if ~isempty(parms.extract_ico_order), p.ico_order = parms.extract_ico_order; end
  if ~isempty(parms.extract_ico_smooth), p.ico_smooth = parms.extract_ico_smooth; end
  p.verbose = 2;
  p.cleanupflag = 0;
  a = mmil_parms2args(p);
  cmig_extract_fmri_tseries(fname_data,fspath,M_reg,outdir,a{:});
  % copy info json to output directory; ensure writable (source may be read-only)
  [~,fstem,fext] = fileparts(fname_info);
  fname_info_out = sprintf('%s/%s%s',outdir,fstem,fext);
  copyfile(fname_info,fname_info_out);
  fileattrib(fname_info_out,'+w');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_flag = check_extract_output(fnames_data,outdir)
  % suffixes of expected output files from cmig_extract_fmri_tseries
  suffixlist = {...
    '-ico3-lh.mgz'...
    '-ico3-rh.mgz'...
    '-ico5-lh.mgz'...
    '-ico5-rh.mgz'...
    '-ico3-mask-lh.mgz'...
    '-ico3-mask-rh.mgz'...
    '-ico5-mask-lh.mgz'...
    '-ico5-mask-rh.mgz'...
    '_aparc_data.mat'...
    '_aseg_data.mat'...
    '_fparc_data.mat'...
    '_roi_data.mat'...
  };
  run_flag = 0;
  % set run_flag=1 if any output file does not exist
  for f=1:length(fnames_data)
    [~,fstem] = cmig_get_fileparts(fnames_data{f});
    for s=1:length(suffixlist)
      fname_out = sprintf('%s/%s%s',outdir,fstem,suffixlist{s});
      if ~exist(fname_out,'file'), run_flag = 1; return; end
    end
  end
return;
