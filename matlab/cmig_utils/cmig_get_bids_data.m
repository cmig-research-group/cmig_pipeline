function [fnames_data,fnames_info,scan_info] = cmig_get_bids_data(participant_id,session_id,rootdir,modality,series_type,task,acq_type,ignore_scans_flag)
%function [fnames_data,fnames_info,scan_info] = cmig_get_bids_data(participant_id,session_id,rootdir,modality,series_type,task,acq_type,ignore_scans_flag)
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   rootdir = BIDS compliant root directory
%
% Optional Input:
%   modality: modality subdirectory string
%     e.g., 'anat', 'dwi', 'func', 'fmap'
%     if empty, return file names and info for all modalities
%     {default = []}
%   series_type: series type string
%     e.g., 'T1w', 'T2w', 'dwi', 'bold', 'epi'
%     if empty, return files for all series types
%     {default = []}
%   task: task string
%     e.g., 'rest', 'mid', 'nback', 'sst'
%     if empty, return files for all tasks
%     {default = []}
%   acq_type: acquisition type string
%     e.g., 'func', 'dwi'
%     if empty, return files for all acquisition types
%     {default = []}
%   ignore_scans_flag: [0|1] ignore the scans.tsv file, so that
%     scans it marks unusable are still returned
%     {default = 0}
%
% Output:
%   fnames_data: cell array of data file names (nii or nii.gz format)
%   fnames_info: cell array of info file name (json format)
%   scan_info: struct array of scan information
%
% Created:  08/01/2025 by Don Hagler
% Prev Mod: 08/04/2025 by Don Hagler
% Prev Mod: 09/04/2025 by Don Hagler
% Prev Mod: 09/25/2025 by Don Hagler
% Prev Mod: 03/11/2026 by Don Hagler
% Prev Mod: 06/19/2026 by Don Hagler
% Prev Mod: 08/10/2026 by Don Hagler
% Prev Mod: 08/20/2026 by Don Hagler
% Prev Mod: 09/22/2026 by Don Hagler
% Last Mod: 09/25/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames_data = []; fnames_info = []; scan_info = [];
if ~mmil_check_nargs(nargin,3), return; end

if ~exist('modality','var'), modality = []; end
if ~exist('series_type','var'), series_type = []; end
if ~exist('task','var'), task = []; end
if ~exist('acq_type','var'), acq_type = []; end
if ~exist('ignore_scans_flag','var') || isempty(ignore_scans_flag), ignore_scans_flag = 0; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sessdir = sprintf('%s/%s/%s',rootdir,participant_id,session_id);
if ~exist(sessdir,'dir')
  fprintf('%s: WARNING: session dir %s not found\n',mfilename,sessdir);
  return;
end

fname_scan_info = sprintf('%s/%s_%s_scans.tsv',sessdir,participant_id,session_id);
scan_info_flag = 0;
% NOTE: ignore_scans_flag takes the same path as a missing scans.tsv, leaving
%   scan_info and fname_scan empty so that nothing is filtered out. Raw QC is
%   what determines usability, so it has to see the scans marked unusable
if ~exist(fname_scan_info,'file') || ignore_scans_flag
%  fprintf('%s: WARNING: scans file %s not found\n',mfilename,fname_scan_info);
  scan_info = [];
  fname_scan = [];
else
  scan_info_flag = 1;
  scan_info = mmil_csv2struct(fname_scan_info,'separator','\t');
  if isfield(scan_info,'usable')
    idx_usable = find([scan_info.usable]==1);
    scan_info = scan_info(idx_usable);
  else
    fprintf('%s: WARNING: no usable column in %s; treating all scans as usable\n',...
      mfilename,fname_scan_info);
  end
  fname_scan = strcat(sessdir,'/',{scan_info.filename});
end

% find scans that meet criteria of input parameters
if isempty(modality)
  indir = sessdir;
else
  indir = sprintf('%s/%s',sessdir,modality);
end

% set pattern for regular expression
pat = [];
if ~isempty(acq_type), pat = sprintf('%s.*acq-%s',pat,acq_type); end
if ~isempty(task), pat = sprintf('%s.*task-%s',pat,task); end
if ~isempty(series_type), pat = sprintf('%s.*%s',pat,series_type); end
mask_pat = sprintf('%s.*_mask\\.nii',pat);
pat = sprintf('%s.*\\.nii',pat);

% find file names that match pattern
% modality specified: search files directly in indir (sessdir/modality/)
% modality empty: search files in each immediate subdir of sessdir
% never recurse deeper than one modality directory level
fnames_data = {};
if isempty(modality)
  d = dir(indir);
  for i = 1:length(d)
    if d(i).isdir && ~ismember(d(i).name,{'.','..'})
      subdir = fullfile(indir,d(i).name);
      d2 = dir(subdir);
      for j = 1:length(d2)
        if ~d2(j).isdir && ~isempty(regexp(d2(j).name,pat,'once'))
          fnames_data{end+1} = fullfile(subdir,d2(j).name);
        end
      end
    end
  end
else
  d = dir(indir);
  for i = 1:length(d)
    if ~d(i).isdir && ~isempty(regexp(d(i).name,pat,'once'))
      fnames_data{end+1} = fullfile(indir,d(i).name);
    end
  end
end
if isempty(fnames_data)
  fprintf('%s: WARNING: no scans in %s match pattern %s\n',mfilename,indir,pat);
  return;
end

% exclude mask files unless that is all there is
idx_notmask = find(cellfun(@isempty,regexp(fnames_data,mask_pat)));
if ~isempty(idx_notmask)
  fnames_data = fnames_data(idx_notmask);
end

% exclude unusable files
% NOTE: gated on whether scans.tsv was read, not on fname_scan being non-empty,
%   so that a session with no usable scans excludes everything rather than nothing
if scan_info_flag
  [fnames_data,~,idx_scan] = intersect(fnames_data,fname_scan);
  if isempty(fnames_data)
    fprintf('%s: WARNING: no usable scans in session dir %s\n',mfilename,sessdir);
    return;
  end
  scan_info = scan_info(idx_scan);
end

% get correponding json files
fnames_info = cell(size(fnames_data));
for f=1:length(fnames_data)
  fname_image = fnames_data{f};
  fname_json = regexprep(fname_image,'\.nii.+','.json');
  if ~exist(fname_json,'file')
    fprintf('%s: WARNING: no json file found for %s\n',mfilename,fname_image);
  else
    fnames_info{f} = fname_json;
  end
end

return;


