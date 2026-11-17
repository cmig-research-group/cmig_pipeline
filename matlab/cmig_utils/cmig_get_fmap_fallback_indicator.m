function fmap_fallback_indicator = cmig_get_fmap_fallback_indicator(fnames_data,fnames_info)
%function fmap_fallback_indicator = cmig_get_fmap_fallback_indicator(fnames_data,fnames_info)
%
% Purpose: determine whether a substitute field map pair was used for a session
%
%   Checks each scan's proc JSON sidecar for a B0UnwarpFieldMapFallback field
%   (written by cmig_proc_bids_func); reports the session as having used a
%   substitute (1) if any scan did, not (0) if every scan explicitly reports
%   it did not
%
%   A scan gets a substitute when its B0FieldSource names no available field
%   map pair, e.g. because that pair was rejected by check_fmap_pairs
%
% Required Input:
%   fnames_data: cell array of per-scan data file names
%   fnames_info: cell array of per-scan JSON sidecar file names
%     (may contain empty entries; falls back to deriving the JSON name
%     from the corresponding data file name)
%
% Output:
%   fmap_fallback_indicator: 1 if any scan used a substitute pair, 0 if none
%     did, NaN if no scan's sidecar carries the field
%
% Created:  08/27/2026 by Don Hagler
% Last Mod: 08/27/2026 by Don Hagler
%

fmap_fallback_indicator = NaN;
if ~mmil_check_nargs(nargin,2), return; end

vals = [];
for f=1:length(fnames_data)
  fname_info = [];
  if f<=length(fnames_info), fname_info = fnames_info{f}; end
  if isempty(fname_info) || ~exist(fname_info,'file')
    fname_info = regexprep(fnames_data{f},'\.nii.*','.json');
  end
  if ~exist(fname_info,'file'), continue; end
  try
    json = jsondecode(fileread(fname_info));
    val = mmil_getfield(json,'B0UnwarpFieldMapFallback',NaN);
    if ~isnan(val)
      vals(end+1) = double(val);
    end
  catch
  end
end
if isempty(vals), return; end
% report the session as having used a substitute if any scan did
fmap_fallback_indicator = max(vals);

return;
