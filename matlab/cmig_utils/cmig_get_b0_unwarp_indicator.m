function b0_unwarp_indicator = cmig_get_b0_unwarp_indicator(fnames_data,fnames_info,proc_rootdir,participant_id,session_id,acq_type)
%function b0_unwarp_indicator = cmig_get_b0_unwarp_indicator(fnames_data,fnames_info,proc_rootdir,participant_id,session_id,acq_type)
%
% Purpose: determine whether B0 unwarp was applied for a session
%
%   Checks each scan's proc JSON sidecar for a B0UnwarpApplied field
%   (written by cmig_proc_bids_func); reports the whole session as
%   unwarped (1) if any scan has it applied, not unwarped (0) if every
%   scan explicitly reports it was not applied
%
%   If no scan's JSON has B0UnwarpApplied (e.g. not yet written by an
%   upstream proc module), falls back to an availability-based check
%   (mirroring source's B0unwarp_flag): whether a B0dx displacement
%   field file exists for this acq_type, regardless of whether unwarp
%   was actually applied
%
% Required Input:
%   fnames_data: cell array of per-scan data file names
%   fnames_info: cell array of per-scan JSON sidecar file names
%     (may contain empty entries; falls back to deriving the JSON name
%     from the corresponding data file name)
%   proc_rootdir: BIDS processed derivatives root directory
%   participant_id: participant identifier string
%   session_id: session identifier string
%   acq_type: acquisition type string used in the B0dx filename's acq-
%     entity (e.g. 'func', 'dwi')
%
% Output:
%   b0_unwarp_indicator: 1 if unwarp applied/available, 0 if not, NaN
%     if indeterminate (no B0UnwarpApplied and no fmap directory)
%
% Created:  07/10/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Last Mod: 10/01/2026 by Don Hagler
%

b0_unwarp_indicator = NaN;
if ~mmil_check_nargs(nargin,6), return; end

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
    val = mmil_getfield(json,'B0UnwarpApplied',NaN);
    if ~isnan(val)
      vals(end+1) = double(val);
    end
  catch
  end
end
if ~isempty(vals)
  % report the whole session as unwarped if any scan had it applied,
  % rather than whichever scan happened to be checked first
  b0_unwarp_indicator = max(vals);
  return;
end

% fall back to an availability-based check (mirroring source's
% B0unwarp_flag) when no scan's JSON has B0UnwarpApplied: were the
% reference scans required for B0 unwarp available at all, regardless
% of whether unwarp was actually applied
fmap_dir = sprintf('%s/%s/%s/fmap',proc_rootdir,participant_id,session_id);
% no fmap directory: the proc tree says nothing about unwarping
if ~exist(fmap_dir,'dir'), return; end
dlist = dir(sprintf('%s/*_acq-%s_desc-AP*PA*_B0dx.nii.gz',fmap_dir,acq_type));
b0_unwarp_indicator = double(~isempty(dlist));

return;
