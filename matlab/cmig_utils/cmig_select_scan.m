function idx = cmig_select_scan(scan_info)
%function idx = cmig_select_scan(scan_info)
%
% Purpose: Select among multiple scans by QC note length then acquisition time
%
% Required Input:
%   scan_info: struct array with fields 'notes' and 'acq_time'
%     notes: cell array of QC note strings for each scan
%     acq_time: acquisition time string (format: 'yyyy-mm-ddTHH:MM:SS')
%
% Output:
%   idx: index of the selected scan
%
% Created:  05/21/2026 by Don Hagler
% Last Mod: 05/21/2026 by Don Hagler
%

idx = 1;
if ~mmil_check_nargs(nargin,1), return; end

if length(scan_info) > 1
  % select scans with shortest QC notes
  note_length = cellfun(@length,{scan_info.notes});
  if any(note_length>0) && length(unique(note_length))>1
    [min_val,idx_min] = min(note_length);
    idx = find(note_length == min_val);
  else
    idx = 1:length(note_length);
  end
  if length(idx) > 1
    % select last scan acquired of those remaining
    dnums = cellfun(@(x) datenum(x,'yyyy-mm-ddTHH:MM:SS'),{scan_info(idx).acq_time});
    [max_val,idx_max] = max(dnums);
    idx = idx(idx_max);
  end
end
