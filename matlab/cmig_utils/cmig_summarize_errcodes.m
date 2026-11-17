function [errcode,errmsg] = cmig_summarize_errcodes(errcodes,errmsgs)
%function [errcode,errmsg] = cmig_summarize_errcodes(errcodes,errmsgs)
%
% Purpose: Combine per-scan error codes and messages into a single result
%
% Required Input:
%   errcodes: vector of error codes (0 = success, 1 = error) for each scan
%   errmsgs:  cell array of error message strings for each scan
%
% Output:
%   errcode: 1 if all scans failed, 0 otherwise
%   errmsg:  concatenated non-empty error messages
%
% Created:  05/21/2026 by Don Hagler
% Last Mod: 05/21/2026 by Don Hagler
%

errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,2), return; end

if any(~cellfun(@isempty,errmsgs))
  errmsg = strtrim(sprintf('%s ',errmsgs{:}));
end
if all(errcodes), errcode = 1; end
