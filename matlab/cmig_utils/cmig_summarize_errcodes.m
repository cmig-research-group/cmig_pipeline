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
%   errcode: 1 if all scans failed or a failed scan's message is a resource
%     error (cmig_is_resource_error), 0 otherwise
%   errmsg:  non-empty error messages, joined with '; '
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,2), return; end

for i=1:length(errmsgs)
  errmsg = cmig_append_errmsg(errmsg,errmsgs{i});
end
if all(errcodes), errcode = 1; end
% a resource error (e.g. out of memory) fails the module even when other
%   scans succeeded, so the job is rerun
for i=1:length(errmsgs)
  if errcodes(i) && ~isempty(errmsgs{i}) && cmig_is_resource_error(errmsgs{i})
    errcode = 1;
  end
end
