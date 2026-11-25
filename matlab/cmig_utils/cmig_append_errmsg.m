function errmsg = cmig_append_errmsg(errmsg,msg)
%function errmsg = cmig_append_errmsg(errmsg,msg)
%
% Purpose: Append a message to an error message, separated by '; '
%
% Required Input:
%   errmsg: existing error message (may be empty)
%   msg: message to append (empty: errmsg is returned unchanged)
%
% Output:
%   errmsg: the joined message
%
% Created:  11/20/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

if ~exist('errmsg','var'), errmsg = []; end
if ~exist('msg','var') || isempty(msg), return; end
msg = strtrim(msg);
if isempty(msg), return; end
if isempty(errmsg)
  errmsg = msg;
elseif isempty(strfind(errmsg,msg))
  errmsg = sprintf('%s; %s',errmsg,msg);
end
