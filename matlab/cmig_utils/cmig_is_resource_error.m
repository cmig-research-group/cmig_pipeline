function flag = cmig_is_resource_error(me)
%function flag = cmig_is_resource_error(me)
%
% Purpose: Check whether an error is a resource failure (out of memory,
%   I/O error, disk full, quota exceeded or stale file handle) rather than
%   a problem with the data
%
% Required Input:
%   me: MException, or an error message string
%
% Output:
%   flag: true for a resource error, false otherwise
%
% Notes:
%   The message patterns follow the rules rs_memory, rs_disk and rs_io in
%   parms/errors_rules.tsv (keep the two consistent)
%
% Created:  11/20/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

flag = false;
if ~exist('me','var') || isempty(me), return; end

ids = {'MATLAB:nomem','MATLAB:array:SizeLimitExceeded'};
patterns = {...
  '[Oo]ut of [Mm]emory',...
  'exceeds maximum array size',...
  'Cannot allocate memory',...
  'std::bad_alloc',...
  'No space left on device',...
  'Disk quota exceeded',...
  'Input/output error',...
  'Stale file handle'};

if isa(me,'MException')
  if ismember(me.identifier,ids)
    flag = true;
    return;
  end
  msg = me.message;
elseif ischar(me)
  msg = me;
else
  return;
end
if any(~cellfun(@isempty,regexp(msg,patterns,'once')))
  flag = true;
end
