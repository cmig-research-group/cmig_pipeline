function val = cmig_getenv(name)
%function val = cmig_getenv(name)
%
% Purpose: return the value of an environment variable the pipeline needs,
%   with an error naming the variable if it is not set
%
% Required Input:
%   name: name of the environment variable
%     e.g. 'CMIG_PIPELINE_DIR', 'CMIG_PIPELINE_ATLASDIR'
%
% Output:
%   val: value of the environment variable
%
% Created:  10/26/2026 by Don Hagler
% Last Mod: 10/26/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

val = getenv(name);
if isempty(val)
  error('environment variable %s not set',name);
end;

return;
