function options = cmig_parse_deployed_args(options)
%function options = cmig_parse_deployed_args(options)
%
% Purpose: Parse options cell array values from strings to appropriate types
%   when running as a compiled (deployed) application; call before
%   mmil_args2parms to ensure typed valid ranges are not violated by char
%   string inputs from compiled deployment
%
% Required Input:
%   options: cell array of name-value pairs (as passed to mmil_args2parms)
%
% Output:
%   options: cell array with even-indexed values coerced from string to
%     typed values (numeric, cell array, or numeric array) via cmig_parseval
%
% Created:  06/24/2026 by Don Hagler
% Last Mod: 06/24/2026 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end

if isdeployed
  for i=2:2:length(options)
    options{i} = cmig_parseval(options{i});
  end
end
