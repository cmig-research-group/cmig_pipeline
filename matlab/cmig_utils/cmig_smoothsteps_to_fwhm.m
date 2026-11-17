function fwhm = cmig_smoothsteps_to_fwhm(smoothsteps)
%function fwhm = cmig_smoothsteps_to_fwhm(smoothsteps)
%
% Purpose: convert surface smoothing steps to the equivalent full width at
%   half maximum, in mm
%
% Required Input:
%   smoothsteps: number of smoothing steps (scalar or vector)
%
% Output:
%   fwhm: full width at half maximum, in mm
%
% Created:  08/28/2026 by Don Hagler
% Last Mod: 08/28/2026 by Don Hagler
%

fwhm = [];
if ~mmil_check_nargs(nargin,1), return; end

fwhm = 1.25 * sqrt(smoothsteps);

return;
