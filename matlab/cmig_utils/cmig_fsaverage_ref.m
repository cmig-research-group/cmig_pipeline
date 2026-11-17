function [spatial_ref,desc] = cmig_fsaverage_ref(ico_order)
%function [spatial_ref,desc] = cmig_fsaverage_ref(ico_order)
%
% Purpose: SpatialReference and Description text of a surface map on the
%   atlas-ico{N} sphere: fsaverage truncated to its first 10*4^N+2 vertices
%
% Required Input:
%   ico_order: icosahedral order (0-7)
%
% Output:
%   spatial_ref: 'fsaverage5' for ico5, 'fsaverage6' for ico6, else 'fsaverage'
%   desc: e.g. 'fsaverage''s first 642 vertices (ico3)'
%
% NOTE: fsaverage5 and fsaverage6 have the same vertices as fsaverage's first
%   10242 and 40962; fsaverage3 and fsaverage4 order their vertices otherwise
%
% Created:  11/12/2026 by Don Hagler
% Last Mod: 11/12/2026 by Don Hagler
%

spatial_ref = []; desc = [];
if ~mmil_check_nargs(nargin,1), return; end

switch ico_order
  case {5,6}, spatial_ref = sprintf('fsaverage%d',ico_order);
  otherwise,  spatial_ref = 'fsaverage';
end
desc = sprintf('fsaverage''s first %d vertices (ico%d)',10*4^ico_order+2,ico_order);

return;
