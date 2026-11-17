function bids_hemi = cmig_get_bids_hemi(hemi)
%function bids_hemi = cmig_get_bids_hemi(hemi)
%
% Purpose: convert a FreeSurfer hemisphere label to the BIDS hemi entity value
%
% Required Input:
%   hemi: hemisphere label, 'lh' or 'rh'
%
% Output:
%   bids_hemi: 'L' or 'R'; any other input is returned unchanged
%
% Created:  08/28/2026 by Don Hagler
% Last Mod: 08/28/2026 by Don Hagler
%

bids_hemi = [];
if ~mmil_check_nargs(nargin,1), return; end

switch hemi
  case 'lh', bids_hemi = 'L';
  case 'rh', bids_hemi = 'R';
  otherwise,  bids_hemi = hemi;
end

return;
