function fullmeas = cmig_get_smri_fullmeas(fullmeas)
%function fullmeas = cmig_get_smri_fullmeas(fullmeas)
%
% Purpose: convert a short sMRI surface measure name to the full measure name
%   used in the file names written by cmig_analyze_bids_anat
%
% Required Input:
%   fullmeas: short measure name, e.g. 'thk', 't1_gwc'
%
% Output:
%   fullmeas: full measure name, e.g. 'thickness', 'T1w_contrast_pdist0.2';
%     a name with no mapping is returned unchanged, e.g. 'area', 'sulc'
%
% Created:  08/28/2026 by Don Hagler
% Last Mod: 08/28/2026 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), fullmeas = []; return; end

switch fullmeas
  case 'thk'
    fullmeas = 'thickness';
  case 't1_gwc'
    fullmeas = 'T1w_contrast_pdist0.2';
  case 't2_gwc'
    fullmeas = 'T2w_contrast_pdist0.2';
  case 't1_gm'
    fullmeas = 'T1w_pdist0.2';
  case 't2_gm'
    fullmeas = 'T2w_pdist0.2';
  case 't1_wm'
    fullmeas = 'T1w_pdist-0.2';
  case 't2_wm'
    fullmeas = 'T2w_pdist-0.2';
end

return;
