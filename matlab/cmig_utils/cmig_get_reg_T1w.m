function M_reg = cmig_get_reg_T1w(json)
%function M_reg = cmig_get_reg_T1w(json)
%
% Purpose: read the registration matrix to T1w from a processed json sidecar,
%   accepting either the current field name or the legacy one
%
%   RegistrationMatrixT1w is the current name. registration_matrix_T1 is the
%   name used before 08/2026, and is the only cmig-added field that shipped in
%   released mproc derivatives, which cannot be reprocessed. Reading both keeps
%   those derivatives usable without writing the old name alongside the new one
%
% Required Input:
%   json: struct loaded from a processed json sidecar
%
% Output:
%   M_reg: 4x4 registration matrix; empty if neither field is present
%
% Created:  08/07/2026 by Don Hagler
% Last Mod: 08/07/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_reg = [];
if ~mmil_check_nargs(nargin,1), return; end
if isempty(json) || ~isstruct(json), return; end

% current name first, legacy name second
M_reg = mmil_getfield(json,'RegistrationMatrixT1w',[]);
if isempty(M_reg)
  M_reg = mmil_getfield(json,'registration_matrix_T1',[]);
end

return;
