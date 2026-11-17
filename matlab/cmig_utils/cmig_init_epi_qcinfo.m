function auto_qcinfo = cmig_init_epi_qcinfo(fname, stype, participant_id, session_id)
%function auto_qcinfo = cmig_init_epi_qcinfo(fname, stype, participant_id, session_id)
%
% Purpose: initialize QC info struct for an EPI (func or dwi) scan
%
% Required Input:
%   fname: input EPI file name (used to derive fstem)
%   stype: scan type string (e.g. 'func', 'dwi')
%   participant_id: participant identifier string (e.g. 'sub-XXXXXXXX')
%   session_id: session identifier string (e.g. 'ses-00A')
%
% Output:
%   auto_qcinfo: initialized QC info struct with NaN/empty placeholders
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Last Mod: 10/05/2026 by Don Hagler
%

auto_qcinfo = [];
if ~mmil_check_nargs(nargin,4), return; end

[~,fstem] = cmig_get_fileparts(fname);
auto_qcinfo.VisitID                = sprintf('%s_%s',participant_id,session_id);
auto_qcinfo.stype                  = stype;
auto_qcinfo.fstem                  = fstem;
auto_qcinfo.StudyDate              = '';
auto_qcinfo.Manufacturer           = '';
auto_qcinfo.ManufacturersModelName = '';
auto_qcinfo.DeviceSerialNumber     = '';
auto_qcinfo.PhaseEncodingDirection = '';
auto_qcinfo.RepetitionTime         = NaN;
auto_qcinfo.registration_matrix_T1 = [];
auto_qcinfo.volsz_raw              = [];
auto_qcinfo.volsz_proc             = [];
auto_qcinfo.regT1_corr             = NaN;
auto_qcinfo.regT1_norm             = NaN;
auto_qcinfo.regT1_warp             = NaN;
auto_qcinfo.regT1_imgcorr1         = NaN;
auto_qcinfo.regT1_imgcorr2         = NaN;
auto_qcinfo.regT1_imgdist3         = NaN;
auto_qcinfo.regT1_imgdist4         = NaN;
auto_qcinfo.regT1_cost0            = NaN;

return;
