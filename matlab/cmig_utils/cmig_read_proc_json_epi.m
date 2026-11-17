function [auto_qcinfo,errcode,errmsg] = cmig_read_proc_json_epi(auto_qcinfo, fname, fname_info)
%function [auto_qcinfo,errcode,errmsg] = cmig_read_proc_json_epi(auto_qcinfo, fname, fname_info)
%
% Purpose: read proc JSON sidecar for an EPI scan and populate qcinfo fields
%   JSON written by cmig_proc_bids_func (or equivalent for dwi)
%   contains: original BIDS fields + registration_matrix_T1,
%     vox2ras_matrix_raw, vox2ras_matrix_proc, volsz_raw, volsz_proc
%
% Required Input:
%   auto_qcinfo: QC info struct to populate
%   fname: EPI file name (used to derive JSON path if fname_info is empty)
%   fname_info: proc JSON sidecar file name; if empty, derived from fname
%
% Output:
%   auto_qcinfo: updated QC info struct
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Last Mod: 08/07/2026 by Don Hagler
%

errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,3), errcode = 1; return; end

if isempty(fname_info) || ~exist(fname_info,'file')
  fname_info = regexprep(fname,'\.nii.*','.json');
end
if ~exist(fname_info,'file')
  errmsg = sprintf('proc JSON sidecar not found: %s',fname_info);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  errcode = 1;
  return;
end
try
  json = jsondecode(fileread(fname_info));
catch me
  errmsg = sprintf('failed to read JSON %s: %s',fname_info,me.message);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  errcode = 1;
  return;
end
% acquisition info
auto_qcinfo.StudyDate              = mmil_getfield(json,'AcquisitionDateTime','');
auto_qcinfo.Manufacturer           = mmil_getfield(json,'Manufacturer','');
auto_qcinfo.ManufacturersModelName = mmil_getfield(json,'ManufacturersModelName','');
auto_qcinfo.DeviceSerialNumber     = mmil_getfield(json,'DeviceSerialNumber','');
auto_qcinfo.PhaseEncodingDirection = mmil_getfield(json,'PhaseEncodingDirection','');
auto_qcinfo.RepetitionTime         = mmil_getfield(json,'RepetitionTime',NaN);
% processing metrics from cmig_proc_bids_func
% registration_matrix_T1 = M_T1_to_EPI: transform from T1 (atlas) space to EPI space
auto_qcinfo.registration_matrix_T1 = cmig_get_reg_T1w(json);
volsz_raw  = mmil_getfield(json,'SourceVolumeSize',[]);
volsz_proc = mmil_getfield(json,'VolumeSize',[]);
if ~isempty(volsz_raw),  auto_qcinfo.volsz_raw  = volsz_raw;  end
if ~isempty(volsz_proc), auto_qcinfo.volsz_proc = volsz_proc; end

return;
