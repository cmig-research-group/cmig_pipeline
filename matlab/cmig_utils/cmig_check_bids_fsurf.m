function [fspath,fs_version,fs_build,fs_complete_flag,errcode,errmsg] = cmig_check_bids_fsurf(participant_id,session_id,fsurf_rootdir)
%function [fspath,fs_version,fs_build,fs_complete_flag,errcode,errmsg] = cmig_check_bids_fsurf(participant_id,session_id,fsurf_rootdir)
%
% Purpose: locate and validate a FreeSurfer recon directory for a
%          given subject/session; return path and version info
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   fsurf_rootdir: root directory containing FreeSurfer recon directories
%
% Output:
%   fspath: full path to the FreeSurfer recon directory
%   fs_version: FreeSurfer major version number
%     empty if the recon directory is missing
%   fs_build: FreeSurfer build stamp string
%     empty if the recon directory is missing
%   fs_complete_flag: [0|1] whether surface outputs (not just volume
%     outputs) are present; 0 only when recon is "volume-only" complete
%     (mri/aseg.mgz exists but surface outputs do not), matching
%     source's MMIL_Get_FSReconStatus status==6 - errcode is NOT set
%     in this case, allowing the caller to proceed with volume/aseg
%     measures while skipping surface-dependent ones
%   errcode: [0|1] error code (1 = fatal error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Last Mod: 08/11/2026 by Don Hagler
%

errcode = 0; errmsg = [];
fs_version = []; fs_build = []; fs_complete_flag = 0;

fspath = sprintf('%s/%s_%s',fsurf_rootdir,participant_id,session_id);

if ~exist(fspath,'dir')
  errcode = 1;
  errmsg = sprintf('FreeSurfer directory %s not found',fspath);
  return;
end

fprintf('%s: FreeSurfer directory:\n %s\n',mfilename,fspath);

% get version and build stamp before the completeness checks,
%   so they are available even for an incomplete recon
[fs_version,fs_build] = cmig_get_fs_ver(fspath);
fprintf('%s: FreeSurfer build stamp: %s\n',mfilename,fs_build);
fprintf('%s: FreeSurfer version: %d\n',mfilename,fs_version);

fs_complete_flag = cmig_check_fsurf_files(fspath);

if ~fs_complete_flag
  % check for volume-only completeness (mri/aseg.mgz present even though
  % surface outputs are missing), matching source's
  % MMIL_Get_FSReconStatus status==6 ("volume-only recon complete")
  fname_aseg = sprintf('%s/mri/aseg.mgz',fspath);
  if exist(fname_aseg,'file')
    errmsg = sprintf('freesurfer recon %s is volume-only (surface outputs missing)',fspath);
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
  else
    errcode = 1;
    errmsg = sprintf('freesurfer recon %s is incomplete (missing required files)',fspath);
    return;
  end
end

return;
