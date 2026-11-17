function [fs_version,fs_build] = cmig_get_fs_ver(fspath)
%function [fs_version,fs_build] = cmig_get_fs_ver(fspath)
%
% Optional Input:
%   fspath: full path of FreeSurfer output directory
%     if supplied, will return version derived from {fspath}/scripts/build_stamp.txt
%     if not, or if build stamp not found, will return FREESURFER_VER environment variable
%     
% Output:
%   fs_version: FreeSurfer version number
%     e.g., 450 corresponds to FreeSurfer version 4.5.0
%   fs_build: FreeSurfer build stamp
%     e.g., freesurfer-linux-centos7_x86_64-7.1.1-20200723-8b40551
%
% Created:  11/04/2025 by Don Hagler
% Last Mod: 04/24/2026 by Don Hagler
%

fs_version = []; fs_build = [];
if ~exist('fspath','var'), fspath = []; end

% get fs_version from build stamp
if ~isempty(fspath)
  if exist(fspath,'dir')
    fname_build = sprintf('%s/scripts/build-stamp.txt',fspath);
    if exist(fname_build,'file')
      try
        tmp = mmil_readtext(fname_build);
        fs_build = tmp{1};
        % e.g., freesurfer-linux-centos7_x86_64-7.1.1-20200723-8b40551
        k = regexp(fs_build,'freesurfer-.+-(?<version>\d\.\d\.\d)-','names');
        if ~isempty(k)
          fs_version = str2num(regexprep(k.version,'\.',''));
        else
          fprintf('%s: WARNING: failed to determine FreeSurfer version from build stamp\n',mfilename);
        end
      catch
        fprintf('%s: WARNING: failed to read build stamp from %s\n',mfilename,fname_build);
      end
    else
      fprintf('%s: WARNING: unable to find FreeSurfer build stamp file %s\n',mfilename,fname_build);
    end
  else
    fprintf('%s: WARNING: fspath %s not found\n',mfilename,fspath);
  end
end

% set FreeSurfer version if not already set
if isempty(fs_version)
  fs_version = str2num(getenv('FREESURFER_VER'));
end

return;


