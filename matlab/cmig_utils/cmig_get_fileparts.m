function [fpath,fstem,fext] = cmig_get_fileparts(fname)
%function [fpath,fstem,fext] = cmig_get_fileparts(fname)
%
% Required Input:
%   fname: file name string (may include .nii, .nii.gz, etc.)
%
% Output:
%   fpath: file path
%   fstem: file stem (for .nii* files, stem retains any non-nii suffix)
%   fext:  file extension (.nii, .nii.gz, etc.)
%
% Created:  05/19/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  k = regexp(fname,'(?<ext>.nii.*)','names');
  if ~isempty(k)
    fname_tmp = regexprep(fname,'.nii.*','');
    [fpath,fstem,fext] = fileparts(fname_tmp);
    fstem = [fstem fext];
    fext = k.ext;
  else
    [fpath,fstem,fext] = fileparts(fname);
  end
return;
