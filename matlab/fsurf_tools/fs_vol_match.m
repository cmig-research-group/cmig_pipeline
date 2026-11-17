function [dimsmatch,orientmatch,voxmatch,rotmatch,offmatch,M1,M2,volsz1,volsz2] = fs_vol_match(fname1,fname2,outdir,forceflag,max_rotation,max_offset)
%function [dimsmatch,orientmatch,voxmatch,rotmatch,offmatch,M1,M2,volsz1,volsz2] = fs_vol_match(fname1,fname2,[outdir],[forceflag],[max_rotation],[max_offset])
%
% Purpose: read header of mgh files and check if volume dimensions match
%
% Required Input:
%   fname1: full path name of first mgh/mgz file
%   fname2: full path name of second mgh/mgz file
%
% Optional Input:
%   outdir: output directory for info mat files
%     if empty, will not save mat files
%     {default = []}
%   forceflag: [0|1] overwrite existing info mat files
%     {default = 0}
%   max_rotation: maximum value in M1(1:3,1:3) - M2(1:3,1:3)
%     {default = 1e-3}
%   max_offset: maximum value in M1(1:3,4) - M2(1:3,4)
%     {default = 10}
%
% Output:
%   dimsmatch: [0|1] whether first three dimensions match
%   orientmatch: [0|1] whether volume orientations match
%   voxmatch: [0|1] whether vox2ras matrices match
%   rotmatch: [0|1] whether rotation components of vox2ras matrices match
%   offmatch: [0|1] whether offset components of vox2ras matrices match
%   M1: vox2ras matrix for fname1
%   M2: vox2ras matrix for fname2
%   volsz1: volume dimensions for fname1
%   volsz2: volume dimensions for fname2
%
% Created:  02/08/11 by Don Hagler
% Prev Mod: 03/09/20 by Don Hagler
% Last Mod: 06/29/20 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('outdir','var'), outdir = []; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end;
if ~exist('max_rotation','var') || isempty(max_rotation)
  max_rotation = 1e-3;
end;
if ~exist('max_offset','var') || isempty(max_offset)
  max_offset = 20;
end;

dimsmatch = 0;
orientmatch = 0;
voxmatch = 0;
rotmatch = 0;
offmatch = 0;
M1=[]; M2=[];
volsz1=[]; volsz2=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that input volumes share dimensions
[M1,volsz1] = mmil_load_mgh_info(fname1,forceflag,outdir);
[M2,volsz2] = mmil_load_mgh_info(fname2,forceflag,outdir);
diff_volsz = volsz1(1:3) - volsz2(1:3);
if all(diff_volsz==0), dimsmatch = 1; end;

if nargout>1
  orient1 = fs_read_orient([],M1);
  orient2 = fs_read_orient([],M2);
  if strcmp(orient1,orient2), orientmatch = 1; end;
  [voxmatch,rotmatch,offmatch,maxrot,maxoff] = fs_vox2ras_match(M1,M2,max_rotation,max_offset);
end;

