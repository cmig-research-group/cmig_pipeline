function [vol,M,mr_parms,volsz] = mmil_load_mgh(fname,slices,frames,headeronly,keepsingle)
%function [vol,M,mr_parms,volsz] = mmil_load_mgh(fname,slices,frames,headeronly,keepsingle)
%
% Purpose: load mgh/mgz file, delete if corrupt
%
% Required Input:
%   fname: full path of mgh/mgz file
% 
% Optional Input:
%   slices: list of one-based slice numbers to load. All
%     slices are loaded if slices is not specified, or
%     if slices is empty, or if slices(1) <= 0.
%     {default = []}
%   frames: list of one-based frame numbers to load. All
%     frames are loaded if frames is not specified, or
%     if frames is empty, or if frames(1) <= 0.
%     {default = []}
%   headeronly: [0|1] whether to load header info only
%     if 1, vol will be empty
%     {default = 0}
%   keepsingle: [0|1] whether to keep float data as single
%     otherwise convert to double
%     {default = 0}
%
% Output:
%   vol: 4D volume (if surface data, size will be [n,1,1,f])
%   M: 4x4 vox2ras transform such that
%     y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
%     indices are 1-based. 
%     If the input has multiple frames,
%     only the first frame is read.
%   mr_parms: [tr flipangle te ti fov]
%   volsz: size(vol). Helpful when using headeronly as vol is [].
%
% See also: fs_load_mgh
%
% Created:  04/14/23 by Don Hagler
% Last Mod: 05/16/25 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on mmil_load_matfile
%   Created:  05/14/2020 by Don Hagler
%   Last Mod: 04/20/2021 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol = []; M = []; volsz = []; mr_parms = [];
if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('slices','var'), slices = []; end
if ~exist('frames','var'), frames = []; end
if ~exist('headeronly','var'), headeronly = []; end
if ~exist('keepsingle','var'), keepsingle = []; end

errmsg = [];
rmflag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = dir(fname);
if isempty(d)
  errmsg = sprintf('file %s not found',fname);
elseif d.bytes==0
  % delete if mgh/mgz file has file size zero
  errmsg = sprintf('file with size 0 found:\n%s\ndeleting empty file... rerun to recreate and proceed',fname);
  rmflag = 1;
else
  % try to load mgh/mgz file and delete if corrupt
  try
    [vol,M,mr_parms,volsz] = fs_load_mgh(fname,slices,frames,headeronly,keepsingle);
  catch me
    if ~isempty(regexp(me.message,'tried to read'))
      errmsg = sprintf('truncated file found:\n%s\ndeleting truncated file... rerun to recreate and proceed',fname);
      rmflag = 1;
    elseif ~isempty(regexp(me.message,'could not read'))
      errmsg = sprintf('corrupt file found:\n%s\ndeleting corrupt file... rerun to recreate and proceed',fname);
      rmflag = 1;
    else
      errmsg = sprintf('failed to load mgh file:\n%s\n%s\ndeleting bad file... rerun to recreate and proceed',...
        fname,me.message);
      rmflag = 1;
    end
  end
end

if rmflag, delete(fname); end;
if ~isempty(errmsg), error(errmsg); end;


