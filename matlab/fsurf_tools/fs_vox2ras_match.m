function [voxmatch,rotmatch,offmatch,maxrot,maxoff] = fs_vox2ras_match(M1,M2,max_rotation,max_offset)
%function [voxmatch,rotmatch,offmatch,maxrot,maxoff] = fs_vox2ras_match(M1,M2,max_rotation,max_offset)
%
% Purpose: compare two vox2ras matrices
%
% Required Input:
%   M1: 4x4 vox2ras matrix
%   M2: 4x4 vox2ras matrix
%
% Optional Input:
%   max_rotation: maximum relative rotation (degrees clockwise)
%     {default = 1e-3}
%   max_offset: maximum relative offset (mm)
%     {default = 20}
%
% Output:
%   voxmatch: [0|1] whether all components match
%   rotmatch: [0|1] whether rotations match
%   offmatch: [0|1] whether offsets match (within max_offset)
%   maxrot: maximum rotation difference absolute value
%   maxoff: maximum offset difference absolute value
%
% Created:  06/25/20 by Don Hagler
% Prev Mod: 06/25/20 by Don Hagler
% Last Mod: 02/04/21 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('max_rotation','var') || isempty(max_rotation)
  max_rotation = 1e-3;
end;
if ~exist('max_offset','var') || isempty(max_offset)
  max_offset = 20;
end;

voxmatch = 0; rotmatch = 0; offmatch = 0; maxrot = 0; maxoff = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vec: [transx, transy, transz, rotx, roty, rotz]
%      translation in mm and rotation in degrees clockwise
%      rotations would be applied in order [x,y,z]
vec1 = mmil_M_mat2vec(M1);
vec2 = mmil_M_mat2vec(M2);

offset = vec1(1:3) - vec2(1:3);
rotation = vec1(4:6) - vec2(4:6);
for r=1:length(rotation)
  if rotation(r)>180, rotation(r) = rotation(r) - 360; end;
  if rotation(r)<-180, rotation(r) = rotation(r) + 360; end;
end;

maxrot = max(abs(rotation(:)));
maxoff = max(abs(offset(:)));

if all(abs(M1(:)-M2(:))<=max_rotation), voxmatch = 1; end;
if maxrot <= max_rotation, rotmatch = 1; end;
if maxoff <= max_offset, offmatch = 1; end;

