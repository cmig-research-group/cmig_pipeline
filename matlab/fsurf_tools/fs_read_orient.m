function [orient,oblique_flag] = fs_read_orient(fname,M)
%function [orient,oblique_flag] = fs_read_orient(fname,[M])
%
% Purpose: read header of mgh file and determine volume orientation
%
% Required Input:
%   fname: full path name of mgh/mgz file
%    If empty, expects M to be supplied
%
% Optional Input:
%   M: vox2ras matrix
%    {default = []}
%
% Output:
%   orient: three character string specifying volume orientation
%     e.g. 'RAS', 'LPI', 'PRI', etc.
%   oblique_flag: [0|1] whether slices are oblique
%
% Created:  02/14/10 by Don Hagler
% Prev Mod: 05/03/10 by Don Hagler
% Prev Mod: 12/12/18 by Don Hagler
% Last Mod: 04/10/25 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
orient = 'XXX';
oblique_flag = 0;
smf = 10*eps;

if ~exist('M','var') || isempty(M)
  [vol,M] = fs_load_mgh(fname,[],[],1);
end;

M = M(1:3,1:3);
voxsize = sqrt(sum(abs(M).^2,1));

Rvec = [1 0 0]*M./voxsize;
Avec = [0 1 0]*M./voxsize;
Svec = [0 0 1]*M./voxsize;

% method below based on FreeSurfer's MRIdircosToOrientationString mri.c by Bruce Fischl
for c=1:3
  R = Rvec(c);
  A = Avec(c);
  S = Svec(c);
  if abs(R) > abs(A) && abs(R) > abs(S)
    if (R > 0)
      orient(c) = 'R';
    else
      orient(c) = 'L';
    end
  elseif abs(A) > abs(S)
    if (A > 0)
      orient(c) = 'A';
    else
      orient(c) = 'P';
    end
  else
    if (S > 0)
      orient(c) = 'S';
    else
      orient(c) = 'I';
    end
  end
end

if nargout>1
  if length(find(abs(Rvec)>smf))>1 ||...
    length(find(abs(Avec)>smf))>1 ||...
    length(find(abs(Svec)>smf))>1
    oblique_flag = 1;
  end
end


