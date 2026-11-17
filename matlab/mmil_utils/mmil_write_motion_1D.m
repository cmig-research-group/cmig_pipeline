function mmil_write_motion_1D(motion_data,fname_motion,varargin)
%function mmil_write_motion_1D(motion_data,fname_motion,varargin)
%
% Purpose: save motion time courses to motion.1D file like AFNI's 3dvolreg
%
% Usage:
%   mmil_write_motion_1D(motion_data,fname_motion,'key1', value1,...)
%
% Required Parameters:
%   motion_data: matrix of motion parameter estimates with size [nframes,npar]
%     order of parameters: dx, dy, dz, rx, ry, rz
%     with units of mm and degrees (rotations are clockwise)
%   fname_motion: output text file containing 6 columns of motion estimates
%     plus index and rmsold and rmsnew columns
%     order of columns: n, rz, rx, ry, dz, dx, dy, rmsold, rmsnew
%     with units of mm and degrees (rotations are counterclockwise)
%     from AFNI's 3dvolreg -help:
%       n     = sub-brick index
%       roll  = rotation about the I-S axis
%       pitch = rotation about the R-L axis
%       yaw   = rotation about the A-P axis
%       dS  = displacement in the Superior direction
%       dL  = displacement in the Left direction
%       dP  = displacement in the Posterior direction
%       rmsold = RMS difference between input brick and base brick
%       rmsnew = RMS difference between output brick and base brick
%
% Optional parameters:
%  'skipTRs': number of initial repetitions to remove
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  02/20/18 Don Hagler
% Last Mod: 02/20/18 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'skipTRs',0,[0 Inf],...
  'forceflag',false,[false true],...
...
  'reorder_vec',[6,4,5,3,1,2],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_motion,'file') || parms.forceflag
  numTRs = size(motion_data,1);
  ncols = size(motion_data,2);
  % change rotations to be counterclockwise
  motion_data(:,4:6) = -motion_data(:,4:6);
  % reorder motion parameter estimates
  %   from dx,dy,dz,rx,ry,rz to rz,rx,ry,dz,dx,dy
  motion_data = motion_data(:,parms.reorder_vec);
  fid = fopen(fname_motion,'w');
  k = 1;
  for t=parms.skipTRs+1:numTRs
    % add index and dummy rmsold and rmsnew columns
    fprintf(fid,'  %d  %s  0.0  0.0\n',k,...
      sprintf('%0.4f  ',squeeze(motion_data(t,:))));
    k = k + 1;
  end;
  fclose(fid);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

