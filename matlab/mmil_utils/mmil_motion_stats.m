function [motion_stats,motion_fd,deriv,motion_fdmat] = mmil_motion_stats(motion_tseries,motion_radius,absflag,thresholds,nodyflag,motion_fd_orig,gradflag)
%function [motion_stats,motion_fd,deriv,motion_fdmat] = mmil_motion_stats(motion_tseries,[motion_radius][absflag],[thresholds],[nodyflag],[motion_fd_orig],[gradflag])
%
% Purpose: calculate statistics from motion time course
%
% Required Parameters:
%   motion_tseries: matrix of motion estimates [nframes,6]
%     with this order: dx,dy,dz,rx,ry,rz
%         also called: dL, dP, dS, pitch, yaw, roll
%                      to match output of mmil_load_motion_1D
%     with translations in mm and rotations in degrees clockwise
%
% Optional parameters:
%  motion_radius: approximate radius of typical head (mm)
%    for calculating distance from an angle of rotation
%    {default = 50}
%  absflag: calculate motion as sum of absolute values of each component
%    otherwise, calculate square root of sum of squares
%    {default = 1}
%  thresholds: vector of motion thresholds (mm)
%    used to calculate percent and time of scan below threshold
%    {default = [0.2,0.3,0.4]}
%  nodyflag: calculate motion ignoring dy (translation in y)
%    {default = 0}
%  motion_fd_orig: vector of framewise displacements in mm [nframes,1]
%    if supplied, used to normalize motion_fd before calculating mean_motion, etc.
%    {default = []}
%  gradflag: calculate framewise displacements using gradient instead of diff
%    {default = 0}
%
% Output:
%   motion_stats: struct containing calculated motion statistics with fields:
%     max_rz:      maximum absolute rotz (relative to baseline) (mm)
%     max_rx:      maximum absolute rotx (relative to baseline) (mm)
%     max_ry:      maximum absolute roty (relative to baseline) (mm)
%     max_dz:      maximum absolute rotz (relative to baseline) (mm)
%     max_dx:      maximum absolute rotx (relative to baseline) (mm)
%     max_dy:      maximum absolute roty (relative to baseline) (mm)
%     mean_trans:  mean frame-to-frame translation (mm/frame)
%     max_trans:   maximum frame-to-frame translation (mm/frame)
%     mean_rot:    mean frame-to-frame rotation (mm/frame)
%     max_rot:     maximum frame-to-frame rotation (mm/frame)
%     mean_motion: mean frame-to-frame motion (including translation and rotation) (mm/frame)
%     std_motion:  standard deviation of frame-to-frame motion (mm/frame)
%     mode_motion: mode of frame-to-frame motion (mm/frame)
%     med_motion:  median of frame-to-frame motion (mm/frame)
%     min_motion:  minimum of  frame-to-frame motion (mm/frame)
%     max_motion:  maximum of  frame-to-frame motion (mm/frame)
%     mean_accel:  mean of derivative of frame-to-frame motion  (mm/frame^2)
%     subthresh_nvols: number of volumes with subthreshold frame-wise displacement
%     subthresh_perc: percent of volumes with subthreshold frame-wise displacement
%     thresholds: vector of motion thresholds
%     nvols: number of volumes
%   motion_fd: vector of framewise displacements in mm [nframes,1]
%   deriv: first derivative of motion time courses [nframes,6]
%   motion_fdmat: matrix of framewise displacements in mm [nframes,6]
%
% Created:  02/12/13 Don Hagler
% Prev Mod: 08/23/17 Dani Cornejo 
% Prev Mod: 10/16/17 Don Hagler
% Prev Mod: 01/25/18 Don Hagler
% Last Mod: 05/05/25 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

motion_stats = [];
motion_fd = [];
deriv = [];
motion_fdmat = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('motion_radius','var') || isempty(motion_radius)
  motion_radius = 50;
end;
if ~exist('absflag','var') || isempty(absflag), absflag = 1; end;
if ~exist('thresholds','var') || isempty(thresholds)
  thresholds = [0.2 0.3 0.4];
end;
if ~exist('nodyflag','var') || isempty(nodyflag), nodyflag = 0; end;
if ~exist('motion_fd_orig','var') || isempty(motion_fd_orig), motion_fd_orig = []; end;
if ~exist('gradflag','var') || isempty(gradflag), gradflag = 0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate scaling factor to convert rotations in degrees
%   to approximate displacements in mm (on sphere about the size of a head)
rot_fact = pi*motion_radius/180;

% calculate maximum translations relative to baseline
motion_stats.max_dx = max(abs(motion_tseries(:,1)));
motion_stats.max_dy = max(abs(motion_tseries(:,2)));
motion_stats.max_dz = max(abs(motion_tseries(:,3)));
motion_stats.max_rx = max(abs(motion_tseries(:,4)))*rot_fact;
motion_stats.max_ry = max(abs(motion_tseries(:,5)))*rot_fact;
motion_stats.max_rz = max(abs(motion_tseries(:,6)))*rot_fact;

% calculate difference in head position or rotation
if size(motion_tseries,1) > 1
  if gradflag
    d = gradient(motion_tseries);
    d2 = gradient(d);
  else
    d = diff(motion_tseries,1,1);
    d = cat(1,zeros(1,6),d);
    d2 = diff(motion_tseries,2,1);
    d2 = cat(1,zeros(2,6),d2);
  end  
else
  d = motion_tseries;
  d2 = 0;
end;

% convert rotations in degrees to approximate displacements in mm
d(:,4:6) = d(:,4:6)*rot_fact;
deriv = d; 

% extract rotation estimates (mm/frame)
drot = d(:,4:6);

% extract translation estimates
if ~nodyflag
  dtrans = d(:,1:3);
else
  % exclude y translation
  dtrans = d(:,[1,3]);
  d = d(:,[1,3:6]);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate absolute value of translation
trans = calc_stat(dtrans,absflag);

% calculate overall rotation
rot = calc_stat(drot,absflag);

% calculate motion_fdmat (absolute displacements)
motion_fdmat = abs(d);

% calculate overall motion
motion = calc_stat(d,absflag);
motion_fd = motion;

% calculate overall acceleration
accel = calc_stat(d2,absflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalize motion_fd
if ~isempty(motion_fd_orig)
  norm_fact = sum(motion_fd_orig)/sum(motion_fd);
  motion_fd = motion_fd * norm_fact;
  motion_fdmat = motion_fdmat * norm_fact;
  motion = motion * norm_fact;
  accel = accel * norm_fact;
  trans = trans * norm_fact;
  rot = rot * norm_fact;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove first frame from vectors used to calculate stats
%   for backward consistency with prior version
if size(motion_tseries,1) > 1
  motion = motion(2:end);
  trans = trans(2:end);
  rot = rot(2:end);
  if size(motion_tseries,1) > 2
    accel = accel(3:end);
  else
    accel = accel(2:end);
  end
end

% calculate mean translation
motion_stats.mean_trans = mean(trans);

% calculate maximum translation
motion_stats.max_trans = max(trans);

% calculate mean rotation
motion_stats.mean_rot = mean(rot);

% calculate maximum rotation
motion_stats.max_rot = max(rot);

% calculate mean motion
motion_stats.mean_motion = mean(motion);

% calculate std of motion
motion_stats.std_motion = std(motion);

% calculate mode motion
motion_stats.mode_motion = mode(motion);

% calculate median motion
motion_stats.med_motion = median(motion);

% calculate minimum motion
motion_stats.min_motion = min(motion);

% calculate maximum motion
motion_stats.max_motion = max(motion);

% calculate mean acceleration
motion_stats.mean_accel = mean(accel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate time and percent of scan below each threshold
nthresh = length(thresholds);
nvols = length(motion_fd);
motion_stats.subthresh_nvols = zeros(nthresh,1);
for i=1:nthresh
  motion_stats.subthresh_nvols(i) = sum(motion_fd<thresholds(i));
end;
motion_stats.subthresh_perc = 100*motion_stats.subthresh_nvols/nvols;
motion_stats.thresholds = thresholds;
motion_stats.nvols = nvols;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stat = calc_stat(vec,absflag)
  if absflag
    stat = sum(abs(vec),2);
  else
    stat = sqrt(sum(vec.^2,2));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
