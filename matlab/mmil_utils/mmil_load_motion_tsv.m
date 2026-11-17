function motion_data = mmil_load_motion_tsv(fname_motion,varargin)
%function motion_data = mmil_load_motion_tsv(fname_motion,varargin)
%
% Purpose: load _desc-motion_timeseries.tsv file from BIDS derived data
%
% Usage:
%  motion_data = mmil_load_motion_tsv(fname_motion,'key1', value1,...); 
%
% Required Parameters:
%   fname_motion: text file containing 6 columns of motion estimates
%     in tab-delimited format
%
% Optional parameters:
%  'skipTRs': number of initial repetitions to remove
%    {default: 0}
%  'nframes': number of frames in corresponding data file
%    If supplied, will check that number of time points in fname_motion matches
%    {default: []}
%
% Output:
%   motion_data: matrix of motion parameter estimates with size [nframes,npar]
%     order of parameters: dx, dy, dz, rx, ry, rz
%     with units of mm and degrees
%
% Created:  09/04/25 Don Hagler
% Prev Mod: 09/04/25 Don Hagler
% Last Mod: 11/11/26 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on mmil_load_motion_1D
%   Created:  03/06/12 Don Hagler
%   Last Mod: 09/28/18 Don Hagler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

motion_data = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'skipTRs',0,[0 Inf],...
  'nframes',[],[],...
...
  'reorder_vec',[5,6,4,2,3,1],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  load fname_motion
if ~exist(fname_motion,'file')
  error('file %s not found',fname_motion);
else
  % read tab-delimited text file
  motion_data = mmil_readtext(fname_motion,'\t');
  % discard first row and first column
  motion_data = cell2mat(motion_data(2:end,2:end));
  
  if size(motion_data,2)>6
    motion_data = motion_data(:,2:7);
  end;
  if isempty(parms.nframes)
    parms.nframes = size(motion_data,1);
  end;
  if size(motion_data,1) ~= parms.nframes
    error('number of time points (%d) in fname_motion %s does not match number of frames (%d)',...
      size(motion_data,1),fname_motion,parms.nframes);
  end;
  if parms.skipTRs>0 && parms.skipTRs<parms.nframes
    motion_data = motion_data(parms.skipTRs+1:end,:);
  end;
end;

% reorder motion parameter estimates
%   from rz,rx,ry,dz,dx,dy to dx,dy,dz,rx,ry,rz
motion_data = motion_data(:,parms.reorder_vec);

% change rotations to be clockwise
motion_data(:,4:6) = -motion_data(:,4:6);


