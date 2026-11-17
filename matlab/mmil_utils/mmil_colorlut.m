function [roicodes, roinames, rgbv] = mmil_colorlut(varargin)
%function [roicodes, roinames, rgbv] = mmil_colorlut([options])
%
% Purpose: reads a freesurfer color lookup table
%   and combines with mmil color lookup table
%
% Optional input ('key',value pairs):
%   'FS_version': FreeSurfer version
%     if specified, used to determine $FREESURFER_HOME
%     {default = []}
%   'fname_fscolorlut': name of FreeSurfer color LUT file
%     {default = $FREESURFER_HOME/FreeSurferColorLUT.txt}
%   'fname_excolorlut': name of MMIL color LUT file
%     {default = $CMIG_PIPELINE_DIR/parms/CMIG_ColorLUT.txt}
%
% Output:
%  roicodes: vector of numeric ROI codes
%  roinames: cell array of ROI names
%  rgbv: matrix of rgb values for each ROI
%
% Created:  12/05/20 by Don Hagler
% Prev Mod: 12/07/20 by Don Hagler
% Last Mod: 10/26/2026 by Don Hagler
%

roicodes = [];
roinames = [];
rgbv = [];

parms = mmil_args2parms(varargin,{...
  'FS_version',[],[],...
  'fname_fscolorlut',[],[],...
  'fname_excolorlut',[],[],...
});

fshome = mmil_get_fshome(parms.FS_version);
if isempty(fshome) || ~exist(fshome,'dir')
  error('FREESURFER_HOME not found');
end;

if isempty(parms.fname_fscolorlut)
  parms.fname_fscolorlut = [fshome '/FreeSurferColorLUT.txt'];
end;
if ~exist(parms.fname_fscolorlut,'file'), error('file %s not found',parms.fname_fscolorlut); end;

if isempty(parms.fname_excolorlut)
  parms.fname_excolorlut = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/CMIG_ColorLUT.txt'];
end;
if ~exist(parms.fname_excolorlut,'file'), error('file %s not found',parms.fname_excolorlut); end;

[fs_roicodes,fs_roinames,fs_rgbv] = fs_colorlut(parms.fname_fscolorlut);
[ex_roicodes,ex_roinames,mmil_rgbv] = fs_colorlut(parms.fname_excolorlut);

roicodes = cat(1,fs_roicodes,ex_roicodes);
roinames = cat(1,fs_roinames,ex_roinames);
rgbv = cat(1,fs_rgbv,mmil_rgbv);

return;
    
    