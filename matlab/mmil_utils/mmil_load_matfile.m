function output = mmil_load_matfile(fname,required_vars)
%function output = mmil_load_matfile(fname,required_vars)
%
% Purpose: load mat file, delete if corrupt
%
% Required input:
%   fname: full path of mat file
%
% Optional input:
%   required_vars: cell array of required variable names
%     {default = []}
%
% Output:
%   output: struct with contents of mat file as fields
%
% Created:  05/14/2020 by Don Hagler
% Prev Mod: 03/30/2021 by Don Hagler
% Last Mod: 04/20/2021 by Don Hagler
%

output = [];
if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('required_vars','var'), required_vars = []; end;

if ~isempty(required_vars) && ~iscell(required_vars)
  required_vars = {required_vars};
end;

errmsg = [];
rmflag = 0;

d = dir(fname);
if isempty(d)
  errmsg = sprintf('file %s not found',fname);
elseif d.bytes==0
  % delete if mat file has file size zero
  errmsg = sprintf('file with size 0 found:\n%s\ndeleting empty file... rerun to recreate and proceed',fname);
  rmflag = 1;
else
  % try to load mat file and delete if corrupt
  try
    output = load(fname);
  catch me
    if ~isempty(regexp(me.message,'corrupt'))
      errmsg = sprintf('corrupt file found:\n%s\ndeleting corrupt file... rerun to recreate and proceed',fname);
      rmflag = 1;
    elseif ~isempty(regexp(me.message,'Cannot read file'))
      errmsg = sprintf('unreadable file found:\n%s\ndeleting unreadable file... rerun to recreate and proceed',fname);
      rmflag = 1;
    elseif ~isempty(regexp(lastwarn,'Unexpected end-of-file while reading compressed data'))
      errmsg = sprintf('truncated file found:\n%s\ndeleting truncated file... rerun to recreate and proceed',fname);
      rmflag = 1;
    else
      errmsg = sprintf('failed to load mat file:\n%s',me.message);
    end;
  end;
end;

% check for required variables in mat file
%   delete if not found
if isempty(errmsg)
  for i=1:length(required_vars)
    varname = required_vars{i};
    if ~isfield(output,varname)
      %% todo: error?
      errmsg = sprintf('mat file missing %s found:\n%s\ndeleting... rerun to recreate and proceed',...
        varname,fname);
%      keyboard
      rmflag = 1;
    end;
  end;
end;

if rmflag, delete(fname); end;
if ~isempty(errmsg), error(errmsg); end;

