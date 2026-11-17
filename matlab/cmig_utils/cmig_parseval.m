function val = cmig_parseval(val)
%function val = cmig_parseval(val)
%
% Purpose: to parse a string to return a
%   number, string, cell array, or numeric array
%
% Created:  12/25/2025 by Don Hagler
% Last Mod: 01/06/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on mmil_parsecells
%   Created:  04/01/09 by Jason Sherfey
%   Prev Mod: 01/30/18 by Don Hagler
%   Last Mod: 07/02/19 by Feng Xue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try if isnan(val), return; end; end
if isnumeric(val), return; end
if ischar(val) || isstring(val)
  % nan
  if strcmpi(val,'nan'), val = nan; return; end
  % numeric
  [tmp status] = isnum(val);
  if status, val = tmp; return; end
  % cell array
  [tmp status] = iscellarr(val);
  if status, val = tmp; return; end
  % numeric array
  [tmp status] = isnumarr(val);
  if status, val = tmp; return; end
  % do not convert to numeric if negative sign is not at beginning
  if ~isempty(regexp(val,'-','ONCE')) && isempty(regexp(strtrim(val),'^-','ONCE')), return; end
  % convert to numeric if all characters are numbers, spaces, brackets, periods, or negative signs
  if length(regexp(val,'[\d\s\[\].-]'))==length(val) &&...
     ~isempty(str2double(val)) &&...
     ~isnan(str2double(val))
    val = str2double(val);
    return;
  end      
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val status] = isnum(val)
  status = 0;
  try
    if any(regexp(val,'[\[\]]'))
      val = eval(val);
    end
  catch
    return;
  end
  if isnumeric(val), status = 1; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val status] = iscellarr(val)
  status = 0;
  if isempty(regexp(val,'^{.*}','match')), return; end
  val = regexp(val,'[^,{}'']+','match');
  if isstring(val), val = cellstr(val); end
  keep = {};
  for n = 1:length(val)
    % nan
    if strcmpi(val{n},'nan'), keep{end+1} = nan; continue; end
    % numeric
    [tmp status] = isnum(val{n});
    if status, keep{end+1} = tmp; continue; end
    % numeric array
    [tmp status] = isnumarr(val{n});
    if status==1, keep{end+1} = tmp; continue; end
    if status==2, keep(end+1:end+length(tmp)) = tmp; continue; end
    if isempty(regexp(val{n},'\w','ONCE')), continue; end
    keep(end+1) = regexp(val{n},'[^'']+','match');
  end
  val = keep;
  status = 2;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [val status] = isnumarr(val)
  status = 0;
  brk = regexp(val,'[\[^\]$]');
  len = length(brk);
  if len < 2 || mod(len,2)~=0, return; end
  if ~(brk(1)==1 && brk(end)==length(val)), return; end
  if len > 2
    tmp = val;
    val = {};
    for n = 1:2:len
      val{end+1} = str2double(tmp(brk(n):brk(n+1)));
    end
  else
    val = str2double(val);
  end
  if isnumeric(val), status = 1; end
  if iscell(val),    status = 2; end
  if isnan(val),     status = 0; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


