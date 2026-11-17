function mmil_write_csv(fname,data,varargin)
%function mmil_write_csv(fname,data,[options])
%
% Required Input:
%  fname: output file name
%  data: 2D data matrix of numbers or cell array of numbers and/or strings
%
% Optional Parameters:
%   'row_labels': cell array of row labels
%     {default: []}
%   'col_labels': cell array of column labels
%     {default: []}
%   'separator': separator between columns
%     {default: ','}
%   'firstcol_label': label for first column
%     {default: []}
%   'quote_flag': [0|1] put double quotes around column labels and strings
%     and a trailing space after each number
%     with 0, no quotes (except a string containing the separator)
%     and no trailing space
%     {default: 1}
%
% Created:  05/23/09 by Don Hagler
% Prev Mod: 05/12/20 by Don Hagler
% Last Mod: 11/11/26 by Don Hagler
%   

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'row_labels',[],[],...
  'col_labels',[],[],...
  'separator',',',[],...
  'firstcol_label',[],[],...
  'quote_flag',true,[false true],...
});

[nrows,ncols] = size(data);

if ~isempty(parms.row_labels) & length(parms.row_labels) ~= nrows
  error('number of row labels (%d) does not match first dim of data (%d)',...
    length(parms.row_labels),nrows);
end;
if ~isempty(parms.col_labels) & length(parms.col_labels) ~= ncols
  error('number of col labels (%d) does not match second dim of data (%d)',...
    length(parms.col_labels),ncols);
end;

% check if separator is \t, change to actual tab character
if strcmp(parms.separator,'\t'), parms.separator = sprintf(parms.separator); end;

data_is_cellarray = iscell(data);

fid = fopen(fname,'w');
if fid<0, error('failed to open %s for writing',fname); end;
if ~isempty(parms.col_labels)
  if ~isempty(parms.row_labels)
    fprintf(fid,'%s%s',parms.firstcol_label,parms.separator);
  end;
  for j=1:ncols
    if j>1, fprintf(fid,'%s',parms.separator); end
    write_string(fid,parms.col_labels{j},parms);
  end;
  fprintf(fid,'\n');
end;
for i = 1:nrows
  if ~isempty(parms.row_labels)
    fprintf(fid,'%s%s',parms.row_labels{i},parms.separator);
  end;
  for j = 1:ncols
    if j>1, fprintf(fid,'%s',parms.separator); end
    if data_is_cellarray, tvar = data{i,j};
    else  tvar = data(i,j);
    end
    if isempty(tvar)
    elseif ischar(tvar)
      write_string(fid,tvar,parms);
    elseif ~parms.quote_flag && isfloat(tvar) && isscalar(tvar)
      write_number(fid,tvar);
    else
      try
        mmil_write_value(fid,tvar);
      catch me
        error('element {%d,%d} of input data matrix is invalid:\n%s',...
          i,j,me.message);
      end;
    end;
  end
  fprintf(fid,'\n');
end
fclose(fid);


return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_string(fid,str,parms)
  if parms.quote_flag || contains(str,parms.separator)
    fprintf(fid,'"%s"',str);
  else
    fprintf(fid,'%s',str);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% same formats as mmil_write_value, without the trailing space
function write_number(fid,value)
  if isinf(value)
    if value<0
      fprintf(fid,'-Inf');
    else
      fprintf(fid,'Inf');
    end;
  elseif mmil_isint(value)
    fprintf(fid,'%d',value);
  else
    fprintf(fid,'%0.8f',value);
  end;
return;

