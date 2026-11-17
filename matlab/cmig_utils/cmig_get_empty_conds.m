function empty_conds = cmig_get_empty_conds(fname_condnames,run_tag)
%function empty_conds = cmig_get_empty_conds(fname_condnames,[run_tag])
%
% Purpose: report which deconvolution conditions had no events, so that a
%   caller writes them as missing rather than as the zero an all-zero
%   regressor returns
%
% Required Input:
%   fname_condnames: full path of the deconv_condnames mat file
%
% Optional Input:
%   run_tag: tag identifying one run, in either the 'run{NN}' or 'r{NN}' form,
%     to ask about that run alone; if empty, asks whether the condition is
%     absent from every run
%     {default = []}
%
% Output:
%   empty_conds: [1,ncondnames] logical, true where the condition has no
%     events; empty when the file or the variables are not there, as for a
%     file written before they were recorded
%
% Created:  09/16/2026 by Don Hagler
% Prev Mod: 09/16/2026 by Don Hagler
% Last Mod: 09/17/2026 by Don Hagler
%

empty_conds = [];
if ~mmil_check_nargs(nargin,1), return; end
if ~exist('run_tag','var'), run_tag = []; end

if isempty(fname_condnames) || ~exist(fname_condnames,'file'), return; end
vars = {whos('-file',fname_condnames).name};

if ~isempty(run_tag) && all(ismember({'empty_conds_run','run_numbers'},vars))
  tmp = load(fname_condnames,'empty_conds_run','run_numbers');
  % the tag is 'run01' from the analyze module and 'r01' in the column names
  rnum = str2double(regexp(run_tag,'\d+','match','once'));
  r = find(tmp.run_numbers == rnum,1);
  if ~isempty(r)
    empty_conds = logical(tmp.empty_conds_run(r,:));
  end
  return;
end

if ismember('empty_conds',vars)
  tmp = load(fname_condnames,'empty_conds');
  empty_conds = logical(tmp.empty_conds(:)');
end

return;
