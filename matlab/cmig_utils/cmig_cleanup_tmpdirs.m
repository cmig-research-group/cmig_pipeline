function errmsg = cmig_cleanup_tmpdirs(outdir, tmpdir)
%function errmsg = cmig_cleanup_tmpdirs(outdir, [tmpdir])
%
% Required Input:
%   outdir: directory containing temporary subdirectories to remove
%
% Optional Input:
%   tmpdir: specific temporary directory to remove; if provided, only this
%     directory is removed rather than all subdirectories of outdir
%
% Output:
%   errmsg: error message if removal failed; empty otherwise
%
% Created:  05/19/2026 by Don Hagler
% Last Mod: 05/29/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  errmsg = [];
  if nargin < 2 || isempty(tmpdir)
    % delete all subdirectories of outdir
    dlist = dir(sprintf('%s/*',outdir)); if isempty(dlist), return; end
    dlist = dlist(find([dlist.isdir] == 1)); if isempty(dlist), return; end
    dlist = dlist(find(~ismember({dlist.name},{'.','..'}))); if isempty(dlist), return; end
    tmpdirs = strcat(outdir,'/',{dlist.name});
  else
    % delete only the specified tmpdir
    if ~exist(tmpdir,'dir'), return; end
    tmpdirs = {tmpdir};
  end
  for t=1:length(tmpdirs)
    fprintf('%s: removing temporary directory %s...\n',mfilename,tmpdirs{t});
    try
      rmdir(tmpdirs{t},'s');
    catch
      errmsg = sprintf('failed to remove temporary directory %s',tmpdirs{t});
      fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    end
  end
return;
