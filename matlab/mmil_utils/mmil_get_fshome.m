function fshome = mmil_get_fshome(FS_version)
%function fshome = mmil_get_fshome(FS_version)
%
% Created:  12/07/20 by Don Hagler
% Last Mod: 12/07/20 by Don Hagler
%

fshome = [];
if ~exist('FS_version','var'), FS_version = []; end;

if ~isempty(FS_version)
  cmd = sprintf('source $PUBSH/bin/SetUpFreeSurfer.csh %d; echo $FREESURFER_HOME\n',FS_version);
  [status,fshome,errmsg] = mmil_unix(cmd);
  if status
    error('failed to get FREESURFER_HOME for FS_version %d:\n%s',...
      FS_version,errmsg);
  end;
else
  fshome = getenv('FREESURFER_HOME');
end;

