function epi_estimateB0(fname_for,fname_rev,varargin)
%function epi_estimateB0(fname_for,fname_rev,[options])
%
% Purpse: estimate distortions caused by B0 susceptibility artifacts
%
% Usage: epi_estimateB0(fname_for,fname_rev,'key','value'...)
%
% Required Parameters:
%   fname_for: input file name with "forward" phase encode polarity
%   fname_rev: input file name with "reverse" phase encode polarity
%
% Optional Parameters:
%   'fname_dx': output file name for displacement field
%     If empty, will construct from fname_for and fname_rev
%     {default = []}
%   'fname_for_out': output file name for forward scan after B0 unwarping
%     If empty, will not save this output
%     {default = []}
%   'fname_rev_out': output file name for reverse scan after B0 unwarping
%     If empty, will not save this output
%     {default = []}
%   'epic_flag': [0|1] whether to use epic binary instead of unwarpB0
%     {default = 0}
%   'epic_inorm_flag': [0|1] whether to rescale image intensities to match
%     Specific to epic, so ignored if epic_flag = 0
%     {default = 0}
%   'fname_for_noJ_out': output file name for forward scan after B0 unwarping
%     without Jacobian intensity correction (only applies if epic_flag = 1)
%     If empty, will set based on fname_for_out
%     {default = []}
%   'fname_rev_noJ_out': output file name for reverse scan after B0 unwarping
%     without Jacobian intensity correction (only applies if epic_flag = 1)
%     If empty, will set based on fname_for_out
%     {default = []}
%   'swapxy_flag': [0|1] whether to swap x and y dimensions before estimating B0
%     {default = 0}
%   'kernelWidthMax': maximum smoothing kernel width (must be odd integer)
%     {default = 25}
%   'lambda2': regularization factor (determines tightness of morph) 
%     {default = 1100}
%   'forceflag': [0|1] whether to overwrite existing fname_dx
%     {default = 0}
%
% Created:  02/02/07 by Don Hagler
% Prev Mod: 06/30/20 by Don Hagler
% Prev Mod: 07/07/20 by Don Hagler
% Last Mod: 07/14/20 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'fname_dx',[],[],...
  'fname_for_out',[],[],...
  'fname_rev_out',[],[],...
  'epic_flag',false,[false true],...
  'epic_inorm_flag',false,[false true],...
  'fname_for_noJ_out',[],[],...
  'fname_rev_noJ_out',[],[],...
  'swapxy_flag',false,[false true],...
  'kernelWidthMax',25,[1:100],...
  'lambda2',1100,[1:10000],...
  'forceflag',false,[false true],...
...
  'cleanupflag',true,[false true],...
  'outext','.mgz',{'.mgz','.mgh'},...
});

% check that kernelWidthMax is odd integer
if ~mmil_isint(parms.kernelWidthMax)
  fprintf('%s: WARNING: rounding kernelWidthMax (%g) to an integer...',...
    mfilename,parms.kernelWidthMax);
  parms.kernelWidthMax = round(parms.kernelWidthMax);
end;
if ~mod(parms.kernelWidthMax,2)
  fprintf('%s: WARNING: adding 1 to make kernelWidthMax (%d) odd...\n',...
    mfilename,parms.kernelWidthMax);
  parms.kernelWidthMax = parms.kernelWidthMax + 1;
end;

% check that input volumes exist
if ~exist(fname_for,'file')
  error('file %s not found',fname_for);
end;
if ~exist(fname_rev,'file')
  error('file %s not found',fname_rev);
end;

% create fname_B0dx from fname_for and fname_ref stems
if isempty(parms.fname_dx)
  [tpath,tstem_rev,text] = fileparts(fname_rev);
  [tpath,tstem_for,text] = fileparts(fname_for);
  parms.fname_dx = [tpath '/' tstem_for 'VS' tstem_rev '_B0dx' parms.outext];
end;
if exist(parms.fname_dx,'file') && ~parms.forceflag, return; end;

% create tmp directory
[tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_dx);
if isempty(tmp_path), tmp_path = pwd; end;
tmpdir = [tmp_path '/tmp_estimateB0'];
mmil_mkdir(tmpdir);

% check that input volumes share dimensions and orientations
[dimsmatch,orientmatch,voxmatch] = ...
  fs_vol_match(fname_for,fname_rev,tmpdir,parms.forceflag);
if ~dimsmatch
  error('mismatch in input volume dimensions for %s and %s',...
    fname_for,fname_rev);
end;
if ~orientmatch
  error('mismatch in volume orientation for %s and %s',...
    fname_for,fname_rev);
end;
if ~voxmatch
  error('mismatch in voxel to RAS transformation for %s and %s',...
    fname_for,fname_rev);
end;

% set fname_for_noJ_out
if isempty(parms.fname_for_noJ_out) && ~isempty(parms.fname_for_out)
  [fdir,fstem,fext] = fileparts(parms.fname_for_out);
  parms.fname_for_noJ_out = sprintf('%s/%s_noJ%s',fdir,fstem,fext);
end;
% set fname_rev_noJ_out
if isempty(parms.fname_rev_noJ_out) && ~isempty(parms.fname_rev_out)
  [fdir,fstem,fext] = fileparts(parms.fname_rev_out);
  parms.fname_rev_noJ_out = sprintf('%s/%s_noJ%s',fdir,fstem,fext);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if swapxy_flag, swap x and y
if parms.swapxy_flag
  orig_orient = fs_read_orient(fname_for);
  pref_orient = orig_orient([2,1,3]);
  fprintf('%s: reorienting input volumes to %s...\n',mfilename,pref_orient);
  tmp_fname_for = [tmpdir '/for.mgz'];
  tmp_fname_rev = [tmpdir '/rev.mgz'];
  epi_reorient_vol(fname_for,tmp_fname_for,pref_orient);
  epi_reorient_vol(fname_rev,tmp_fname_rev,pref_orient);
else
  tmp_fname_for = fname_for;
  tmp_fname_rev = fname_rev;
end;

% set tmp output file names
tmp_fname_dx = 'dx.mgz';
tmp_fname_for_out = 'for_B0uw.mgz';
tmp_fname_rev_out = 'rev_B0uw.mgz';
if parms.epic_flag
  tmp_fname_for_noJ_out = 'fB0uwNoJ.mgz';
  tmp_fname_rev_noJ_out = 'rB0uwNoJ.mgz';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate command string
if parms.epic_flag
  cmd = 'epic';
  cmd = sprintf('%s -averageOverNbrs -resample -xs 5 -ys 5 -zs 5 -smoothDfieldKernelWidth 5\\\n',cmd);
  if parms.epic_inorm_flag
    cmd = sprintf('%s -rescaleForwardImage\\\n',cmd);
  end;
else
  cmd = 'unwarpB0';
end;
cmd = sprintf('%s -f %s -r %s\\\n',cmd,tmp_fname_for,tmp_fname_rev);
cmd = sprintf('%s -fo %s -ro %s\\\n',cmd,tmp_fname_for_out,tmp_fname_rev_out);
cmd = sprintf('%s -do %s -od %s\\\n',cmd,tmp_fname_dx,tmpdir);
cmd = sprintf('%s -kernelWidthMax %g -lambda2 %g -lambda2P %g\\\n',...
              cmd,parms.kernelWidthMax,parms.lambda2,parms.lambda2);

% run command
fprintf('%s: cmd = %s\n',mfilename,cmd);
[status,result] = mmil_unix(cmd);
if status
  error('cmd %s failed:\n%s\n',cmd,result);
else
  disp(result);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% copy tmp output files or reorient as necessary
tmp_fname_dx = [tmpdir '/' tmp_fname_dx];
tmp_fname_for_out = [tmpdir '/' tmp_fname_for_out];
tmp_fname_rev_out = [tmpdir '/' tmp_fname_rev_out];
if parms.epic_flag
  tmp_fname_for_noJ_out = [tmpdir '/' tmp_fname_for_noJ_out];
  tmp_fname_rev_noJ_out = [tmpdir '/' tmp_fname_rev_noJ_out];
end;
flist1 = {tmp_fname_dx,tmp_fname_for_out,tmp_fname_rev_out};
flist2 = {parms.fname_dx,parms.fname_for_out,parms.fname_rev_out};
if parms.epic_flag
  flist1 = cat(2,flist1,{tmp_fname_for_noJ_out,tmp_fname_rev_noJ_out});
  flist2 = cat(2,flist2,{parms.fname_for_noJ_out,parms.fname_rev_noJ_out});
end;
if parms.swapxy_flag
  fprintf('%s: reorienting output volumes to %s...\n',mfilename,orig_orient);
  for f=1:length(flist1)
    if isempty(flist2{f}), continue; end;
    epi_reorient_vol(flist1{f},flist2{f},orig_orient);
  end;
else
  for f=1:length(flist1)
    if isempty(flist2{f}), continue; end;
    [tmp,tmp,ext1] = fileparts(flist1{f});
    [tmp,tmp,ext2] = fileparts(flist2{f});
    if strcmp(ext1,ext2)
      cmd = sprintf('mv %s %s',flist1{f},flist2{f});
      [status,result,errmsg] = mmil_unix(cmd);
      if status
        error('cmd %s failed:\n%s\n%s',cmd,result,errmsg);
      end;
    else
      fs_copy_mgh(flist1{f},flist2{f});
    end;
  end;
end;

% remove tmp dir
if parms.cleanupflag
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result,errmsg] = mmil_unix(cmd);
  if status
    fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n%s\n%s\n',mfilename,cmd,result,errmsg);
  end;
end;
