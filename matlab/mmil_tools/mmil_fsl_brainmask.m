function fname_mask = mmil_fsl_brainmask(fname,varargin)
%function fname_mask = mmil_fsl_brainmask(fname,[options])
%
% Purpose:
%   create brain mask using FSL's brain extraction tool (bet)
%
% Required Input:
%   fname: full path of input file name (mgz or nii)
%
% Optional Parmaeters:
%   'fname_mask': output file name
%     if not supplied, will generate a name based on file stem of fname and outdir
%     {default = []}
%   'outdir': output directory
%     if not supplied, will use path of fname_mask
%     will use pwd if both outdir and fname_mask not supplied
%     {default = []}
%   'binarize_flag': [0|1] output binary mask
%     otherwise return masked brain image
%     {default = 0}
%   'cleanup_flag': [0|1] remove intermediate files
%     {default = 1}
%   'verbose': [0|1} display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_mask: output file containing brain mask volume
%
% Created:  02/20/16 by Don Hagler
% Prev Mod: 10/01/18 by Sean Hatton
% Prev Mod: 07/07/20 by Don Hagler
% Prev Mod: 02/24/21 by Don Hagler
% Prev Mod: 01/31/22 by Don Hagler
% Last Mod: 06/25/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mask = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'outdir',[],[],...
  'binarize_flag',false,[false true],...
  'cleanup_flag',true,[false true],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
... % hidden
  'out_orient',[],[],...
  'ext','.nii',{'.nii','.nii.gz'},...
  'convert_tags',{'out_orient','verbose','forceflag'},[],...
});

if ~exist(fname,'file')
  error('file %s not found',fname);
end;


if isempty(parms.outdir)
  if isempty(parms.fname_mask)
    parms.outdir = pwd;
  else
    parms.outdir = fileparts(parms.fname_mask);
    if isempty(parms.outdir), parms.outdir = pwd; end;
  end;
end;

[fpath,fstem,fext] = cmig_get_fileparts(fname);
suffix_brain = '_brainmask';
if parms.binarize_flag
  suffix = '_brainmask_mask';
else
  suffix = suffix_brain;
end;

% set output file name
if isempty(parms.fname_mask)
  fname_mask = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,fext);
else
  fname_mask = parms.fname_mask;
end;

% quit if output already exists
if exist(fname_mask,'file') && ~parms.forceflag
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

% convert to nii if necessary
switch fext
  case {'.mgh','.mgz'}
    fname_in = sprintf('%s/%s%s',parms.outdir,fstem,parms.ext);
    args = mmil_parms2args(parms,parms.convert_tags);
    fs_mri_convert(fname,fname_in,args{:});
  case {'.nii','.nii.gz'}
    fname_in = fname;
  otherwise
    error('invalid file extension: %s\n',mfilename,fext);
end;

% create brain mask
fname_tmp = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
fname_tmp_out = [];
if ~exist(fname_tmp,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: creating brainmask...\n',mfilename);
    tic;
  end;
  if strcmp(parms.ext,'.nii')
    cmd = sprintf('setenv FSLOUTPUTTYPE NIFTI; ');
  elseif strcmp(parms.ext,'.nii.gz')
    cmd = sprintf('setenv FSLOUTPUTTYPE NIFTI_GZ; ');
  else
    cmd = [];
  end;
  fname_tmp_out = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix_brain,parms.ext);
  cmd = sprintf('%sbet %s %s -R',cmd,fname_in,fname_tmp_out);
  if parms.binarize_flag
    cmd = sprintf('%s -m',cmd);
  end;
  [result,errmsg] = run_cmd(cmd);
  if parms.verbose, toc; end;

	% check that output file exists
	if ~exist(fname_tmp,'file')
		error('failed to create brain mask... check input data');
	end
end;

[~,~,fext] = cmig_get_fileparts(fname_mask);

% convert back to mgz
if ismember(fext,{'.mgh','.mgz'}) || ~isempty(parms.fname_mask)
  args = mmil_parms2args(parms,parms.convert_tags);
  fs_mri_convert(fname_tmp,fname_mask,args{:});
else
  fname_mask = fname_tmp;
end;

if parms.cleanup_flag
  if ~isempty(fname_tmp_out) && ~strcmp(fname_mask,fname_tmp_out) && exist(fname_tmp_out,'file'), delete(fname_tmp_out); end;
  if ~strcmp(fname_mask,fname_tmp) && exist(fname_tmp,'file'), delete(fname_tmp); end;
  if ~strcmp(fname_in,fname), delete(fname_in); end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result,errmsg] = run_cmd(cmd)
  [s,result,errmsg] = mmil_unix(cmd);
  if s, error('cmd %s failed:\n%s%s',cmd,result,errmsg); end;
return;

