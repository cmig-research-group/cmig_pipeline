function fname_mask = mmil_synthstrip_brainmask(fname,varargin)
%function fname_mask = mmil_synthstrip_brainmask(fname,[options])
%
% Purpose:
%   create brain mask using FreeSurfer's synthstrip
%
% Required Input:
%   fname: full path of input file name (mgz or nii)
%
% Optional Parameters:
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
%   'nthreads': number of PyTorch CPU threads (mri_synthstrip -t)
%     fixed so the mask does not depend on the number of CPUs on the host:
%     a voxel at the border threshold can change with the thread count
%     {default = 8}
%   'verbose': [0|1} display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fname_mask: output file containing brain mask volume
%
% Created:  03/11/2025 by Don Hagler
% Prev Mod: 03/11/2025 by Don Hagler
% Prev Mod: 06/04/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Last Mod: 10/30/2026 by Don Hagler
%

% Based on mmil_fsl_brainmask
% Created:  02/20/16 by Don Hagler
% Prev Mod: 10/01/18 by Sean Hatton
% Last Mod: 01/31/22 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mask = [];
if ~mmil_check_nargs(nargin,1), return; end

parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'outdir',[],[],...
  'binarize_flag',false,[],...
  'cleanup_flag',true,[],...
  'nthreads',8,[1,Inf],...
  'verbose',false,[],...
  'forceflag',false,[],...
... % hidden variables for mri_convert
  'options','-f 0',[],...
  'convert_tags',{'options','verbose','forceflag'},[],...
});

if ~exist(fname,'file')
  error('file %s not found',fname);
end

if isempty(parms.outdir)
  if isempty(parms.fname_mask)
    parms.outdir = pwd;
  else
    parms.outdir = fileparts(parms.fname_mask);
    if isempty(parms.outdir), parms.outdir = pwd; end
  end
end

[fpath,fstem,fext] = cmig_get_fileparts(fname);

% set output file name
if isempty(parms.fname_mask)
  if parms.binarize_flag
    fname_mask = sprintf('%s/%s_brainmask_mask%s',parms.outdir,fstem,fext);
  else
    fname_mask = sprintf('%s/%s_brainmask%s',parms.outdir,fstem,fext);
  end
else
  fname_mask = parms.fname_mask;
end

% quit if output already exists
if exist(fname_mask,'file') && ~parms.forceflag
  return;
end

% parse string inputs for compiled deployment
if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

% extract first frame if input volume is multi-frame
%   fs_read_header only reads mgh/mgz, so pick the reader from fext;
%   both return volsz with four elements
switch lower(fext)
  case {'.mgh','.mgz'}
    [~,volsz] = fs_read_header(fname);
  case {'.nii','.nii.gz'}
    [~,~,volsz] = cmig_read_nifti(fname,1);
  otherwise
    error('unsupported file extension %s for %s',fext,fname);
end
if volsz(4)>1
  fprintf('%s: extracting first frame...\n',mfilename);
  tic;
  fname_in = sprintf('%s/%s_f0%s',parms.outdir,fstem,fext);
  args = mmil_parms2args(parms,parms.convert_tags);
  fs_mri_convert(fname,fname_in,args{:});
  toc;
else
  fname_in = fname;
end

% create brain mask
if parms.verbose
  fprintf('%s: creating brainmask...\n',mfilename);
  tic;
end
[s,~] = unix('which mri_synthstrip');
if s==0
  synthstrip_cmd = 'mri_synthstrip';
else
  synthstrip_cmd = 'synthstrip-singularity';
end
cmd = sprintf('%s -i %s -t %d',synthstrip_cmd,fname_in,parms.nthreads);
if parms.binarize_flag
  cmd = sprintf('%s -m %s',cmd,fname_mask);
else
  cmd = sprintf('%s -o %s',cmd,fname_mask);
end
[s,result,errmsg] = mmil_unix(cmd);
if s, error('cmd %s failed:\n%s%s',cmd,result,errmsg); end
if parms.verbose, toc; end

% check that output file exists
if ~exist(fname_mask,'file')
	error('failed to create brain mask... check input data');
end

if parms.cleanup_flag && ~strcmp(fname_in,fname), delete(fname_in); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

