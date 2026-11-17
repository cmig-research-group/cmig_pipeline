function [fwhm,fwhm_frames] = cmig_3dFWHM(fname,varargin)
%function [fwhm,fwhm_frames] = cmig_3dFWHM(fname,[options])
%
% Purpose:
%   calculate spatial smoothness using AFNI's 3dFWHMx
%
% Usage:
%   [fwhm,fwhm_frames] = cmig_3dFWHM(fname,'key1',value1,...)
%
% Required Input:
%   fname: full path of input file name (mgz or nii)
%
% Optional Parameters:
%   'fname_mask': binary brain mask
%     if not supplied, will use 3dFWHMx -automask option
%     {default = []}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': stem for the output files, so two measurements of the same input
%     can share one directory (e.g. detrended and undetrended arms)
%     if empty, the input file stem is used
%     {default = []}
%   'detrend_flag': [0|1|2] how signal is removed before estimating smoothness
%     0: nothing removed, so the estimate reflects anatomical structure
%          plus noise
%     1: 3dFWHMx -detrend, which also applies -unif
%     2: 3dFWHMx -demed, subtracting each voxel's median across frames
%     AFNI ignores 1 and 2 for single frame input, which cannot be detrended
%     {default = 1}
%   'detrend_order': polynomial order used when detrend_flag = 1
%     if empty, AFNI chooses q = nframes/30, which is 0 for fewer than 30
%       frames, i.e. no detrending at all
%     {default = []}
%   'geom_flag': [0|1] combine the per-frame estimates with the geometric mean
%     0: arithmetic mean (3dFWHMx -arith)
%     1: geometric mean (3dFWHMx -geom)
%     AFNI excludes frames it could not estimate from either mean
%     {default = 0}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   fwhm: vector of FWHM for x, y, and z
%     an axis AFNI could not estimate is returned as NaN
%   fwhm_frames: matrix of per-frame FWHM values, size = [3,nframes]
%     frames AFNI could not estimate are returned as NaN
%
% Created:  08/17/26 by Don Hagler
% Prev Mod: 08/17/26 by Don Hagler
% Prev Mod: 10/09/2026 by Don Hagler
% Last Mod: 10/15/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize output
fwhm = nan(3,1);
fwhm_frames = [];

% check input
if ~mmil_check_nargs(nargin,1), return; end;
parms = check_input(fname,varargin);

% calculate smoothness for brain voxels
%   NOTE: the input conversion happens inside calculate_fwhm, in the branch
%   that decides whether 3dFWHMx runs, rather than here. Converting here would
%   convert an input for a call a cached re-run never makes
[fwhm,fwhm_frames] = calculate_fwhm(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
  ...
    'fname_mask',[],[],...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'detrend_flag',1,[0,1,2],...
    'detrend_order',[],[],...
    'geom_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'out_orient',[],[],...
    'ext','.nii',{'.nii','.nii.gz'},...
    'convert_tags',{'out_orient','verbose','forceflag'},[],...
  });
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  % NOTE: created here rather than as a side effect of convert_input, so the
  %   output directory does not depend on whether a conversion happens
  mmil_mkdir(parms.outdir);
  if ~isempty(parms.fname_mask) && ~exist(parms.fname_mask,'file')
    error('file %s not found',parms.fname_mask);
  end;
  if ~isempty(parms.detrend_order)
    if parms.detrend_order<0 || parms.detrend_order~=round(parms.detrend_order)
      error('detrend_order must be a non-negative integer');
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = convert_input(fname,parms)
  [fpath,fstem,fext] = cmig_get_fileparts(fname);
  switch fext
    case {'.mgh','.mgz'}
      fname_in = sprintf('%s/%s%s',parms.outdir,fstem,parms.ext);
      args = mmil_parms2args(parms,parms.convert_tags);
      fs_mri_convert(fname,fname_in,args{:});
    case {'.nii','.nii.gz'}
      fname_in = fname;
    otherwise
      error('invalid file extension: %s',fext);
  end;
  fname = fname_in;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = run_cmd(cmd)
  [errcode,r] = mmil_unix(cmd);
  if errcode
    fprintf('%s: WARNING: cmd %s failed:\n%s',mfilename,cmd,r);
    errcode = 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cmd = set_cmd(parms,fname_out,fname_frames)
  % mask
  if isempty(parms.fname_mask)
    mask_opt = '-automask';
  else
    mask_opt = sprintf('-mask %s',parms.fname_mask);
  end;
  % signal removal
  switch parms.detrend_flag
    case 0
      detrend_opt = '';
    case 1
      if isempty(parms.detrend_order)
        detrend_opt = '-detrend';
      else
        detrend_opt = sprintf('-detrend %d',parms.detrend_order);
      end;
    case 2
      detrend_opt = '-demed';
  end;
  % combination across frames
  if parms.geom_flag
    mean_opt = '-geom';
  else
    mean_opt = '-arith';
  end;
  % classic FWHM estimates
  % NOTE: 3dFWHMx from AFNI_17 on prints the ACF parameters unless given
  %   -ShowMeClassicFWHM, which older versions reject as unknown
  fwhmx_help = get_help;
  if ~isempty(strfind(fwhmx_help,'-ShowMeClassicFWHM'))
    classic_opt = '-ShowMeClassicFWHM';
  else
    classic_opt = '';
  end;
  % NOTE: newer 3dFWHMx (e.g. AFNI 26.2.06) always computes the ACF and
  %   writes it to 3dFWHMx.1D in the working directory unless given
  %   -acf NULL, which older versions reject as unknown
  if ~isempty(strfind(fwhmx_help,'-acf'))
    acf_opt = '-acf NULL';
  else
    acf_opt = '';
  end;
  % per-frame values go to fname_frames, the mean across frames to stdout
  cmd = sprintf('3dFWHMx -input %s %s %s %s %s %s -out %s > %s',...
    parms.fname,mask_opt,detrend_opt,mean_opt,classic_opt,acf_opt,...
    fname_frames,fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = get_help()
  % help text of the 3dFWHMx on PATH, empty if it fails
  [status,r] = mmil_unix('3dFWHMx -help');
  if status, r = ''; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = read_vals(fname)
  vals = [];
  fid = fopen(fname,'r');
  if fid<0, return; end;
  vals = fscanf(fid,'%f',[3,Inf]);
  fclose(fid);
  % AFNI writes a non-positive value for an estimate it could not make
  vals(vals<=0) = NaN;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fwhm,fwhm_frames] = calculate_fwhm(parms)
  fwhm = nan(3,1);
  fwhm_frames = [];
  % set output file names
  if isempty(parms.outstem)
    [fpath,fstem] = cmig_get_fileparts(parms.fname);
  else
    fstem = parms.outstem;
  end;
  outstem = sprintf('%s/%s',parms.outdir,fstem);
  fname_out = sprintf('%s_fwhm.txt',outstem);
  fname_frames = sprintf('%s_fwhm.1D',outstem);
  % check for existing output, treating an empty file as missing
  runflag = parms.forceflag;
  for fname_tmp = {fname_out,fname_frames}
    k = dir(fname_tmp{1});
    if isempty(k) || k.bytes==0, runflag = 1; end;
  end;
  if runflag
    if parms.verbose
      fprintf('%s: calculating FWHM...\n',mfilename);
    end;
    % convert to nii if necessary, only now that 3dFWHMx is going to run
    %   the stem above comes from the original name, which convert_input
    %   preserves, so no output name moves
    parms.fname = convert_input(parms.fname,parms);
    if ~isempty(parms.fname_mask)
      parms.fname_mask = convert_input(parms.fname_mask,parms);
    end;
    delete_output(fname_out,fname_frames,parms);
    errcode = run_cmd(set_cmd(parms,fname_out,fname_frames));
    if errcode, return; end;
  end;
  % mean across frames, calculated by AFNI
  vals = read_vals(fname_out);
  if numel(vals)<3
    fprintf('%s: WARNING: failed to read FWHM values from %s\n',...
      mfilename,fname_out);
    return;
  end;
  fwhm = vals(1:3,1);
  % per-frame values
  fwhm_frames = read_vals(fname_frames);
  if parms.verbose
    fprintf('%s: FWHM = %0.2f, %0.2f, %0.2f\n',...
      mfilename,fwhm(1),fwhm(2),fwhm(3));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delete_output(fname_out,fname_frames,parms)
  for fname_tmp = {fname_out,fname_frames}
    if exist(fname_tmp{1},'file')
      if parms.verbose
        fprintf('%s: removing old FWHM output %s...\n',mfilename,fname_tmp{1});
      end;
      delete(fname_tmp{1});
    end;
  end;
return;

