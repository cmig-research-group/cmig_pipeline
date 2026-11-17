function [fname_snr,fname_mean,fname_std,fname_motion] = mmil_calc_tSNR(fname,varargin)
%function [fname_snr,fname_mean,fname_std,fname_motion] = mmil_calc_tSNR(fname,[options])
%
% Purpose:
% 	calculate temporal SNR
%
% Required Input:
% 	fname: full path of input file name (mgz or nii)
%
% Optional Parmaeters:
% 	'outdir': output directory
% 		{default = pwd}
%   'mcflag': [0|1] apply motion correction
%     {default = 1}
% 	'hp_sigma': high pass filter sigma (seconds)
% 		{default = 100}
% 	'TR': repetition time (seconds)
% 		{default = 1}
%   'skipTRs': number of TRs at beginning of scan to be ignored
%      for calculation of mean, stdev, and tSNR
%     {default = 8}
%   'method': method used for head motion estimation and correction
%     allowed values: 'fslmean', 'fsl', 'afni'
%       fslmean: FSL's mcflirt with -meanvol
%       fsl: FSL's mcflirt with -r fname_ref
%       afni: AFNI's 3dvolreg
%    {default = 'fslmean'}
%   'snr_suffix': suffix appended to the input file stem for the tSNR volume
%     {default = '_snr'}
% 	'verbose': [0|1] display status messages
% 		{default = 0}
% 	'forceflag': [0|1] overwrite existing output
% 		{default = 0}
%
% Output:
% 	fname_snr: output file containing tSNR volume
% 	fname_mean: output file containing mean volume
% 	fname_std: output file containing standard deviation volume
% 	fname_motion: output file containing motion estimates
%
% Created:  02/20/16 by Don Hagler
% Prev Mod: 07/07/20 by Don Hagler
% Prev Mod: 02/24/21 by Don Hagler
% Prev Mod: 05/30/23 by Don Hagler
% Prev Mod: 05/31/23 by Don Hagler
% Prev Mod: 03/19/25 by Don Hagler
% Prev Mod: 10/13/2026 by Don Hagler
% Prev Mod: 10/20/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_snr = []; fname_mean = []; fname_std = []; fname_motion = [];
if ~mmil_check_nargs(nargin,1), return; end

parms = mmil_args2parms(varargin,{...
	'outdir',[],[],...
  'mcflag',true,[false true],...
	'hp_sigma',100,[10,1000],...
	'TR',1,[0.01,100],...
  'skipTRs',8,[0,Inf],...
  'method','fslmean',{'fslmean','fsl','afni'},...
  'snr_suffix','_snr',[],...
	'verbose',false,[false true],...
	'forceflag',false,[false true],...
... % hidden
	'ext','.nii',{'.nii','.nii.gz'},...
	'convert_tags',{'verbose','forceflag'},[],...
});

if ~exist(fname,'file')
	error('file %s not found',fname);
end

if strcmp(parms.ext,'.nii')
  parms.precmd = sprintf('setenv FSLOUTPUTTYPE NIFTI; ');
elseif strcmp(parms.ext,'.nii.gz')
  parms.precmd = sprintf('setenv FSLOUTPUTTYPE NIFTI_GZ; ');
else
  parms.precmd = [];
end

% calculate highpass sigma in volumes
parms.hp_sigma = parms.hp_sigma / parms.TR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

[fpath,fstem,fext] = fileparts(fname);
if parms.skipTRs
  % remove initial skipTRs
  suffix = sprintf('_skip%dTRs',parms.skipTRs);
  outstem = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix);
  fname_tmp = sprintf('%s%s',outstem,parms.ext);
  if ~exist(fname_tmp,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: removing initial %d frames...\n',mfilename,parms.skipTRs);
      tic;
    end
    args = mmil_parms2args(parms,parms.convert_tags);
    fs_mri_convert(fname,fname_tmp,...
      'options',sprintf('--nskip %d',parms.skipTRs),args{:});
    if parms.verbose, toc; end
  end
  fname_in = fname_tmp;
else
  % convert to nii if necessary
  suffix = [];
  switch fext
    case {'.mgh','.mgz'}
      fname_in = sprintf('%s/%s%s',parms.outdir,fstem,parms.ext);
      args = mmil_parms2args(parms,parms.convert_tags);
      fs_mri_convert(fname,fname_in,args{:});
    case {'.nii','.nii.gz'}
      fname_in = fname;
    otherwise
      error('invalid file extension: %s\n',mfilename,fext);
  end
end

if ~strcmp(parms.method,'fslmean')
  % calculate voxel-wise mean
  suffix_ref = '_mean_ref';
  fname_ref = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix_ref,parms.ext);
  if ~exist(fname_ref,'file') || parms.forceflag
	  if parms.verbose
		  fprintf('%s: calculating voxel-wise mean...\n',mfilename);
		  tic;
	  end
    cmd = sprintf('%sfslmaths %s -Tmean %s',...
							    parms.precmd,fname_in,fname_ref);
	  run_cmd(cmd);
	  if parms.verbose, toc; end
  end
end

% motion correction
if parms.mcflag
  suffix = [suffix '_mcf'];
  outstem = sprintf('%s/%s%s',parms.outdir,fstem,suffix);
  fname_tmp = sprintf('%s%s',outstem,parms.ext);
  if strcmp(parms.method,'afni')
    fname_motion = [outstem '_motion.1D'];
  else
    fname_motion = [outstem '.par'];
  end
  if ~exist(fname_tmp,'file') || parms.forceflag
	  if parms.verbose
		  fprintf('%s: performing motion correction...\n',mfilename);
		  tic;
	  end
    if strcmp(parms.method,'afni')
      % AFNI_DONT_LOGFILE keeps 3dvolreg from appending to ~/.afni.log
      cmd = sprintf('env AFNI_DONT_LOGFILE=YES 3dvolreg -prefix %s',fname_tmp);
      cmd = sprintf('%s -base ''%s''',cmd,fname_ref);
      cmd = sprintf('%s -dfile %s',cmd,fname_motion);
      cmd = sprintf('%s %s',cmd,fname_in);
    elseif strcmp(parms.method,'fsl')
      cmd = sprintf('%smcflirt -in %s -out %s -r %s -plots',...
            parms.precmd,fname_in,outstem,fname_ref);
    else
      cmd = sprintf('%smcflirt -in %s -out %s -meanvol -plots',...
					  parms.precmd,fname_in,outstem);
    end
	  run_cmd(cmd);
	  if parms.verbose, toc; end
  end
  fname_in = fname_tmp;
end

% high pass temporal filtering
%   fslmaths -bptf keeps the temporal mean in FSL 5 and removes it in FSL 6,
%   so the mean is subtracted before the filter and added back after: the
%   result is the same under either version
suffix = [suffix '_hpf'];
fname_tmp = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_tmp,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: applying high pass temporal filter...\n',mfilename);
		tic;
	end
  fname_tmean = sprintf('%s/%s%s_tmean%s',parms.outdir,fstem,suffix,parms.ext);
  cmd = sprintf('%sfslmaths %s -Tmean %s',parms.precmd,fname_in,fname_tmean);
  run_cmd(cmd);
  cmd = sprintf('%sfslmaths %s -sub %s -bptf %0.1f -1 -add %s %s',...
                parms.precmd,fname_in,fname_tmean,parms.hp_sigma,fname_tmean,fname_tmp);
	run_cmd(cmd);
  delete(fname_tmean);
	if parms.verbose, toc; end
end
fname_in = fname_tmp;

% calculate voxel-wise mean
suffix = '_mean';
fname_mean = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_mean,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: calculating voxel-wise mean...\n',mfilename);
		tic;
	end
  cmd = sprintf('%sfslmaths %s -Tmean %s',...
							  parms.precmd,fname_in,fname_mean);
	run_cmd(cmd);
	if parms.verbose, toc; end
end

% calculate voxel-wise standard deviation
suffix = '_std';
fname_std = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_std,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: calculating voxel-wise standard deviation...\n',mfilename);
		tic;
	end
  cmd = sprintf('%sfslmaths %s -Tstd %s',...
								parms.precmd,fname_in,fname_std);
	run_cmd(cmd);
	if parms.verbose, toc; end
end

% calculate temporal SNR
suffix = parms.snr_suffix;
fname_tmp = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
if ~exist(fname_tmp,'file') || parms.forceflag
	if parms.verbose
		fprintf('%s: calculating temporal SNR...\n',mfilename);
		tic;
	end
  cmd = sprintf('%sfslmaths %s -div %s %s',...
							  parms.precmd,fname_mean,fname_std,fname_tmp);
	run_cmd(cmd);
	if parms.verbose, toc; end
end

% convert back to mgz
if ismember(fext,{'.mgh','.mgz'})
	args = mmil_parms2args(parms,parms.convert_tags);
  suffix_list = {'_mean','_std',parms.snr_suffix};
  for i=1:length(suffix_list)
    suffix = suffix_list{i};
    fname_in = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,parms.ext);
    fname_out = sprintf('%s/%s%s%s',parms.outdir,fstem,suffix,fext);
    fs_mri_convert(fname_in,fname_out,args{:});
    switch i
      case 1
        fname_mean = fname_out;
      case 2
        fname_std = fname_out;
      case 3
        fname_snr = fname_out;
    end
  end
else
  fname_snr  = fname_tmp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_cmd(cmd)
	[s,r] = mmil_unix(cmd);
	if s, error('cmd %s failed:\n%s',cmd,r); end
return;

