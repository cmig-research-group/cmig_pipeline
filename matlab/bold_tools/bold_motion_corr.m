function [fname_out,fname_motion_mat,fname_motion_txt,start_info] = bold_motion_corr(fname,varargin)
%function [fname_out,fname_motion_mat,fname_motion_txt,start_info] = bold_motion_corr(fname,[options])
%
% Purpose: Correct head motion for BOLD scans using
%   rigid-body volume registration between data and reference
%
% Usage:
%  [fname_out,fname_motion_mat,fname_motion_txt,start_info] = bold_motion_corr(fname,'key1', value1,...);
%
% Required Input
%  fname: full path file name of input 4D volume (mgh/mgz format)
%
% Optional Parameters:
%  'outdir': output directory
%    if empty, will set to directory containing fname
%    {default = []}
%  'outstem': output file stem
%    if empty, will set to file stem of input file anme
%    {default = []}
%  'fname_out': output file name
%    if empty, will add '_mc' to outstem
%    {default = []}
%  'fname_motion_mat': output file name of mat file containing registration
%    matrices M_ref_to_orig for each frame
%    if empty, will append '_motion.mat' to outstem
%    {default = []}
%  'fname_motion_txt': output file name of text file containing motion parameters
%    if empty, will append '_motion.1D' or '_motion.par' to outstem
%    {default = []}
%  'fname_ref': full path file name of reference volume (mgh/mgz format)
%    if empty, will use fname as reference
%    {default = []}
%  'fname_mask': name of mask volume
%    if empty, will generate one from fname_ref
%    {default = []}
%  'ref_frame': reference frame number (1-based)
%    if ref_frame = 0, will use mean across frames, excluding first skipTRs
%    {default = 1}
%  'skipTRs': number of TRs at beginning of scan to be ignored in calculating mean
%    only used if ref_frame = 0
%    {default = 0}
%  'method': method used for head motion estimation and correction
%     allowed values: 'rbreg', 'afni', 'fsl'
%       rbreg: bold_estimate_motion and rbreg_vol2vol_icorr_mask
%       afni: AFNI's 3dvolreg
%       fsl: FSL's mcflirt
%    {default = 'afni'}
%  'interpm': interpolation method
%     0 = nearest neighbor, 1 = linear, 2 = cubic
%     3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%     if resolution or slice position differs for fname_in and fname_ref
%       and method is 'afni' or 'fsl', will resample volume before motion correction
%    { default: 2 }
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': overwrite existing output files
%     {default: 0}
%
% Optional Parameters for motion estimation
%     Used only if method = 'rbreg'
%  'intra_flag': perform within-scan motion correction
%     {default = 1}
%  'min_trans': minimum translation (mm)
%     {default = 0.05}
%  'min_rot': minimum rotation (degrees)
%     {default = 0.05}
%
% Output:
%   fname_out: output file name of motion corrected 4D volume
%   fname_motion_mat: output file name of mat file with registration matrices
%   fname_motion_txt: output file name of text file containing motion parameters
%   start_info: start of the first registration and distance between the
%     two starts' results, from bold_estimate_motion (saved in fname_motion_mat)
%     empty if method is not 'rbreg' or fname_motion_mat lacks it
%
% Created:  04/26/10 by Don Hagler
% Prev Mod: 01/31/21 by Don Hagler
% Prev Mod: 02/04/21 by Don Hagler
% Prev Mod: 05/26/23 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 10/13/2026 by Don Hagler
% Last Mod: 11/12/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_motion_mat = []; fname_motion_txt = []; start_info = [];
if ~mmil_check_nargs(nargin,2), return; end

parms = check_input(fname,varargin);

if output_exists(parms)
  fname_out = parms.fname_out;
  fname_motion_mat = parms.fname_motion_mat;
  fname_motion_txt = parms.fname_motion_txt;
  start_info = load_start_info(parms.fname_motion_mat);
  return;
end

if strcmp(parms.method,'rbreg')
  [parms,fname_motion_mat,fname_motion_txt] = estimate_motion(parms);
  fname_out = correct_estimated_motion(parms);
  start_info = load_start_info(fname_motion_mat);
else
  parms = check_mismatch(parms);
  [fname_out,fname_motion_mat,fname_motion_txt] = correct_motion(parms);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
  ... % optional
    'outdir',[],[],...
    'outstem',[],[],...
    'fname_out',[],[],...
    'fname_motion_mat',[],[],...
    'fname_motion_txt',[],[],...
    'fname_ref',fname,[],...
    'ref_frame',1,[0 Inf],...
    'skipTRs',0,[0,Inf],...
    'fname_mask',[],[],...
    'interpm',2,[0:5],...
    'method','afni',{'rbreg','afni','fsl'},...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % hidden
    'tmpdir','tmp_mc',[],...
    'tmpext','.mgh',[],...
    'cleanupflag',true,[false true],...
    'motion_radius',50,[],... % for calculating distance from angle
    'vox2ras_maxdiff',1e-3',[1e-16,1],...
    'offset_maxdiff',20',[1e-16,1e16],...
  ... % specific to method = 'fsl'
   	'ext','.nii',{'.nii','.nii.gz'},... % for temporary nifti
  ... % specific to method = 'rbreg'
    'intra_flag',true,[false true],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'mstep',4,[],...
    'scales',[0 83 49 27 16 9 5 3 2 1],[],...
  ...
    'estimate_tags',{'intra_flag','fname_ref','ref_frame','skipTRs',...
                     'fname_mask','min_trans','min_rot','verbose',...
                     'mstep','scales'},[],...
    'resamp_tags',{'fname_out','save_reg_flag','fname_reg','qmat',...
      'fname_qmat','native_flag','nvoxels','resolution','deoblique_flag',...
      'std_orient','M_reg','M_ref','volsz_ref','regT1flag','fname_T1',...
      'M_T1_to_EPI','smooth','rot','trans','interpm','bclamp','EPI_type',...
      'ext','verbose','forceflag'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  if ~exist(fname,'file')
    error('file %s not found',fname);
  end
  if ~exist(parms.fname_ref,'file')
    error('file %s not found',parms.fname_ref);
  end

  [fpath,fstem,fext] = fileparts(parms.fname_in);
  if isempty(parms.outdir)
    if ~isempty(parms.fname_out)
      parms.outdir = fileparts(parms.fname_out);
    else
      parms.outdir = fpath;
    end
  end
  if mmil_isrelative(parms.tmpdir)
    parms.tmpdir = [parms.outdir '/' parms.tmpdir];
  end
  if isempty(parms.outstem)
    parms.outstem = fstem;
  end
  if isempty(parms.fname_out)
    parms.fname_out = [parms.outdir '/' parms.outstem '_mc' fext];
  end
  if isempty(parms.fname_motion_mat)
    parms.fname_motion_mat = [parms.outdir '/' parms.outstem '_motion.mat'];
  end
  if isempty(parms.fname_motion_txt)
    if strcmp(parms.method,'fsl')
      parms.fname_motion_txt = [parms.outdir '/' parms.outstem '_motion.par'];
    else
      parms.fname_motion_txt = [parms.outdir '/' parms.outstem '_motion.1D'];
    end
  end
  parms.fname_res = [parms.tmpdir '/' parms.outstem '_res' parms.tmpext];

  % get info about input volume
  [parms.M,parms.volsz] = ...
    mmil_load_mgh_info(parms.fname_in,parms.forceflag);
  parms.nframes = parms.volsz(4);

  % get info about reference volume
  [parms.M_ref,parms.volsz_ref] = ...
    mmil_load_mgh_info(parms.fname_ref,parms.forceflag);
    
  % check that ref_frame is not greater than num frames
  parms.ref_frame = min(parms.ref_frame,parms.volsz_ref(4));
  % if ref_frame=0, check that skipTRs < num frames - 1
  if ~parms.ref_frame
    parms.skipTRs = min(parms.skipTRs,parms.volsz_ref(4)-1);
  end
    
  if strcmp(parms.method,'fsl')
    if strcmp(parms.ext,'.nii')
      parms.precmd = sprintf('setenv FSLOUTPUTTYPE NIFTI; ');
    elseif strcmp(parms.ext,'.nii.gz')
      parms.precmd = sprintf('setenv FSLOUTPUTTYPE NIFTI_GZ; ');
    else
      parms.precmd = [];
    end
  else
    parms.precmd = [];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end
  if ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end
  if ~exist(parms.fname_motion_mat,'file')
    exist_flag = 0;
  end
  if ~exist(parms.fname_motion_txt,'file')
    exist_flag = 0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_mismatch(parms)
  % compare volsz and volsz_ref
  vszmatch = 1;
  if ~all(parms.volsz(1:3)==parms.volsz_ref(1:3))
    fprintf('%s: WARNING: volume dimensions do not match for %s and %s\n',...
      mfilename,parms.fname_in,parms.fname_ref);
    vszmatch = 0;
  end
  % compare M and M_ref, volsz and volsz_ref
  [voxmatch,rotmatch,offmatch,maxrot,maxoff] = fs_vox2ras_match(parms.M,parms.M_ref,...
                                                                parms.vox2ras_maxdiff,parms.offset_maxdiff);
  if ~voxmatch
    fprintf('%s: WARNING: vox2ras matrices do not match for:\n',mfilename);
    fprintf('%s\n',parms.fname_in);
    disp(parms.M);
    fprintf('%s\n',parms.fname_ref);
    disp(parms.M_ref);
    fprintf('    max rotation = %0.4f\n',maxrot);
    fprintf('    max offset = %0.4f\n',maxoff);
  end
  if vszmatch && voxmatch, return; end
  
  % resample input volume to match ref
  if ~exist(parms.fname_res,'file') || parms.forceflag
    mmil_mkdir(parms.tmpdir);
    if ~vszmatch || ~rotmatch
      fprintf('%s: resampling to reference...\n',mfilename);
      % resample to reference
      tparms = parms;
      tparms.fname_out = parms.fname_res;
      tparms.nvoxels = parms.volsz_ref(1:3);
      tparms.resolution = fs_voxsize(parms.M_ref);
      tparms.deoblique_flag = 0;
      % adjust for difference in offset
      if ~offmatch
        fprintf('%s: adjusting for offset...\n',mfilename);
        in_ctr_vox = (parms.volsz(1:3))/2+1;
        in_ctr = parms.M(1:3,:)*[in_ctr_vox 1]';
        ref_ctr_vox = (parms.volsz_ref(1:3))/2+1;
        ref_ctr = parms.M_ref(1:3,:)*[ref_ctr_vox 1]';
        offset = ref_ctr-in_ctr;
        M_adj = parms.M;
        M_adj(1:3,4) = M_adj(1:3,4) + offset;
        M_init = M_RAS_TO_LPH*parms.M*inv(M_RAS_TO_LPH*M_adj);
        tparms.M_reg = M_init;
      end
      args = mmil_parms2args(tparms,parms.resamp_tags);
      epi_resample_data(parms.fname_in,args{:});
    elseif ~voxmatch
      % replace M with M_ref so afni output will have M that matches M_ref
      fprintf('%s: replacing M with M_ref...\n',mfilename);
      vol = fs_load_mgh(parms.fname_in);
      fs_save_mgh(vol,parms.fname_res,parms.M_ref);
    end
  end
  parms.fname_in = parms.fname_res;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,fname_motion_mat,fname_motion_txt] = correct_motion(parms)
  fname_out = [];
  fname_motion_mat = [];
  fname_motion_txt = [];
  if ~exist(parms.fname_out,'file') ||...
     ~exist(parms.fname_motion_txt,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    mmil_mkdir(parms.tmpdir);
    if parms.volsz(4)==1 && strcmp(parms.fname_in,parms.fname_ref)
      if parms.verbose
        fprintf('%s: copying %s to %s\n',...
          mfilename,parms.fname_in,parms.fname_out);
      end
      % copy from mgz to mgh
      fs_copy_mgh(parms.fname_in,parms.fname_out);
      % create motion txt file with 0s
      fid = fopen(parms.fname_motion_txt,'wt');
      if fid<0
        error('failed to open %s for writing',parms.fname_motion_txt);
      end
      fprintf(fid,'0 0 0 0 0 0 0 1 1\n');
      fclose(fid);
    else
      % convert reference scan to nii
      [fpath_ref,fstem_ref,fext_ref] = fileparts(parms.fname_ref);
      fname_ref_tmp = [parms.tmpdir '/' fstem_ref '_ref.nii'];
      if parms.ref_frame>0
        if parms.verbose
          fprintf('%s: using frame %d from %s as motion reference for %s...\n',...
            mfilename,parms.ref_frame,parms.fname_ref,parms.fname_in);
        end
        conv_options = sprintf('-f %d',parms.ref_frame-1);
        fs_mri_convert(parms.fname_ref,fname_ref_tmp,'options',conv_options,...
          'forceflag',parms.forceflag);
      else
        if parms.verbose
          fprintf('%s: using mean volume of %s (excluding initial %d frames) as motion reference for %s...\n',...
            mfilename,parms.fname_ref,parms.skipTRs,parms.fname_in);
        end
        %% todo: also save mean reference in outdir as mgz?
        %%       otherwise waste time for each scan reg to ref?
        if ~exist(fname_ref_tmp,'file') || parms.forceflag
          [vol,M] = fs_load_mgh(parms.fname_ref);
          if parms.skipTRs>0
            vol = vol(:,:,:,parms.skipTRs+1:end);
          end
          vol_ref = mean(vol,4);
          fs_save_nifti(vol_ref,fname_ref_tmp,M);
        end
      end
      % convert input image to nii
      [fpath,fstem,fext] = fileparts(parms.fname_in);
      fname_in_tmp = [parms.tmpdir '/' fstem '.nii'];
      fs_mri_convert(parms.fname_in,fname_in_tmp,'forceflag',parms.forceflag);
      % correct for head motion
      fstem_tmp = [parms.tmpdir '/' fstem '_mc'];
      fname_out_tmp = [fstem_tmp '.nii'];
      
      if ~exist(fname_out_tmp,'file') | parms.forceflag
        if parms.verbose
          fprintf('%s: correcting motion for %s\n',mfilename,parms.fname_in);
        end
        if strcmp(parms.method,'afni')
          % AFNI_DONT_LOGFILE keeps 3dvolreg from appending to ~/.afni.log
          cmd = sprintf('env AFNI_DONT_LOGFILE=YES 3dvolreg -prefix %s',fname_out_tmp);
          % kernel named: the default changed to heptic in newer AFNI;
          %   -zpad 4 keeps Fourier shifts from wrapping across the edge slices
          cmd = sprintf('%s -Fourier -zpad 4',cmd);
          cmd = sprintf('%s -base ''%s''',cmd,fname_ref_tmp);
          cmd = sprintf('%s -dfile %s',cmd,parms.fname_motion_txt);
          cmd = sprintf('%s %s',cmd,fname_in_tmp);
        elseif strcmp(parms.method,'fsl')
          cmd = sprintf('%smcflirt -in %s -out %s -r %s -plots',...
                parms.precmd,fname_in_tmp,fstem_tmp,fname_ref_tmp);
        end
        if parms.verbose
          fprintf('%s: cmd = %s\n',mfilename,cmd);
        end
        [status,result,errmsg] = mmil_unix(cmd);
        if status
          error('registration failed:\n%s\n%s',result,errmsg);
        end
      end
      % convert nifti to mgh
      fs_mri_convert(fname_out_tmp,parms.fname_out);
      
      % copy par file to fname_motion_txt
      if strcmp(parms.method,'fsl')
        fname_par_tmp = [fstem_tmp '.par'];
        if ~exist(fname_par_tmp,'file')
          error('cannot find par file %s\n',fname_par_tmp);
        end
        if ~strcmp(fname_par_tmp,parms.fname_motion_txt)
          cmd = sprintf('cp -p %s %s',fname_par_tmp,parms.fname_motion_txt);
          [status,result] = mmil_unix(cmd);
          if status
            error('cp failed:\n%s\n%s',cmd,result);
          end
        end
      end
    end
  end

  % create M_motion from fname_motion_txt
  if ~exist(parms.fname_motion_mat,'file') || parms.forceflag
    if strcmp(parms.method,'afni')
      motion_tseries_orig = mmil_load_motion_1D(parms.fname_motion_txt);
    elseif strcmp(parms.method,'fsl')
      motion_tseries_orig = mmil_load_motion_par(parms.fname_motion_txt);
    end
    nframes = size(motion_tseries_orig,1);
    M_motion = cell(nframes,1);
    motion_tseries = zeros(nframes,6);
    for f=1:size(motion_tseries_orig,1)
      tmp_motion = motion_tseries_orig(f,:);
      trans = tmp_motion(1:3);
      rot   = tmp_motion(4:6);
      M_reg = mmil_construct_M('trans',trans,'rot',rot,...
        'order',{'rot','trans'},'rot_order',{'z','x','y'});
      M_motion{f} = inv(M_reg);
      % recalculate motion_tseries
      tmp_motion = mmil_M_mat2vec(M_reg);
      motion_tseries(f,:) = tmp_motion;
    end
    
    % calculate absolute difference in head position or rotation
    [mean_motion,mean_trans,mean_rot] = ...
      mmil_mean_motion(motion_tseries,parms.motion_radius);
    % calculate framewise displacement and motion stats
    [motion_stats,motion_fd] = mmil_motion_stats(motion_tseries,parms.motion_radius);      
    % save motion data to mat file
    save(parms.fname_motion_mat,'M_motion',...
      'motion_tseries_orig','motion_tseries',...
      'mean_motion','mean_trans','mean_rot',...
      'motion_fd','motion_stats');
  end
  
  % delete temporary files
  if parms.cleanupflag & exist(parms.tmpdir,'dir')
    cmd = sprintf('rm -r %s\n',parms.tmpdir);
    [status,result,errmsg] = mmil_unix(cmd);
    if status
      fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n%s\n',parms.tmpdir,result,errmsg);
    end
  end
  
  fname_out = parms.fname_out;
  fname_motion_txt = parms.fname_motion_txt;
  fname_motion_mat = parms.fname_motion_mat;  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,fname_motion_mat,fname_motion_txt] = estimate_motion(parms)
  fname_motion_mat = [];
  fname_motion_txt = [];
  
  % estimate motion for each frame
  if ~exist(parms.fname_motion_mat,'file') || parms.forceflag
    mmil_mkdir(parms.outdir);
    if ~isempty(parms.fname_mask)
      if ~exist(parms.fname_mask) || parms.forceflag
        % create brain mask
        if parms.verbose
          fprintf('%s: creating brain mask %s...\n',mfilename,parms.fname_mask);
        end
        cmig_create_brainmask(parms.fname_ref,...
            'fname_mask',parms.fname_mask,...
            'outdir',parms.outdir,'binarize_flag',1,...
            'verbose',parms.verbose,'forceflag',parms.forceflag);
      end
    end
    if parms.verbose
      fprintf('%s: estimating motion in %s...\n',mfilename,parms.fname_in);
      tic;
    end
    args = mmil_parms2args(parms,parms.estimate_tags);
    [M_motion,motion_tseries,start_info] = ...
      bold_estimate_motion(parms.fname_in,args{:});
    if parms.verbose, toc; end

    % calculate absolute difference in head position or rotation
    [mean_motion,mean_trans,mean_rot] = ...
      mmil_mean_motion(motion_tseries,parms.motion_radius);
    % calculate framewise displacement and motion stats
    [motion_stats,motion_fd] = mmil_motion_stats(motion_tseries,parms.motion_radius);      
    % save motion data to mat file
    save(parms.fname_motion_mat,'M_motion',...
      'motion_tseries',...
      'mean_motion','mean_trans','mean_rot',...
      'motion_fd','motion_stats','start_info');
  else
    load(parms.fname_motion_mat);
  end
  parms.M_motion = M_motion;
  fname_motion_mat = parms.fname_motion_mat;

  % save motion_tseries to motion_1D file
  if ~exist(parms.fname_motion_txt,'file') || parms.forceflag
    fid = fopen(parms.fname_motion_txt,'wt');
    if fid==-1
      error('failed to open %s for writing',parms.fname_motion_txt);
    end
    for f=1:size(motion_tseries,1)
      fprintf(fid,'%4d%s\n',f-1,sprintf('%10.4f',motion_tseries(f,:)));
    end
    fclose(fid);
  end
  fname_motion_txt = parms.fname_motion_txt;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = correct_estimated_motion(parms)
  if ~exist(parms.fname_out,'file') || parms.forceflag
    [vol_orig,M] = fs_load_mgh(parms.fname_in);
    % initialize output
    vol = zeros([parms.volsz_ref(1:3),parms.nframes],'single');
    if parms.verbose
      if parms.nframes>1
        plurstr = 's';
      else
        plurstr = [];
      end;      
      fprintf('%s: resampling %d frame%s...\n',...
        mfilename,parms.nframes,plurstr);
    end
    for f = 1:parms.nframes
      vol_tmp = squeeze(vol_orig(:,:,:,f));
      vol_tmp = mmil_resample_vol(vol_tmp,M,...
        'M_ref',parms.M_ref,'nvox_ref',parms.volsz_ref,...
        'M_reg',parms.M_motion{f},...
        'interpm',parms.interpm);
      vol(:,:,:,f) = vol_tmp;
    end
    fs_save_mgh(vol,parms.fname_out,parms.M_ref);
  end
  fname_out = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function start_info = load_start_info(fname_motion_mat)
  % start_info saved by estimate_motion, empty if absent
  start_info = [];
  if isempty(fname_motion_mat) || ~exist(fname_motion_mat,'file'), return; end
  vars = whos('-file',fname_motion_mat);
  if ~ismember('start_info',{vars.name}), return; end
  tmp = load(fname_motion_mat,'start_info');
  start_info = tmp.start_info;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
