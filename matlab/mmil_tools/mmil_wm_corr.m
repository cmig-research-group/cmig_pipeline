function fname_out = mmil_wm_corr(fname_in,varargin)
%function fname_out = mmil_wm_corr(fname_in,[options])
%
% Required Input:
%   fname_in: full path of input file (mgh/mgz format)
%
% Optional Parameters:
%  'fname_T1': input file name for T1-weighted image
%    if not supplied, fname_in is assumed to be T1-weighted
%    if supplied, fname_in is assumed to be T2-weighted
%      requiring registration to T1 and resampling of segmentation masks
%    {default = []}
%  'M_T1_to_T2': registration matrix
%    if supplied, will skip registration to T1
%    {default = []}
%  'fname_out': output file name for bias corrected image (mgh/mgz format)
%    {default = 'corr.mgz'}
%  'fname_wm': output file for white matter segmentation
%    {default = 'wm.mgz'}
%  'fname_bm': output file for brain mask
%    {default = 'bm.mgz'}
%  'fname_gm': output file for gray matter mask
%    {default = 'gm.mgz'}
%  'fname_wm_T1': output file for white matter segmentation
%    in T1-weighted space; only relevant if fname_T1 is supplied
%    {default = 'wm_T1.mgz'}
%  'fname_bm_T1': output file for brain mask
%    in T1-weighted space; only relevant if fname_T1 is supplied
%    {default = 'bm_T1.mgz'}
%  'fname_gm_T1': output file for gray matter mask
%    in T1-weighted space; only relevant if fname_T1 is supplied
%    {default = 'gm_T1.mgz'}
%  'fname_reg_T1': output file for registration matrix M_T1_to_T2
%    only relevant if fname_T1 is supplied and M_T1_to_T2 is not
%    {default = 'regT2toT1.mat'}
%  'erode_wm_flag': [0|1] erode white matter mask
%    only applies if fname_wm does not already exist (unless forceflag=1)
%    {default = 0}
%  'dilate_bm_flag': [0|1] dilate brain mask
%    only applies if fname_bm does not already exist (unless forceflag=1)
%    {default = 0}
%  'n3flag': [0|1] apply N3 bias correction as pre-normalization step
%    {default = 0}
%  'fname_n3': output file for N3 bias corrected image
%    only applies if n3flag=1
%    {default = 'n3.mgz'}
%  'vol_target': target volume for white matter intensity
%    may be a scalar value
%    if empty, will use median white matter intensity value
%    {default = []}
%  'ratio_range': range of values to constrain ratio
%    if empty, no clipping
%    {default = []}
%  'image_range': range of values to clip output image
%    if empty, no clipping
%    {default = []}
%  'outlier_thresh': exclude outliers from wm mask
%    difference from local max/min relative to mean wm intensity
%    {default = 0.1}
%  'outlier_edge_flag': [0|1] only exclude outlier voxels from
%    edge of white matter mask
%    {default = 1}
%  'lambda1': smoothing factor inside mask
%    {default = 20}
%  'lambda2': smoothing factor outside mask
%    {default = 1000}
%  'postsmooth_flag': [0|1] whether to apply post-smoothing
%      to the bias field ratio
%    {default = 1}
%  'atlasdir': full path of atlas directory
%    {default =  [getenv('MMPS_DIR') '/atlases']}}
%  'atlasname': name of atlas file (omit .mat extension)
%     full path or relative to atlasdir
%    {default =  'T1_Atlas/T1_atlas'}
%  'atlasdir_T1T2': name of atlas directory
%     for pre-registration to T1-T2 atlas
%     full path or relative to atlasdir
%    {default =  'T1T2'}
%  'tmpdir': temporary directory containing intermediate output
%    created if fname_wm or fname_bm or fname_gm are empty
%    {default = 'tmp_wm_corr'}
%  'cleanup_flag': [0|1] remove temporary files before quitting
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 1}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  11/27/15 by Don Hagler
% Prev Mod: 09/23/20 by Don Hagler
% Prev Mod: 10/08/20 by Don Hagler
% Prev Mod: 10/14/20 by Don Hagler
% Prev Mod: 10/20/20 by Don Hagler
% Last Mod: 08/13/25 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname_in,varargin);
fname_out = parms.fname_out;

% check output file
if check_output(parms), return; end;

% white matter segmentation using T1
create_wmseg(parms);

if parms.T2_flag
  % register T2 to T1
  if isempty(parms.M_T1_to_T2), parms = regT2toT1(parms); end;
  
  % resample segmentation to T2
  resT2_wmseg(parms);
end;

% perform bias correction
bias_corr(parms);

% optionally remove temporary directory
cleanup_tmpdir(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_in,options)
  % parse options
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
  ...
    'fname_T1',[],[],...
    'M_T1_to_T2',[],[],...
    'fname_out','corr.mgz',[],...
    'fname_wm',[],[],...
    'fname_bm',[],[],...
    'fname_gm',[],[],...
    'fname_wm_T1',[],[],...
    'fname_bm_T1',[],[],...
    'fname_gm_T1',[],[],...
    'fname_reg_T1',[],[],...
    'erode_wm_flag',false,[false true],...
    'dilate_bm_flag',false,[false true],...
    'n3flag',false,[false true],...
    'fname_n3',[],[],...
    'vol_target',[],[],...
    'ratio_range',[],[],...
    'image_range',[],[],...
    'outlier_thresh',0.1,[0,1],...
    'outlier_edge_flag',true,[false true],...
    'lambda1',20,[1e-10,1e10],...
    'lambda2',1e3,[1e-10,1e10],...
    'postsmooth_flag',true,[false true],...
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas',[],... 
    'atlasdir_T1T2','T1T2',[],...
    'tmpdir','tmp_wm_corr',[],...
    'cleanup_flag',true,[false true],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'fname_wm_masked',[],[],...
  ... % for N3 bias correction
    'estop',0.001,[],...
    'maxiter',100,[],...
  ... % for dilation of brain mask
    'smooth1',30,[],...
    'thresh1',0.05,[],...
    'smooth2',30,[],...
    'thresh2',0.05,[],...
    'smooth3',0,[],...
    'thresh3',0,[],...  
  ...
    'ext','.mgz',[],...
  ...
    'wmseg_tags',{'fname_wm','fname_bm','fname_gm','erode_wm_flag','dilate_bm_flag',...
                  'verbose','forceflag','estop','maxiter','nvoxels',...
                  'sigma','thresh_init','thresh','smooth1','thresh1',...
                  'smooth2','thresh2','smooth3','thresh3'},[],...
    'bias_tags',{'vol_target','ratio_range','image_range',...
                 'outlier_thresh','outlier_edge_flag',...
                 'lambda1','lambda2','postsmooth_flag'},[],...
  });

  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;

  % check fname_T1 (used for segmentation)
  parms.T2_flag = ~isempty(parms.fname_T1);
  if ~parms.T2_flag
    parms.fname_T1 = parms.fname_in;
  elseif ~exist(parms.fname_T1,'file')
    error('file %s not found',parms.fname_T1);
  end;
  
  % create temporary output directory
  if isempty(parms.fname_wm) || isempty(parms.fname_bm) ||...
     isempty(parms.fname_gm) ||...
   ((isempty(parms.fname_wm_T1) || isempty(parms.fname_bm_T1) ||...
     isempty(parms.fname_gm_T1) ||...
     (isempty(parms.M_T1_to_T2) &&...
      (isempty(parms.fname_reg_T1) || ~exist(parms.fname_reg_T1,'file'))))...
       && parms.T2_flag) ||...
     (isempty(parms.fname_n3) && parms.n3flag)
    if mmil_isrelative(parms.tmpdir)
      parms.tmpdir = [pwd '/' parms.tmpdir];
    end;
    mmil_mkdir(parms.tmpdir);
  else
    parms.tmpdir = [];
  end;
  
  % set names of intermediate files
  if isempty(parms.fname_wm)
    parms.fname_wm = sprintf('%s/wm%s',parms.tmpdir,parms.ext);
  end;
  if isempty(parms.fname_bm)
    parms.fname_bm = sprintf('%s/bm%s',parms.tmpdir,parms.ext);
  end;
  if isempty(parms.fname_gm)
    parms.fname_gm = sprintf('%s/gm%s',parms.tmpdir,parms.ext);
  end;
  if isempty(parms.fname_n3) && parms.n3flag
    parms.fname_n3 = sprintf('%s/n3%s',parms.tmpdir,parms.ext);
  end;
  if parms.T2_flag
    if isempty(parms.fname_wm_T1)
      parms.fname_wm_T1 = sprintf('%s/wm_T1%s',parms.tmpdir,parms.ext);
    end;
    if isempty(parms.fname_bm_T1)
      parms.fname_bm_T1 = sprintf('%s/bm_T1%s',parms.tmpdir,parms.ext);
    end;
    if isempty(parms.fname_gm_T1)
      parms.fname_gm_T1 = sprintf('%s/gm_T1%s',parms.tmpdir,parms.ext);
    end;
    if isempty(parms.fname_reg_T1)
      parms.fname_reg_T1 = sprintf('%s/regT2toT1.mat',parms.tmpdir);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output_exists = check_output(parms)
  output_exists = 0;
  if parms.forceflag, return; end;
  if ~exist(parms.fname_out,'file'), return; end;
  if ~isempty(parms.fname_wm_masked) &&...
     ~exist(parms.fname_wm_masked,'file')
    return;
  end;
  output_exists = 1;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_wmseg(parms)
  % create white matter segmentation
  if parms.T2_flag
    parms.fname_wm = parms.fname_wm_T1;
    parms.fname_bm = parms.fname_bm_T1;
    parms.fname_gm = parms.fname_gm_T1;
  end;
  args = mmil_parms2args(parms,parms.wmseg_tags);
  mmil_wmseg(parms.fname_T1,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = regT2toT1(parms)
  if ~exist(parms.fname_reg_T1,'file') || parms.forceflag
    fname_T1 = parms.fname_T1;
    fname_T2 = parms.fname_in;
    p = [];
    p.fname_T1_mask = parms.fname_bm_T1;
    p.outdir = parms.tmpdir;
    p.atlasdir = parms.atlasdir;
    p.atlasname = parms.atlasname;
    p.atlasdir_T1T2 = parms.atlasdir_T1T2;
    p.cleanup_flag = parms.cleanup_flag;
    p.forceflag = parms.forceflag;
    fprintf('%s: registering %s to %s...\n',...
      mfilename,fname_T2,fname_T1);
    a = mmil_parms2args(p);
    M_T1_to_T2 = mmil_reg_T2toT1(fname_T2,fname_T1,a{:});
    save(parms.fname_reg_T1,'M_T1_to_T2','fname_T1','fname_T2');
    parms.M_T1_to_T2 = M_T1_to_T2;
  else
    reginfo = load(parms.fname_reg_T1);
    parms.M_T1_to_T2 = reginfo.M_T1_to_T2;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resT2_wmseg(parms)
  segnames = {'wm','bm','gm'};
  M_reg = parms.M_T1_to_T2;
  M_T2 = []; volsz_T2 = [];
  % resample each segmentation volume from T1 to T2
  for i=1:length(segnames)
    segname = segnames{i};
    fname_seg_T1 = parms.(['fname_' segname '_T1']);
    fname_seg = parms.(['fname_' segname]);
    if ~exist(fname_seg,'file') || parms.forceflag
      fprintf('%s: resampling %s into alignment with %s...\n',...
        mfilename,fname_seg_T1,parms.fname_in);
      [vol,M] = fs_load_mgh(fname_seg_T1);
      if isempty(M_T2)
        [M_T2,volsz_T2] = mmil_load_mgh_info(parms.fname_in,...
          parms.forceflag,parms.tmpdir);
      end;
      [vol,M] = mmil_resample_vol(vol,M,...
        'M_ref',M_T2,'nvox_ref',volsz_T2(1:3),...
        'M_reg',inv(M_reg));
      vol = 1.0*(vol>0.5);
      fs_save_mgh(vol,fname_seg,M);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bias_corr(parms)
  if parms.n3flag
    if parms.verbose
      fprintf('%s: performing N3 bias correction...\n',mfilename);
    end;
    if ~exist(parms.fname_n3,'file') || parms.forceflag
      vol_brain = ctx_load_mgh(parms.fname_in);
      vol_brain = vol_correct_bias_field_n3(vol_brain,parms.estop,parms.maxiter);
      ctx_save_mgh(vol_brain,parms.fname_n3);
    end;
    parms.fname_in = parms.fname_n3;
  end;

  % perform bias correction using white matter segmentation
  if parms.verbose
    fprintf('%s: performing wm bias correction...\n',mfilename);
  end;
  [vol_in,M] = fs_load_mgh(parms.fname_in);
  vol_wm = fs_load_mgh(parms.fname_wm);
  vol_bm = fs_load_mgh(parms.fname_bm);
  vol_gm = fs_load_mgh(parms.fname_gm);
  vol_bm = 1.0 * (vol_bm > 0);
  args = mmil_parms2args(parms,parms.bias_tags);
  [vol_out,vol_ratio,vol_pred,vol_wm_masked] = ...
    mmil_corr_wm_bias_spsm(vol_in,vol_wm,vol_bm,vol_gm,args{:});
  fs_save_mgh(vol_out,parms.fname_out,M);
  if ~isempty(parms.fname_wm_masked)
    fs_save_mgh(vol_wm_masked,parms.fname_wm_masked,M);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanup_flag && ~isempty(parms.tmpdir) && exist(parms.tmpdir,'dir')
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,result] = mmil_unix(cmd);
    if status
      fprintf('%s: WARNING: failed to remove tmpdir:\n%s\n%s\n',...
        mfilename,cmd,result);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

