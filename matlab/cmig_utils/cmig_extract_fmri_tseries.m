function cmig_extract_fmri_tseries(fname_data,fspath,M_reg,outdir,varargin)
%function cmig_extract_fmri_tseries(fname_data,fspath,M_reg,outdir,[options])
%
% Purpose: extract average time series data from 3D+time fMRI volume
%   for subcortical ROIs, cortical parcels, and cortical vertices
%
% Usage:
%  cmig_extract_rsfmri_tseries(fname_data,fspath,M_reg,outdir)
%
% Required input:
%  fname_data: 4D volume in .nii or .nii.gz format
%  fspath: full path of FreeSurfer recon directory
%  M_reg: 4x4 registration matrix from T1 scan used for FreeSurfer to fMRI scan
%  outdir: output directory
%
% Optional input:
%  fs_version: FreeSurfer version number
%    e.g., 450 corresponds to FreeSurfer version 4.5.0
%  ico_coverage_thresh: fraction of frames in which an ico7 vertex must be
%    nonzero to be included in the coverage label used for ico smoothing
%    {default = 0.9}

% Created:  08/27/2025 by Don Hagler
% Prev Mod: 08/28/2025 by Don Hagler
% Prev Mod: 09/05/2025 by Don Hagler
% Prev Mod: 09/09/2025 by Don Hagler
% Prev Mod: 04/21/2026 by Don Hagler
% Prev Mod: 06/03/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 09/01/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on ABCD_Analyze_rsBOLD_Exam
%   Created:  07/03/17 by Don Hagler
%   Prev Mod: 03/11/19 by Feng Xue
%   Last Mod: 10/03/23 by Don Hagler
% Based on BOLD_MMIL_Resting_Analysis
%   Created:  03/06/12 by Don Hagler
%   Last Mod: 06/30/17 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,4), return; end

% parse input parameters and check for problems
parms = check_input(fname_data,fspath,M_reg,outdir,varargin);

% import data (convert to mgz)
if parms.verbose==2
  fprintf('%s: import_data...\n',mfilename);
end
parms = import_data(parms);

% write register.dat file (and scripts to check registration)
parms = write_regdat(parms);

% sample data to cortical surface
if parms.aparc_flag || parms.aparc2_flag || parms.fparc_flag || parms.ico_flag
  if parms.verbose==2
    fprintf('%s: paint_data...\n',mfilename);
  end
  parms = paint_data(parms);
end

% resample to ico
if parms.ico_flag || parms.fparc_flag
  if parms.verbose==2
    fprintf('%s: resample_ico...\n',mfilename);
  end
  parms = resample_ico(parms);
end

% downsample to lower ico orders
if parms.ico_flag && parms.ico_downsample_flag
  if parms.verbose==2
    fprintf('%s: downsample_ico...\n',mfilename);
  end
  parms = downsample_ico(parms);
end

% extract time series for surface and/or volume ROIs
if parms.aparc_flag || parms.aparc2_flag || parms.fparc_flag || parms.aseg_flag
  if parms.verbose==2
    fprintf('%s: extract_tseries...\n',mfilename);
  end
  parms = extract_tseries(parms);
end

% cleanup
if parms.cleanupflag, cleanup_tmpdir(parms); end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_data,fspath,M_reg,outdir,args)
  parms_filter = {...
    'fname_data',fname_data,[],...
    'fspath',fspath,[],...
    'M_reg',M_reg,[],...
    'outdir',outdir,[],...
  ...
    'fs_version',[],[],...
    'outstem',[],[],...
    'verbose',1,[0:2],...
    'cleanupflag',true,[false true],...
    'forceflag',false,[false true],...
  ... % for sampling to surface
    'projfrac_flag',false,[false true],...
    'projdist',1,[-10,10],...
    'projfrac',0.5,[-2,2],...
    'projdist_avg',[],[],...
    'projfrac_avg',[],[],...
    'mask_midbrain_flag',false,[false true],...
    'surfname','white',[],...
    'ico_flag',true,[false true],...
    'ico_downsample_flag',true,[false true],...
    'ico_order',[3,5],[0 7],...
    'ico_smooth',[64,16],[0,1000],...
    'ico_coverage_thresh',0.9,[0 1],...
  ... % for extracting ROI time courses
    'roi_outfix',[],[],...
    'aparc_flag',true,[false true],...
    'aparc_infix','aparc',[],...
    'aparc_roilist',[1001:1003,1005:1035,2001:2003,2005:2035],[1,Inf],...
    'aparc2_flag',true,[false true],...
    'aparc2_infix','aparc.a2009s',[],...
    'aparc2_roilist',[11101:11141,11143:11175,12101:12141,12143:12175],[1,Inf],...
    'fparc_flag',true,[false true],...
    'fnames_fparc',[],[],...
    'aseg_flag',true,[false true],...
    'aseg_infix','aseg',[],...
    'aseg_roilist',[2:5,7,8,10:18,24,26,28,41:44,46,47,49:54,58,60,77],[1,Inf],...
    'aseg_erode_flag',false,[false true],...
    'aseg_erode_nvoxels',1,[1:100],...
  ...
    'hemilist',{'lh','rh'},{'lh','rh'},...
    'outtype','mgz',{'mgh','mgz'},...
    'ico_nverts',[42 162 642 2562 10242 40962 163842],[],...
    'icosubj','fsaverage',[],...
    'icosubjdir',[getenv('FREESURFER_HOME') '/subjects'],[],...
  ...
    'paint_tags',{'regfile' 'projfrac_flag' 'projdist' 'projfrac'...
                  'projdist_avg' 'projfrac_avg' 'subjdir' 'surfname'...
                  'mask_midbrain_flag' 'outtype' 'forceflag'},[],...
    'extract_aseg_tags',{'outdir' 'regfile' 'subjdir' 'forceflag',...
                        'fname_aseg','aseg_roilist' 'FS_version'},[],...
    'extract_aparc_tags',{'outdir' 'subjdir' 'forceflag'...
                         'fnames_aparc' 'hemilist' 'fnames_surf' 'verbose'...
                         'projfrac_flag' 'projdist' 'projfrac' 'projdist_avg'...
                         'projfrac_avg' 'surfname'},[],...
  };
  parms = mmil_args2parms(args,parms_filter);

  if length(parms.ico_order) ~= length(parms.ico_smooth)
    error('ico_order and ico_smooth must have the same number of elements');
  end

  parms.nhemi = length(parms.hemilist);

  % check fname_data
  if isempty(parms.fname_data)
    error('fname_data is empty');
  elseif ~exist(parms.fname_data,'file')
    error('fname_data %s not found',parms.fname_data);
  end
  
  % check registration matrix
  if isempty(parms.M_reg)
    fprintf('%s: WARNING: M_reg is empty, assuming identity matrix\n',mfilename);
    parms.M_reg = eye(4);
  end
  
  % check outdir
  if mmil_isrelative(parms.outdir)
    parms.outdir = [pwd '/' parms.outdir];
  end
    
  % create tmpdir for intermediate files
  [~,parms.fstem] = get_fileparts(parms.fname_data);
  parms.tmpdir = sprintf('%s/%s',parms.outdir,parms.fstem);
  parms.tmpstem = sprintf('%s/%s',parms.tmpdir,parms.fstem);
  mmil_mkdir(parms.tmpdir)
  
  % set outstem if empty
  if isempty(parms.outstem)
    parms.outstem = parms.fstem;
  end
  if mmil_isrelative(parms.outstem)
    parms.outstem = [parms.outdir '/' parms.outstem];
  end

  % check fspath
  if isempty(parms.fspath)
    error('fspath is empty');
  elseif ~exist(parms.fspath,'file')
    error('fspath %s not found',parms.fspath);
  end
  % set subj name
  [parms.subjdir,parms.subj,text] = fileparts(parms.fspath);
  parms.subj = [parms.subj text];
  % check T1 file
  parms.fname_T1 = sprintf('%s/mri/nu.mgz',parms.fspath);
  if ~exist(parms.fname_T1,'file')
    error('fname_T1 %s not found',parms.fname_T1);
  end  
  
  % check fname_aseg
  if parms.aseg_flag
    parms.fname_aseg = sprintf('%s/mri/%s.mgz',parms.fspath,parms.aseg_infix);
    if ~exist(parms.fname_aseg,'file')
      error('aseg file %s not found',parms.fname_aseg);
    end
  else
    parms.fname_aseg = [];
  end

  % check aparc files
  if parms.aparc_flag
    parms.fnames_aparc = cell(1,parms.nhemi);
    for h=1:parms.nhemi
      parms.fnames_aparc{h} = sprintf('%s/label/%s.%s.annot',...
          parms.fspath,parms.hemilist{h},parms.aparc_infix);
      if ~exist(parms.fnames_aparc{h},'file')
        error('aparc annot file %s not found',parms.fnames_aparc{h});
      end
    end
  end

  % check aparc2 files
  if parms.aparc2_flag
    parms.fnames_aparc2 = cell(1,parms.nhemi);
    for h=1:parms.nhemi
      parms.fnames_aparc2{h} = sprintf('%s/label/%s.%s.annot',...
          parms.fspath,parms.hemilist{h},parms.aparc2_infix);
      if ~exist(parms.fnames_aparc2{h},'file')
        error('aparc annot file %s not found',parms.fnames_aparc2{h});
      end
    end
  end
    
  % check fparc files
  if parms.fparc_flag
    if isempty(parms.fnames_fparc)
      fparcdir = sprintf('%s/parms/rsfMRI_networks/network_gordon',cmig_getenv('CMIG_PIPELINE_DIR'));
      fparcname = 'gordon';
      parms.fnames_fparc = {sprintf('%s/lh.%s.annot',fparcdir,fparcname)...
                            sprintf('%s/rh.%s.annot',fparcdir,fparcname)};
    end
    if ~iscell(parms.fnames_fparc)
      parms.fnames_fparc = {parms.fnames_fparc};
    end
    parms.num_fparcs = length(parms.fnames_fparc);
    for f=1:parms.num_fparcs
      if ~exist(parms.fnames_fparc{f},'file')
        error('fparc annot file %s not found\n',parms.fnames_fparc{f});
      end
      n = regexp(parms.fnames_fparc{f},'(?<hemi>[lr]h)\.(?<name>.+)\.annot$','names');
      if isempty(n)
        error('unexpected naming convention for fparc file %s\n',parms.fnames_fparc{f});
      end
      parms.fparc_hemis{f} = n.hemi;
      parms.fparc_names{f} = n.name;
    end
  end

  if ~isempty(parms.roi_outfix)
    parms.roi_outstem = [parms.outstem '_' parms.roi_outfix];
  else
    parms.roi_outstem = parms.outstem;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = write_regdat(parms)
  % create fname_reg with RegInfo struct
  fname_regmat = sprintf('%s/%s_regT1.mat',parms.tmpdir,parms.fstem);
  if ~exist(fname_regmat,'file') || parms.forceflag
    [M_T2,volsz_T2] = mmil_load_mgh_info(parms.fname_data,parms.forceflag,parms.tmpdir);
    [M_T1,volsz_T1] = mmil_load_mgh_info(parms.fname_T1,parms.forceflag,parms.tmpdir);
    RegInfo = [];
    RegInfo.fname_T1 = parms.fname_T1;
    RegInfo.fname_T2 = parms.fname_data;
    RegInfo.M_T1_to_T2 = parms.M_reg;
    RegInfo.M_T1 = M_T1;
    RegInfo.M_T2 = M_T2;
    RegInfo.volsz_T1 = volsz_T1;
    RegInfo.volsz_T2 = volsz_T2;
    save(fname_regmat,'RegInfo');
  end
  % create register.dat file and scripts to check registration
  mmil_write_regscripts(fname_regmat,'outdir',parms.tmpdir,'forceflag',parms.forceflag);
  % set regfile
  parms.regfile = sprintf('%s/%s_register.dat',parms.tmpdir,parms.fstem);
  if ~exist(parms.regfile,'file')
    error('regfile %s not found',parms.regfile);
  end
  % convert register.dat to LTA for FreeSurfer 7.x+ compatibility
  % (mri_vol2vol fails with --reg on non-COR geometry volumes in FS 7.x and 8.x)
  % skip if lta_convert is not in PATH (e.g. FS 5.3.0)
  % NOTE: used for aseg resampling only; fs_paint parses the regfile itself
  %   with fs_read_regdat, which cannot read an LTA
  parms.regfile_aseg = parms.regfile;
  [s_which,~] = mmil_unix('which lta_convert');
  if s_which == 0
    fname_lta = regexprep(parms.regfile,'register\.dat$','register.lta');
    if ~exist(fname_lta,'file') || parms.forceflag
      cmd = sprintf('lta_convert --inreg %s --outlta %s --src %s --trg %s',...
        parms.regfile,fname_lta,parms.fname_data,parms.fname_T1);
      [s,r] = mmil_unix(cmd);
      if s, error('lta_convert failed:\n%s',r); end
    end
    parms.regfile_aseg = fname_lta;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = paint_data(parms)
  % sample 4D data to surface
  args = mmil_parms2args(parms,parms.paint_tags);
  parms.fname_data_surf = fs_paint(parms.subj,parms.fname_data,...
    'outstem',parms.tmpstem,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_aparc_tseries(parms,aparc2_flag)
  if ~exist('aparc2_flag','var') || isempty(aparc2_flag), aparc2_flag = 0; end
  if ~aparc2_flag
    suffix = '_aparc_data';
    aparc_roilist = parms.aparc_roilist;
    fnames_aparc = parms.fnames_aparc;
    aparc_name = 'aparc';
  else
    suffix = '_aparc2_data';
    aparc_roilist = parms.aparc2_roilist;
    fnames_aparc = parms.fnames_aparc2;
    aparc_name = 'aparc.a2009s';
  end
  parms.fname_aparc_data = sprintf('%s%s.mat',parms.roi_outstem,suffix);
  if ~exist(parms.fname_aparc_data,'file') || parms.forceflag
    % extract time series for aparc ROIs
    tparms = parms;
    tparms.outdir = parms.tmpdir;
    tparms.hemilist = parms.hemilist;
    tparms.fnames_surf = parms.fname_data_surf;
    tparms.fnames_aparc = fnames_aparc;
%    tparms.verbose = (parms.verbose==2);
    tparms.verbose = 0;
    if parms.verbose==2
      fprintf('%s: extracting time series from %s and %s for surface ROIs from %s and %s...\n',...
        mfilename,parms.fname_data_surf{1},parms.fname_data_surf{2},...
                  parms.fnames_aparc{1},parms.fnames_aparc{2});
    end
    args = mmil_parms2args(tparms,parms.extract_aparc_tags);
    [data,roinames,roicodes] = mmil_extract_aparc_tseries(parms.subj,[],args{:});
    if ~isempty(aparc_roilist)
      % include specified ROIs only
      ind_keep = find(ismember(roicodes,aparc_roilist));
    else
      % exclude ROIs with name 'unknown'
      ind_keep = find(cellfun(@isempty,regexp(lower(roinames),'unknown')));
    end
    data = data(:,ind_keep);
    roinames = roinames(ind_keep);
    roicodes = roicodes(ind_keep);
    % save results to mat file
    save(parms.fname_aparc_data,'data','roinames','roicodes','-v7');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_fparc_tseries(parms)
  parms.fname_fparc_data = sprintf('%s_fparc_data.mat',parms.roi_outstem);
  if ~exist(parms.fname_fparc_data,'file') || parms.forceflag
    % extract time series for fparc ROIs
    data = [];
    roinames = {};
    roicodes = [];
    for i=1:parms.num_fparcs
      hemi = parms.fparc_hemis{i};
      h = find(strcmp(hemi,parms.hemilist));
      tparms = parms;
      tparms.outdir = parms.tmpdir;
      tparms.hemilist = parms.hemilist(h);
      tparms.fnames_surf = parms.fname_data_ico(h);
      tparms.fnames_aparc = parms.fnames_fparc(i);
      tparms.subjdir = parms.icosubjdir;
      tparms.verbose = 0;
      if parms.verbose==2
        fprintf('%s: extracting time series from %s for surface ROIs from %s...\n',...
          mfilename,parms.fname_data_ico{h},parms.fnames_fparc{i});
      end
      args = mmil_parms2args(tparms,parms.extract_aparc_tags);
      [tmp_data,tmp_roinames,tmp_roicodes] = mmil_extract_aparc_tseries(parms.icosubj,[],args{:});
      % exclude ROIs with name 'unknown'
      ind_keep = find(cellfun(@isempty,regexp(lower(tmp_roinames),'unknown')));
      tmp_data = tmp_data(:,ind_keep);
      tmp_roinames = tmp_roinames(ind_keep);
      tmp_roicodes = tmp_roicodes(ind_keep);
      % combine data for multiple parcellations
      data = cat(2,data,tmp_data);
      roinames = cat(2,roinames,tmp_roinames);
      roicodes = cat(2,roicodes,tmp_roicodes);
    end
    save(parms.fname_fparc_data,'data','roinames','roicodes','-v7');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_aseg_tseries(parms)
  parms.fname_aseg_data = sprintf('%s_aseg_data.mat',parms.roi_outstem);
  if ~exist(parms.fname_aseg_data,'file') || parms.forceflag
    % extract time series for aseg ROIs
    if parms.verbose==2
      fprintf('%s: extracting volume ROI time series from %s...\n',...
        mfilename,parms.fname_data);
    end
    tparms = parms;
    tparms.outdir = parms.tmpdir;
    tparms.regfile = parms.regfile_aseg;
    tparms.erode_flag = parms.aseg_erode_flag;
    tparms.erode_nvoxels = parms.aseg_erode_nvoxels;
    tparms.FS_version = parms.fs_version;
    args = mmil_parms2args(tparms,parms.extract_aseg_tags);
    [data,roinames] = mmil_extract_aseg_tseries(parms.subj,parms.fname_data,args{:});
    %% NOTE: ROIs could possibly have no voxels left after eroding if small
    save(parms.fname_aseg_data,'data','roinames','-v7');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% todo: change this to produce concatenated format?
%%       or be run as part of concatenation?
function parms = combine_roi_tseries(parms)
  parms.fname_roi_data = sprintf('%s_roi_data.mat',parms.roi_outstem);
  if ~exist(parms.fname_roi_data,'file') || parms.forceflag
    if parms.verbose==2
      fprintf('%s: saving combined ROI time series data ...\n',mfilename);
    end
    % initialize combined data and roinames
    data = []; roinames = {};
    if parms.aparc_flag
      parms.fname_aparc_data = sprintf('%s_aparc_data.mat',parms.roi_outstem);
      aparc_data = mmil_load_matfile(parms.fname_aparc_data);
      data = cat(2,data,aparc_data.data);
      roinames = cat(2,roinames,aparc_data.roinames);
    end
    if parms.aparc2_flag
      parms.fname_aparc2_data = sprintf('%s_aparc2_data.mat',parms.roi_outstem);
      aparc2_data = mmil_load_matfile(parms.fname_aparc2_data);
      data = cat(2,data,aparc2_data.data);
      roinames = cat(2,roinames,aparc2_data.roinames);
    end
    if parms.fparc_flag
      parms.fname_fparc_data = sprintf('%s_fparc_data.mat',parms.roi_outstem);
      fparc_data = mmil_load_matfile(parms.fname_fparc_data);
      data = cat(2,data,fparc_data.data);
      roinames = cat(2,roinames,fparc_data.roinames);
    end
    if parms.aseg_flag
      parms.fname_aseg_data = sprintf('%s_aseg_data.mat',parms.roi_outstem);
      aseg_data = mmil_load_matfile(parms.fname_aseg_data);
      data = cat(2,data,aseg_data.data);
      roinames = cat(2,roinames,aseg_data.roinames);
    end
    % write mat file containing time series for each ROI  
    save(parms.fname_roi_data,'data','roinames','-v7');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = extract_tseries(parms)
  % extract time series for surface aparc ROIs (subject-space) (Desikan)
  if parms.aparc_flag
    parms = extract_aparc_tseries(parms);
  end
  % extract time series for surface aparc ROIs (subject-space) (Destrieux)
  if parms.aparc2_flag
    parms = extract_aparc_tseries(parms,parms.aparc2_flag);
  end
  % extract time series for surface fparc ROIs (atlas-space)
  if parms.fparc_flag
    parms = extract_fparc_tseries(parms);
  end
  % extract time series for volume ROIs
  if parms.aseg_flag
    parms = extract_aseg_tseries(parms);
  end
  % combine ROI time series data and save
  parms = combine_roi_tseries(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = resample_ico(parms)
  tparms = [];
  tparms.trgsubj = 'ico';
  tparms.icolevel = 7;
  tparms.smooth_in = 0;
  tparms.subjdir = parms.subjdir;
  tparms.intype = parms.outtype;
  tparms.forceflag = parms.forceflag;
%  tparms.verbose = (parms.verbose==2);
  tparms.verbose = 0;
  for h=1:parms.nhemi
    tparms.hemi = parms.hemilist{h};
    tparms.fname_out = sprintf('%s-ico%d-%s.%s',...
      parms.tmpstem,tparms.icolevel,tparms.hemi,parms.outtype);
    if ~exist(tparms.fname_out,'file') || parms.forceflag
      if parms.verbose==2
        fprintf('%s: resampling %s to ico...\n',...
          mfilename,parms.fname_data_surf{h});
      end
      args = mmil_parms2args(tparms);
      fs_surf2surf(parms.fname_data_surf{h},parms.subj,args{:});
    end
    parms.fname_data_ico{h} = tparms.fname_out;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = downsample_ico(parms)
  % coverage label for each hemisphere, from the unsmoothed ico7 time series
  fnames_label = cell(1,parms.nhemi);
  coverage = cell(1,parms.nhemi);
  for h=1:parms.nhemi
    [fnames_label{h},coverage{h}] = make_coverage_label(parms,h);
  end
  for i=1:length(parms.ico_order)
    ico_order = parms.ico_order(i);
    ico_smooth = parms.ico_smooth(i);
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fname_in = parms.fname_data_ico{h};

      % apply smoothing to ico tseries        
      tparms = [];
      tparms.hemi = hemi;
      tparms.trgsubj = parms.icosubj;
      tparms.icolevel = 7;
      tparms.smooth_in = ico_smooth;
      % smooth only within the coverage label (mri_surf2surf --label-src):
      %   each step averages only in-label neighbours, and values outside
      %   the label are set to 0; cortex_flag stays off, since --cortex
      %   replaces --label-src
      tparms.cortex_flag = 0;
      tparms.options = sprintf('--label-src %s',fnames_label{h});
      tparms.subjdir = parms.icosubjdir;
      tparms.intype = parms.outtype;
      tparms.forceflag = parms.forceflag;
      tparms.verbose = 0;
      tparms.fname_out = sprintf('%s-ico%d-sm%d-%s.%s',...
        parms.tmpstem,tparms.icolevel,ico_smooth,tparms.hemi,parms.outtype);
      if ~exist(tparms.fname_out,'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: smoothing %s with %d smoothing steps...\n',...
            mfilename,fname_in,ico_smooth);
        end
        args = mmil_parms2args(tparms);
        fs_surf2surf(fname_in,parms.icosubj,args{:});
      end
      fname = sprintf('fname_data_ico_sm%d',ico_order);
      parms.(fname){h} = tparms.fname_out;

      % create truncated copy with ico_order number of vertices
      fname_in = tparms.fname_out;
      fname_out = sprintf('%s-ico%d-%s.%s',...
        parms.outstem,ico_order,hemi,parms.outtype);
      if ~exist(fname_out,'file') || parms.forceflag
        if parms.verbose==2
          fprintf('%s: truncating %s to ico%d...\n',...
            mfilename,fname_in,ico_order);
        end
        nverts = parms.ico_nverts(ico_order);
        vals = fs_load_mgh(fname_in);
        vals = vals(1:nverts,:,:,:);
        fs_save_mgh(vals,fname_out);
      end
      fname = sprintf('fname_data_ico%d',ico_order);
      parms.(fname){h} = fname_out;

      % save coverage mask (1 covered, 0 not) truncated to ico_order
      fname_out = sprintf('%s-ico%d-mask-%s.%s',...
        parms.outstem,ico_order,hemi,parms.outtype);
      if ~exist(fname_out,'file') || parms.forceflag
        nverts = parms.ico_nverts(ico_order);
        fs_save_mgh(single(coverage{h}(1:nverts)),fname_out);
      end
      fname = sprintf('fname_mask_ico%d',ico_order);
      parms.(fname){h} = fname_out;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_label,coverage] = make_coverage_label(parms,h)
  % a vertex is covered when it is nonzero in at least ico_coverage_thresh of
  %   the frames of the unsmoothed ico7 time series and is in the icosubj
  %   ?h.cortex.label; vertices outside the bold slab are 0 in every frame
  hemi = parms.hemilist{h};
  fname_label = sprintf('%s-ico7-coverage-%s.label',parms.tmpstem,hemi);
  vals = fs_load_mgh(parms.fname_data_ico{h});
  nverts = size(vals,1);
  vals = reshape(vals,[nverts,numel(vals)/nverts]);
  frac = sum(vals~=0,2) / size(vals,2);
  clear vals;
  fname_cortex = sprintf('%s/%s/label/%s.cortex.label',...
    parms.icosubjdir,parms.icosubj,hemi);
  if ~exist(fname_cortex,'file')
    error('cortex label %s not found',fname_cortex);
  end
  cortex = false(nverts,1);
  cortex(fs_read_label(fname_cortex)) = true;
  coverage = (frac >= parms.ico_coverage_thresh) & cortex;
  fs_write_label(find(coverage),fname_label,parms.icosubj);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = import_data(parms)
  % get file stem
  [~,fstem,fext] = get_fileparts(parms.fname_data);
  if ~ismember(fext,{'.mgh','.mgz'})
    % set new name
    fname_tmp = sprintf('%s/%s.mgz',parms.tmpdir,fstem);
    % convert from nifti to mgz
    fs_mri_convert(parms.fname_data,fname_tmp);
    parms.fname_data = fname_tmp;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fpath,fstem,fext] = get_fileparts(fname)
  k = regexp(fname,'(?<ext>.nii.*)','names');
  if ~isempty(k)
    fname_tmp = regexprep(fname,'.nii.*','');
    [fpath,fstem,fext] = fileparts(fname_tmp);
    fstem = [fstem fext];
    fext = k.ext;
  else
    [fpath,fstem,fext] = fileparts(fname);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  fprintf('%s: removing temporary directory %s...\n',mfilename,parms.tmpdir);
  rmdir(parms.tmpdir,'s');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


