function [fname_out,fname_B0dx,fname_motion_txt,start_info] = bold_correct_data(fname,varargin)
%function [fname_out,fname_B0dx,fname_motion_txt,start_info] = bold_correct_data(fname,[options])
%
% Purpose: correct BOLD images for head motion, B0 inhomogeneity, and grad warp
%
% Required Parameters:
%   fname: full path name of mgh/mgz file containing 4D BOLD volume
%
% Optional Control Parameters:
%   'B0unwarp_flag': [0|1|2] whether to correct for B0 distortions
%     if 2, processing aborted if B0 unwarping fails (e.g. no rev scans)
%     {default = 1}
%   'regref_B0uw_flag': [0|1] whether to register reference image (fname_ref)
%      to fname_for (if revflag = 0) or fname_rev (if revflag = 1)
%      and apply transformation to B0dx
%     {default = 1}
%   'motion_B0uw_flag': [0|1] whether to estimate head motion
%     and apply to B0dx field
%     {default = 0}
%   'tshift_flag': [0|1] whether to correct for differences in slice timing
%     {default = 0}
%   'mc_flag': [0|1] whether to correct for within scan head motion
%     {default = 1}
%   'gradunwarp_flag': [0|1|2] whether to correct for grad warp
%     {default = 1}
%
% Optional Parameters:
%   'fname_out': output file fname for corrected data
%     if empty, will append '_corr.mgz' to file stem of fname
%     {default = []}
%   'verbose': [0|1] display status messages
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Optional Parameters for B0 unwarping:
%   'topup_flag': [0|1] whether to use FSL's topup to estimate distortion
%     otherwise uses Dominic Holland's unwarpB0 or epic
%     {default = 0}
%   'epic_flag': [0|1] whether to use epic binary instead of unwarpB0
%     {default = 0}
%   'epic_inorm_flag': [0|1] whether to rescale image intensities to match
%     Specific to epic, so ignored if epic_flag = 0
%     {default = 0}
%   'optimize_B0uw_flag': [0|1] search for optimal B0 unwarp parameters
%     kernelWidthMax and lambda2; ignored if topup_flag = 1
%     {default = 0}
%   'fname_B0dx'; name of input/output B0dx file with estimated displacements
%     if not found, will create using fname_for and fname_rev
%     {default = []}
%   'fname_for'; name of input b=0 image with "forward" phase-encode polarity
%     required for B0 unwarping if fname_B0dx not supplied or does not exist
%     {default = []}
%   'fname_rev'; name of input b=0 image with "reverse" phase-encode polarity
%      Required for B0 unwarping if fname_B0dx not supplied or does not exist
%     {default = []}
%   'fname_B0uw_ref': full path file name of reference volume (mgh/mgz format)
%     for registering to fname_B0dx
%     {default = []}
%   'fname_B0uw_mask': name of mask volume for registration to fname_B0uw_ref
%     If empty, will generate one from fname_ref using mmil_quick_brainmask
%     {default = []}
%   'EchoSpacing': echo spacing for fname
%     if supplied and different from EchoSpacing_for, will scale dx field
%     {default = []}
%   'AcquisitionColumns': acquisition columns for fname
%     if supplied and different from AcquisitionColumns_for, will scale dx field
%     {default = []}
%   'EchoSpacing_for': echo spacing for fname_for
%     {default = []}
%   'AcquisitionColumns_for': acquisition columns for fname_for
%     {default = []}
%   'EchoSpacing_rev': echo spacing for fname_rev
%     if supplied and different from EchoSpacing_for, will give warning
%     {default = []}
%   'AcquisitionColumns_rev': acquisition columns for fname_rev
%     if supplied and differenct from AcquisitionColumns_for, will give warning
%     {default = []}
%   'TotalReadoutTime': total readout time (s) for fname, for applytopup
%     if empty, will use TotalReadoutTime_for
%     {default = []}
%   'TotalReadoutTime_for': total readout time (s) for fname_for and
%     fname_rev, for topup; if empty, topup and applytopup are given 1 s
%     {default = []}
%   'PhaseDir': phase encode direction (obtained from dicom header)
%     'COL' or 'ROW'
%     {default = 'COL'}
%   'revflag': [0|1] whether to treat input data as "reverse" 
%     phase-encode polarity for applying B0 unwarping
%     {default = 0}
%
% Optional Parameters for estimating motion as part of B0uw:
%   'motion_B0uw_iters': number of iterations to estimate motion
%     and B0 displacement
%     {default = 2}
%   'min_trans': minimum translation (mm) for estimating motion for B0uw
%      {default = 0.05}
%   'min_rot': minimum rotation (degrees) for estimating motion for B0uw
%      {default = 0.05}
%   'fname_motion_B0uw': output file name of mat file containing estimated motion
%     if empty, will append '_motion_B0uw.mat' to file stem of fname
%     {default = []}
%   'fname_regref_B0uw': output file name of mat file containing registration
%      between fname_ref and fname_for or fname_ref
%     Ignored if regref_B0uw_flag = 0
%     If empty, will append '_regreg_B0uw.mat' to file stem of fname
%     {default = []}
%
% Optional Parameters for slice timing correction (uses AFNI's 3dTshift):
%   'tpattern': slice time pattern
%     allowed values: {'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z'}
%     {default = 'alt+z'}
%   'skipTRs': number of TRs at beginning of scan to be ignored in time shifting
%     {default = 0}
% 
% Optional Parameters for motion correction (uses AFNI's 3dvolreg):
%   'fname_ref': full path file name of reference volume (mgh/mgz format)
%     if empty, will use fname as reference
%     used if regref_B0uw_flag=1 or motion_B0uw_flag=1 or if mc_flag=1
%     {default = []}
%   'ref_frame': reference frame number (1-based)
%     used if motion_B0uw_flag=1 or if mc_flag=1
%     if ref_frame = 0, will use mean across frames, excluding first skipTRs
%     {default = 1}
%   'skipTRs': number of TRs at beginning of scan to be ignored in calculating mean
%     only used if ref_frame = 0
%    {default = 0}
%   'fname_motion_mat': output file name of mat file containing registration
%     matrices M_ref_to_orig for each frame
%     if empty, will append '_motion.mat' to file stem of fname
%     {default = []}
%   'fname_motion_txt': output file name of text file containing motion parameters
%     if empty, will append '_motion.1D' to file stem of fname
%     {default = []}
% 
% Optional Parameters for grad unwarping:
%   'gruw_type': gradient type number
%     0:  Siemens Sonata
%     1:  Siemens Allegra
%     2:  GE BRM
%     3:  GE CRM
%     4:  Siemens Avanto
%     5:  Siemens AXXess/Espree
%     6:  Siemens Quantum/Symphony
%     7:  GE Twin Speed Whole Body
%     8:  GE Twin Speed Zoom
%     9:  GE mr450 or mr750
%     10: GE MR750W
%     11: Siemens Skyra
%     12: Siemens Connectome Skyra
%     13: Siemens Prisma
%     14: GE Signa Premier
%     15: GE Signa UHP
%     16: Siemens Verio/Biograph
%     17: Siemens Vida
%     18: Siemens CimaX
%     {default = 8}
%   'gruw_unwarpflag': [0|1|2]
%      0: unwarp 3D
%      1: unwarp through plan only
%      2: unwarp inplane only
%      3: no unwarp needed
%      {default = 0}
%   'gruw_isoctrflag': [0|1] whether to adjust for isocenter coordinates
%      {default = 1}
%   'interpm': interpolation method for grad unwarp
%      0 = nearest neighbor, 1 = linear, 2 = cubic
%      3 = key's spline, 4 = cubic spline, 5 = hamming sinc
%      {default = 2}
% 
% Output:
%   fname_out: mgz file containing 4D BOLD volume corrected for
%     slice timing, head motion, B0 distortion, and gradient warping
%   fname_motion_txt: output file name of text file containing motion parameters
%   fname_B0dx: mgz file containing volume of estimated displacement
%     along phase encode direction
%   start_info: start of the between-scan registration and distance between
%     the two starts' results (see bold_estimate_motion)
%     empty if there was no rbreg between-scan registration
%
% Created:  04/27/2010 by Don Hagler
% Prev Mod: 11/06/2018 by Feng Xue
% Prev Mod: 10/01/2021 by Don Hagler
% Prev Mod: 05/30/2023 by Don Hagler
% Prev Mod: 12/09/2024 by Don Hagler
% Prev Mod: 06/04/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 08/20/2026 by Don Hagler
% Prev Mod: 10/09/2026 by Don Hagler
% Last Mod: 11/12/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = []; fname_B0dx = []; fname_motion_txt = []; start_info = [];

if ~mmil_check_nargs(nargin, 1), return; end;

parms = check_input(fname,varargin);

fname_out = parms.fname_out;
fname_B0dx = parms.fname_B0dx;
fname_motion_txt = parms.fname_motion_txt;

if output_exists(parms)
  start_info = load_start_info(parms);
  return;
end;

% estimate distortions due to magnetic field susceptibility inhomogeneities
if parms.B0unwarp_flag
  estimate_B0_distortion(parms);
end;

% correct for slice timing differences
if parms.tshift_flag
  parms = correct_slice_timing(parms);
end;

% correct for head motion
if parms.mc_flag && ~parms.motion_B0uw_flag
  parms = correct_motion(parms,1);
end;

% correct for distortion due to magnetic field susceptibility inhomogeneities
if parms.B0unwarp_flag
  parms = correct_B0_distortion(parms);
end;

% correct for distortion due to gradient coil inhomogeneities
if parms.gradunwarp_flag && parms.gruw_unwarpflag~=3
  parms = correct_gradwarp(parms);
end;

% correct for head motion
if parms.mc_flag || ~isempty(parms.fname_ref)
  parms = correct_motion(parms);
end;

% copy last file to fname_out
if ~exist(parms.fname_out,'file') || parms.forceflag
  fs_copy_mgh(parms.fname_in,parms.fname_out);
end;

% start of the between-scan registration
start_info = load_start_info(parms);

% remove temporary files
cleanup_tmpdir(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms_filter = {...
    'fname_in',fname,[],...
    'fname_orig',fname,[],...
  ...
    'B0unwarp_flag',1,[0:2],...
    'regref_B0uw_flag',true,[false true],...
    'mask_thresh',0.82,[0,1],...
    'motion_B0uw_flag',false,[false true],...
    'tshift_flag',false,[false true],...
    'mc_flag',true,[false true],...
    'gradunwarp_flag',1,[0:2],...
  ...
    'fname_out',[],[],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ... % B0 unwarping
    'topup_flag',false,[false true],...
    'epic_flag',false,[false true],...
    'epic_inorm_flag',false,[false true],...
    'optimize_B0uw_flag',false,[false true],...
    'fname_B0dx',[],[],...
    'fname_for',[],[],...
    'fname_rev',[],[],...
    'fname_B0uw_ref',[],[],...
    'fname_B0uw_mask',[],[],...
    'revflag',false,[false,true],...
    'EchoSpacing',[],[],...
    'AcquisitionColumns',[],[],...
    'AcquisitionRows',[],[],...
    'PhaseDir','COL',{'ROW','COL'},...
    'EchoSpacing_for',[],[],...
    'AcquisitionColumns_for',[],[],...
    'AcquisitionRows_for',[],[],...
    'EchoSpacing_rev',[],[],...
    'AcquisitionColumns_rev',[],[],...
    'AcquisitionRows_rev',[],[],...
    'TotalReadoutTime',[],[],...
    'TotalReadoutTime_for',[],[],...
  ... % motion correction for B0uw
    'motion_B0uw_iters',2,[1:10],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'fname_motion_B0uw',[],[],...
    'fname_regref_B0uw',[],[],...
  ... % slice timing correction
    'tpattern','alt+z',{'alt+z', 'alt+z2', 'alt-z', 'alt-z2', 'seq+z', 'seq-z','unknown'},...
    'skipTRs',0,[0,Inf],...
    'ts_suffix','ts',[],...
  ... % motion correction
    'fname_ref',[],[],...
    'ref_frame',1,[0 Inf],...
    'fname_motion_mat',[],[],...
    'fname_motion_txt',[],[],...
    'mc_suffix','mc',[],...
  ... % grad unwarp
    'gruw_type',8,[0:18],...
    'gruw_unwarpflag',0,[],...
    'gruw_isoctrflag',true,[false true],...
    'interpm',2,[0:5],...
    'gruw_jacobian_flag',false,[false true],...
  ... % B0uw optimization
    'kernelWidthMax',25,[1:100],...
    'lambda2',1100,[1:10000],...
    'kernelWidthMax_vec',[25,31,35],[1:100],...
    'lambda2_vec',[1100,1500,1900],[1:10000],...
    'multi_opt_flag',false,[false true],...
  ... % other
    'outfix','corr',[],...
    'outext','.mgz',{'.mgh','.mgz'},...
    'tmpext','.mgh',{'.mgh','.mgz'},...
    'cleanupflag',true,[false true],...
  ...
    'unwarpB0_tags',{'fname_out','verbose','forceflag',...
      'topup_flag','epic_flag','epic_inorm_flag','optimize_B0uw_flag',...
      'fname_B0dx','fname_for','fname_rev','revflag','EchoSpacing',...
      'AcquisitionColumns','AcquisitionRows','PhaseDir','EchoSpacing_for',...
      'AcquisitionColumns_for','AcquisitionRows_for','EchoSpacing_rev',...
      'AcquisitionColumns_rev','AcquisitionRows_rev',...
      'TotalReadoutTime','TotalReadoutTime_for',...
      'regref_B0uw_flag','mask_thresh','motion_B0uw_flag',...
      'motion_B0uw_iters','fname_ref','ref_frame','min_trans','min_rot',...
      'fname_motion_B0uw','fname_regref_B0uw',...
      'kernelWidthMax','lambda2','kernelWidthMax_vec',...
      'lambda2_vec','multi_opt_flag','outfix','outext','tmpext',...
      'cleanupflag'},[],...
    'tshift_tags',{'tpattern','skipTRs','cleanupflag','forceflag'},[],...
    'motion_tags',{'fname_out','fname_motion_mat','fname_motion_txt',...
                   'method','interpm','verbose','forceflag',...
                   'intra_flag','fname_ref','ref_frame','skipTRs','fname_mask',...
                   'min_trans','min_rot','tmpdir','tmpext','cleanupflag',...
                   'mstep','scales','motion_radius'},[],...
  };
  parms = mmil_args2parms(options,parms_filter);
  % check input file
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  % create output file fnames
  [tpath,fstem,~] = fileparts(parms.fname_in);
  if isempty(parms.fname_out)
    parms.fname_out = [tpath '/' fstem '_' parms.outfix parms.outext];
  end;
  parms.outdir = fileparts(parms.fname_out);
  parms.tmpdir = [parms.outdir '/tmp_BOLD_corr'];
  if isempty(parms.fname_motion_txt)
    parms.fname_motion_txt = [parms.outdir '/' fstem '_' parms.outfix '_motion.1D'];
  end;
  if isempty(parms.fname_motion_mat)
    parms.fname_motion_mat = [parms.outdir '/' fstem '_' parms.outfix '_motion.mat'];
  end;
  % disable tshift if number of TRs < skipTRs
  [M,volsz] = mmil_load_mgh_info(parms.fname_in,parms.forceflag,parms.outdir);
  if volsz(4) <= parms.skipTRs + 2, parms.tshift_flag = false; end;
  % check input files for B0 distortion correction
  if parms.B0unwarp_flag
    if isempty(parms.fname_B0dx) && ...
       ~isempty(parms.fname_for) && ~isempty(parms.fname_rev)
      [~,fstem_rev,~] = fileparts(parms.fname_rev);
      [~,fstem_for,~] = fileparts(parms.fname_for);
      parms.fname_B0dx = ...
        [parms.outdir '/' fstem_for 'VS' fstem_rev '_B0dx' parms.outext];
    end;
    if isempty(parms.fname_B0dx)
      if parms.B0unwarp_flag==2
        error('files required for B0 distortion correction not supplied');
      else
        fprintf('%s: WARNING: skipping B0 distortion correction\n',mfilename);
        parms.B0unwarp_flag = 0;
      end;
    else
      % check that fname_for and fname_rev exist
      if ~exist(parms.fname_B0dx,'file')
        if ~exist(parms.fname_for,'file')
          error('file %s not found',parms.fname_for);
        end;
        if ~exist(parms.fname_rev,'file')
          error('file %s not found',parms.fname_rev);
        end;
      end;
    end;
  end;
  if parms.verbose
    fprintf('%s: input file %s\n',mfilename,parms.fname_in);
    if parms.tshift_flag
      if strcmp(parms.tpattern,'unknown')
        fprintf('%s: skipping slice timing correction between tpattern is %s\n',...
          mfilename,parms.tpattern);
      else
        fprintf('%s: correcting for slice timing using tpattern %s\n',...
          mfilename,parms.tpattern);
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exist_flag = output_exists(parms)
  if parms.forceflag
    exist_flag = 0;
    return;
  else
    exist_flag = 1;
  end;
  if ~isempty(parms.fname_out) && ~exist(parms.fname_out,'file')
    exist_flag = 0;
  end;
  if parms.B0unwarp_flag && ~epi_topup_check_B0dx(parms.fname_B0dx,parms.topup_flag)
    exist_flag = 0;
  end;
  if parms.mc_flag && (~exist(parms.fname_motion_txt,'file') ||...
                       ~exist(parms.fname_motion_mat,'file'))
    exist_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estimate_B0_distortion(parms)
  args = mmil_parms2args(parms,parms.unwarpB0_tags);
  bold_unwarpB0([],args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_B0_distortion(parms)
  [~,fstem,~] = fileparts(parms.fname_in);
  tparms = parms;
  tparms.fname_out = [parms.tmpdir '/' fstem '_B0uw' parms.tmpext];
  if ~exist(tparms.fname_out,'file')
    if parms.verbose
      fprintf('%s: correcting for B0 distortions...\n',mfilename);
    end;
    if parms.motion_B0uw_flag && isempty(parms.fname_motion_B0uw)
      tparms.fname_motion_B0uw = [parms.outdir '/' fstem '_motion_B0uw.mat'];
    end;
    if parms.regref_B0uw_flag && isempty(parms.fname_regref_B0uw)
      tparms.fname_regref_B0uw = [parms.outdir '/' fstem '_regref_B0uw.mat'];
    end;
    if parms.regref_B0uw_flag
      % register to within-scan reference
      tparms.fname_ref = [];
      tparms.fname_mask = [];
    elseif parms.motion_B0uw_flag
      % register to scans used to estimate B0 distortion
      [tparms.fname_ref,tparms.fname_mask] = ...
        set_fname_ref(parms);
    end;
    % do not do between-scan registration for same scan
    if strcmp(tparms.fname_ref,parms.fname_orig)
      tparms.fname_ref = tparms.fname_in;
    end;
    args = mmil_parms2args(tparms,parms.unwarpB0_tags);
    bold_unwarpB0(parms.fname_in,args{:});
  end;
  parms.fname_in = tparms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_ref,fname_mask] = set_fname_ref(parms)
  if ~isempty(parms.fname_B0uw_ref)
    fname_ref = parms.fname_B0uw_ref;
    fname_mask = parms.fname_B0uw_mask;
  else
    if ~parms.revflag
      fname_ref = parms.fname_for;
    else
      fname_ref = parms.fname_rev;
    end;
    fname_mask = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_slice_timing(parms)
  [~,fstem,~] = fileparts(parms.fname_in);
  fname_out = [parms.tmpdir '/' fstem '_' parms.ts_suffix parms.tmpext];
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: correcting for slice timing...\n',mfilename);
      tic;
    end;
    args = mmil_parms2args(parms,parms.tshift_tags);
    fname_out = bold_tshift(parms.fname_in,'fname_out',fname_out,args{:});
    if parms.verbose, toc; end;
  end;
  parms.fname_in = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_motion(parms,preB0uw_flag)
  if ~exist('preB0uw_flag','var') || isempty(preB0uw_flag)
    preB0uw_flag = 0;
  end;
  %[tpath,fstem,text] = fileparts(parms.fname_in);
  [~,fstem,~] = fileparts(parms.fname_in);
  tparms = parms;
  tparms.fname_out = [parms.tmpdir '/' fstem '_' parms.mc_suffix parms.tmpext];
  tparms.tmpdir = [parms.tmpdir '/tmp_mc'];
  if ~exist(tparms.fname_out,'file') || parms.forceflag
    if preB0uw_flag
      if parms.regref_B0uw_flag
        % register to within-scan reference
        tparms.fname_ref = [];
        tparms.fname_mask = [];
      else
        % register to scans used for estimating B0 distortion
        [tparms.fname_ref,tparms.fname_mask] = ...
          set_fname_ref(parms);
      end;
    elseif ~parms.motion_B0uw_flag
      %[tpath,fstem_motion,text] = fileparts(parms.fname_motion_mat);
      [~,fstem_motion,~] = fileparts(parms.fname_motion_mat);
      tparms.fname_motion_mat = [parms.outdir '/' fstem_motion '_inter.mat'];
      tparms.fname_motion_txt = [parms.outdir '/' fstem_motion '_inter.1D'];
    end;
    % do not do between-scan registration for same scan
    if strcmp(tparms.fname_ref,parms.fname_orig)
      tparms.fname_ref = tparms.fname_in;
    end;
    if ~parms.mc_flag || (~parms.motion_B0uw_flag && ~preB0uw_flag)
      % between-scan registration and resampling only
      tparms.intra_flag = 0;
      tparms.method = 'rbreg';
      % NOTE: afni's 3dvolreg does not handle well substantial inter-scan motion
    else
      tparms.intra_flag = 1;
      tparms.method = 'afni';
      % NOTE: rbreg is very slow
    end;
    if parms.verbose
      if tparms.intra_flag
        fprintf('%s: correcting for head motion...\n',mfilename);
      else
        fprintf('%s: correcting for head motion (between scan only)...\n',...
          mfilename);
      end;
      tic;
    end;
    args = mmil_parms2args(tparms,parms.motion_tags);
    bold_motion_corr(parms.fname_in,args{:});
    if parms.verbose, toc; end;
  end;
  parms.fname_in = tparms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function start_info = load_start_info(parms)
  % start_info saved by bold_motion_corr for the final motion correction
  start_info = [];
  if ~parms.mc_flag && isempty(parms.fname_ref), return; end;
  if parms.motion_B0uw_flag
    fname_mat = parms.fname_motion_mat;
  else
    [~,fstem_motion,~] = fileparts(parms.fname_motion_mat);
    fname_mat = [parms.outdir '/' fstem_motion '_inter.mat'];
  end;
  if ~exist(fname_mat,'file'), return; end;
  vars = whos('-file',fname_mat);
  if ~ismember('start_info',{vars.name}), return; end;
  tmp = load(fname_mat,'start_info');
  start_info = tmp.start_info;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = correct_gradwarp(parms)
  %[tpath,fstem,text] = fileparts(parms.fname_in);
  [~,fstem,~] = fileparts(parms.fname_in);
  fname_out = [parms.tmpdir '/' fstem '_gruw' parms.tmpext];
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: correcting for gradient warping...\n',mfilename);
      tic
    end;
    epi_gradunwarp(parms.fname_in,...
      'fname_out',fname_out,...
      'gwtype',parms.gruw_type,...
      'unwarpflag',parms.gruw_unwarpflag,...
      'isoctrflag',parms.gruw_isoctrflag,...
      'jacobian_flag',parms.gruw_jacobian_flag,...
      'interpm',parms.interpm,...
      'forceflag',parms.forceflag);
    if parms.verbose, toc; end;
  end;
  parms.fname_in = fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanupflag && exist(parms.tmpdir,'dir')
    % delete tmp dir
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,result,errmsg] = jsystem(cmd);
    if status
      fprintf('%s: WARNING: deletion of temporary directory %s failed:\n%s\n',...
        mfilename,parms.tmpdir,errmsg);
    end;
  end;
return;

