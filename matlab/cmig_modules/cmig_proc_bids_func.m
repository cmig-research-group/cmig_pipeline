function [outdir,outdir_fmap,errcode,errmsg] = cmig_proc_bids_func(participant_id,session_id,input_rootdir,output_rootdir,varargin)
%function [outdir,outdir_fmap,errcode,errmsg] = cmig_proc_bids_func(participant_id,session_id,input_rootdir,output_rootdir,[options])
%
% Purpose: Process functional MRI and field map data from BIDS format
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   input_rootdir: BIDS rawdata root directory
%   output_rootdir: BIDS output root directory
%
% Optional Input ('key',value):
%   'B0unwarp_flag': [0|1|2] whether to correct for B0 distortions
%     if 0, do not perform B0 unwarping
%     if 1, perform B0unwarp if field map scans are available
%     if 2, quit if field map scans are unavailable
%     {default = 1}
%   'gradunwarp_flag': [0|1|2] whether to correct for gradient warping
%     if 0, do not perform gradwarp
%     if 1, perform gradwarp if gradwarp type available and supported
%     if 2, quit if gradwarp type unavailable or unsupported
%     {default = 2}
%   'native_flag': keep original nvoxels and resolution
%     {default = 0}
%   'nvoxels': vector of number of resampled voxels [nx,ny,nz]
%     ignored if native_flag = 1
%     {default = [90 90 60]}
%   'resolution': desired voxel sizes of resampled data (x y z)
%     ignored if native_flag = 1
%     {default = [2.4 2.4 2.4]}
%   'std_orient': output orientation of the func and fmap outputs, e.g.
%     'LAS' (dcm2niix's, for ABCD) or 'LPS' (MMPS's); empty keeps the input's.
%     PhaseEncodingDirection and SliceEncodingDirection in the output json
%     are restated for the output's voxel axes
%     ignored if native_flag = 1
%     {default = 'LAS'}
%   'regT1flag': [0,1] whether to register func data to T1
%     only compute registration, does not apply registration to output
%     set to 0, with a WARNING, when the session has no processed T1w
%       (cmig_proc_bids_anat is not run from here)
%     {default = 1}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory for processed functional scans
%     will have form of {output_rootdir}/{participant_id}/{session_id}/func
%     each scan's output is named as in rawdata with desc-preproc before the
%       suffix (_desc-preproc_bold.nii.gz, .json), its motion estimates
%       _desc-motion_timeseries.tsv
%   outdir_fmap: output directory for processed field map scans
%     will have form of {output_rootdir}/{participant_id}/{session_id}/fmap
%     named as the func outputs (_desc-preproc_epi.nii.gz, .json)
%     also contains one B0dx displacement field (.nii.gz, .json) per
%       forward/reverse pair, in mm, in the reference field map's frame at
%       output resolution, signed along that frame's phase-encode array axis
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  09/23/2025 by Don Hagler
% Prev Mod: 10/10/2025 by Don Hagler
% Prev Mod: 10/14/2025 by Don Hagler
% Prev Mod: 03/12/2026 by Don Hagler
% Prev Mod: 03/19/2026 by Don Hagler
% Prev Mod: 04/24/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 06/01/2026 by Don Hagler
% Prev Mod: 06/02/2026 by Don Hagler
% Prev Mod: 06/04/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Prev Mod: 06/12/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/12/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/19/2026 by Don Hagler
% Prev Mod: 08/26/2026 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/25/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/09/2026 by Don Hagler
% Prev Mod: 10/19/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Prev Mod: 11/09/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/12/2026 by Don Hagler
% Prev Mod: 11/13/2026 by Don Hagler
% Last Mod: 11/14/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on cmig_proc_bids_anat         
%   Created:  08/04/2025 by Don Hagler 
%   Last Mod: 09/04/2025 by Don Hagler 
% and MMIL_Process_BOLD
%   Created:  06/25/08 by Don Hagler
%   Last Mod: 05/30/23 by Don Hagler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   'regT1flag': [0,1] whether to register func data to T1
%     0=do not register, 1=register only
%     {default = 1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; outdir_fmap = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,4), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

% check parameters
parms = check_parms(participant_id,session_id,input_rootdir,output_rootdir,varargin);

% locate input data
[parms,errcode,errmsg] = check_input(parms);
if errcode
  % write status file
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-proc',errcode,errmsg,[],startdate);
  return;
end

% check whether processing required
[parms,procflag] = check_output(parms);

% process data
if procflag, [parms,errcode,errmsg] = proc_data(parms); end

% cleanup tmpdirs
if parms.cleanup_flag
  cleanup_errmsg = cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir);
  if isempty(errmsg), errmsg = cleanup_errmsg; end
end
if parms.cleanup_flag
  cleanup_errmsg = cmig_cleanup_tmpdirs(parms.outdir_fmap,parms.tmpdir_fmap);
  if isempty(errmsg), errmsg = cleanup_errmsg; end
end

% write status file
% an events-only copy leaves an existing status as it is, so its procdate
%   does not make postqc_func and postqc_fmap stale; analyze_func, which reads
%   the events copies, is checked against their mtime instead
statusflag = procflag && (errcode || parms.proc_fmap_flag || parms.proc_func_flag || parms.proc_b0dx_flag);
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-proc',errcode,errmsg,[],startdate,statusflag);

outdir = parms.outdir;
outdir_fmap = parms.outdir_fmap;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_parms(participant_id,session_id,input_rootdir,output_rootdir,options)
  options = cmig_parse_deployed_args(options);
  parms = mmil_args2parms(options,{...
    'participant_id',participant_id,[],...
    'session_id',session_id,[],...
    'input_rootdir',input_rootdir,[],...
    'output_rootdir',output_rootdir,[],...
  ...
    'B0unwarp_flag',1,[0 1 2],...
    'gradunwarp_flag',2,[0 1 2],...
    'native_flag',false,[false true],...
    'nvoxels',[90 90 60],[],...
    'resolution',[2.4 2.4 2.4],[],...
    'std_orient','LAS',[],...
    'regT1flag',1,[0 1],...
    'cleanup_flag',true,[false true],... % whether tmpdirs are deleted at end
    'forceflag',false,[false true],...
  ... % parameters for cmig_proc_bids_anat
    'anat_select_flag',true,[false true],...
  ... % parameters with value that differ from defaults for bold_correct_data
    'outfix','corr',[],...
    'topup_flag',true,[false true],...
    'ref_frame',0,[0 Inf],...
    'skipTRs',16,[0,Inf],...
  ... % parameters for registration and resampling
    'atlasdir',[],[],...
    'atlasname','T1_Atlas/T1_atlas_ABCD',[],...
    'atlasdir_T1T2','T1T2',[],...
    'deoblique_flag',true,[false true],...
  ... % parameters with values that are the same as defaults for bold_correct_data
    'regref_B0uw_flag',true,[false true],...
    'motion_B0uw_flag',false,[false true],...
    'mc_flag',true,[false true],...
    'motion_B0uw_iters',2,[1:10],...
    'min_trans',0.05,[0,1],... % mm
    'min_rot',0.05,[0,1],... % degrees
    'outext','.mgz',{'.mgh','.mgz'},...
    'cleanupflag',true,[false true],... % NOTE: whether bold_correct_data deletes tmp_BOLD_corr;
    'verbose',true,[false true],...     %   cleanup_flag deletes it anyway, with tmpdir
  ... % other hidden parameters
    'fnamestem','BOLD',[],...
    'ext','.mgz',{'.mgh','.mgz'},...
  ... % parameters for comparing field map scans
    'vox2ras_maxdiff',1e-3,[1e-16,1],...
    'offset_maxdiff',20,[1e-16,1e16],...
  ...
    'correct_tags',{'B0unwarp_flag','regref_B0uw_flag',...
                    'motion_B0uw_flag','tshift_flag',...
                    'mc_flag','gradunwarp_flag','fname_out','verbose',...
                    'forceflag',...
                    'topup_flag','epic_flag','epic_inorm_flag','optimize_B0uw_flag',...
                    'fname_B0dx','fname_for','fname_rev',...
                    'fname_B0uw_ref','fname_B0uw_mask','revflag',...
                    'EchoSpacing','AcquisitionColumns','AcquisitionRows',...
                    'PhaseDir','EchoSpacing_for','AcquisitionColumns_for',...
                    'AcquisitionRows_for','EchoSpacing_rev',...
                    'AcquisitionColumns_rev','AcquisitionRows_rev',...
                    'TotalReadoutTime','TotalReadoutTime_for',...
                    'motion_B0uw_iters','min_trans','min_rot',...
                    'fname_motion_B0uw','fname_regref_B0uw',...
                    'tpattern','skipTRs','ts_suffix',...
                    'fname_ref','ref_frame','fname_motion_mat',...
                    'fname_motion_txt','mc_suffix','gruw_type',...
                    'gruw_unwarpflag','gruw_isoctrflag','interpm',...
                    'gruw_jacobian_flag','kernelWidthMax','lambda2',...
                    'kernelWidthMax_vec','lambda2_vec','multi_opt_flag',...
                    'outfix','outext','tmpext','cleanupflag'},[],...
    'regT1_tags',{'fname_T1','T1ContainerPath','T1type','FSContainerPath',...
                  'fname_reg','T2_type','interpm','ext',...
                  'forceflag','atlasdir','atlasname','atlasdir_T1T2','cleanup_flag'},[],...
    'resample_tags',{'fname_out','save_reg_flag','fname_reg',...
                     'native_flag','nvoxels','resolution',...
                     'deoblique_flag','std_orient','M_reg','M_ref',...
                     'volsz_ref','regT1flag','fname_T1','M_T1_to_EPI',...
                     'smooth','rot','trans','EPI_type','ext','forceflag',...
                     'interpm','smf'},[],...
  });
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];

  % sets output directory
  parms.outdir = sprintf('%s/%s/%s/func',parms.output_rootdir,parms.participant_id,parms.session_id);

  % get input file names for func scans
  [parms.fnames_data_func,parms.fnames_info_func] = get_input_data(parms.participant_id,parms.session_id,parms.input_rootdir,'func');

  % get input file names for fmap scans
  [parms.fnames_data_fmap,parms.fnames_info_fmap] = get_input_data(parms.participant_id,parms.session_id,parms.input_rootdir,'fmap');

  % quit if no input files to process
  if isempty(parms.fnames_data_func)
    errmsg = 'not processed: no usable func scans to process';
    errcode = 1;
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % check field map pairs
  [parms,errmsg] = check_fmap_pairs(parms);

  % check whether there are usable field map pairs for B0 unwarp
  if parms.B0unwarp_flag>0 && isempty(parms.fnames_data_fmap)
    errmsg = 'no usable fmap scans necessary for B0 unwarping';
    if parms.B0unwarp_flag>1
      errcode = 1;
      errmsg = ['not processed: ' errmsg];
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    else
      fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    end
    parms.B0unwarp_flag = 0;
  end

  % check whether gradwarp type is missing
  % check all func scans (not just one) before aborting; individual scan
  %   failures are already handled gracefully later in proc_scan/cmig_summarize_errcodes
  gradunwarp_ok = false;
  for s=1:length(parms.fnames_info_func)
    func_scan_info = get_scan_info(parms.fnames_info_func(s));
    if mmil_Philips(func_scan_info)
      gradunwarp_ok = true;
    else
      [gwtype,unwarpflag,isoctrflag] = cmig_get_bids_gradwarpinfo(func_scan_info);
      if ~isnan(gwtype), gradunwarp_ok = true; end
    end
  end
  if ~gradunwarp_ok
    errmsg = 'missing gradwarp type necessary for grad unwarp';
    if parms.gradunwarp_flag>1
      errcode = 1;
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    else
      fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    end
    parms.gradunwarp_flag = 0;
  end

  % check functional scan phase-encode polarities
  parms = check_func_revflags(parms);

  % check T1w scans for registration
  if parms.regT1flag
    % get file names for processed T1w scans
    [parms.fnames_data_t1_proc,parms.fnames_info_t1_proc] = get_input_data(parms.participant_id,parms.session_id,parms.output_rootdir,'T1w',parms.anat_select_flag,1);

    % without a processed T1w (no T1w scan, or cmig_proc_bids_anat failed or
    %   was not run, which is not done from here), process without registration
    if isempty(parms.fnames_data_t1_proc)
      fprintf('%s: WARNING: no processed T1w; func processed without registration to T1\n',mfilename);
      parms.regT1flag = 0;
    end
  end

  % check events files
  parms.fnames_events = cell(size(parms.fnames_data_func));
  for s=1:length(parms.fnames_data_func)
    fname_events = regexprep(parms.fnames_data_func{s},'_bold\.nii\.gz','_events.tsv');
    if exist(fname_events,'file'), parms.fnames_events{s} = fname_events; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;

  % set names of fmap output directory
  parms.outdir_fmap = sprintf('%s/%s/%s/fmap',parms.output_rootdir,parms.participant_id,parms.session_id);

  % set parent tmpdirs so per-scan/per-reference tmpdirs nest under them,
  % enabling selective cleanup (cmig_cleanup_tmpdirs(outdir,tmpdir))
  parms.tmpdir = sprintf('%s/proc_func',parms.outdir);
  parms.tmpdir_fmap = sprintf('%s/proc_fmap',parms.outdir_fmap);

  % check expected output files
  parms.proc_func_flag = check_output_files(parms.fnames_data_func,parms.outdir,parms.forceflag);
  parms.proc_fmap_flag = check_output_files(parms.fnames_data_fmap,parms.outdir_fmap,parms.forceflag);

  % check output events files
  parms.copy_events_flag = check_events_files(parms.fnames_events,parms.outdir,parms.forceflag);

  % check if B0dx output files are missing
  parms.proc_b0dx_flag = check_b0dx_file(parms);

  % set procflag true if any processing needed
  procflag = (parms.proc_fmap_flag || parms.proc_func_flag || parms.copy_events_flag || parms.proc_b0dx_flag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_data,fnames_info] = get_input_data(participant_id,session_id,input_rootdir,modality,select_flag,proc_flag)
  if ~exist('modality','var'), modality = []; end
  if ~exist('select_flag','var') || isempty(select_flag), select_flag = 1; end
  if ~exist('proc_flag','var') || isempty(proc_flag), proc_flag = 0; end
  switch modality
    case 'func'
      series_type = 'bold';
      acq_type = [];
    case 'fmap'
      series_type = 'epi';
      acq_type = 'func';
    case 'T1w'
      modality = 'anat';
      series_type = 'T1w';
      acq_type = [];
    otherwise
      series_type = [];
      acq_type = [];
  end      
  task = [];
  [fnames_data,fnames_info,scan_info] = cmig_get_bids_data(participant_id,session_id,input_rootdir,modality,series_type,task,acq_type);
  if isempty(fnames_data), return; end
  if proc_flag
    % the proc outputs only (desc-preproc), not a rawdata-named file
    %   left in the proc tree
    idx = find(~cellfun(@isempty,regexp(fnames_data,...
      sprintf('_desc-preproc_%s\\.nii',series_type),'once')));
    if ~isempty(scan_info) && length(scan_info)==length(fnames_data)
      scan_info = scan_info(idx);
    end
    fnames_data = fnames_data(idx);
    fnames_info = fnames_info(idx);
    if isempty(fnames_data), return; end
  end
  if select_flag && strcmp(modality,'anat')
    % choose which scan to process based on scans_info
    idx = cmig_select_scan(scan_info);
    fnames_data = fnames_data(idx);
    fnames_info = fnames_info(idx);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = proc_data(parms)
  errcode = 0; errmsg = [];

  % create output directories
  mmil_mkdir(parms.outdir);
  mmil_mkdir(parms.outdir_fmap);

  % copy events files
  % NOTE: copied first, so a failed session still has its events files
  [parms,errcode,errmsg] = copy_events(parms);
  if errcode, return; end

  
  % identify the frame-defining scan, so each sidecar can name its frame
  % NOTE: not in proc_fmap, which returns early when there are no field maps
  parms = set_spatial_reference(parms);

  % process field map scans
  [parms,errcode,errmsg] = proc_fmap(parms);
  if errcode, return; end
  
  % process all functional scans
  [parms,errcode,errmsg] = proc_func(parms);

  % IntendedFor of the field maps, from the bold sidecars, once all are written
  % NOTE: also after a failed bold run, whose sidecar is then absent or old
  set_fmap_intended_for(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_fmap_ref,errcode,errmsg] = proc_fmap_ref(s_fmap,parms)
  fname_fmap_ref = [];
  errcode = 0; errmsg = [];
  s_fmap_ref = select_ref(s_fmap);
  if ~isempty(s_fmap_ref)
    fname_data = parms.fnames_data_fmap{s_fmap_ref};
    fname_info = parms.fnames_info_fmap{s_fmap_ref};
    [fname_fmap_ref,~,errcode,errmsg] = proc_scan(fname_data,fname_info,parms,1,1);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_func_ref,errcode,errmsg] = proc_func_ref(s_func,s_fmap,parms);
  fname_func_ref = [];
  errcode = 0; errmsg = [];
  % process "for" reference functional scan
  s_func_ref = select_func_ref(s_func,s_fmap,parms);
  if ~isempty(s_func_ref)
    fname_data = parms.fnames_data_func{s_func_ref};
    fname_info = parms.fnames_info_func{s_func_ref};
    [fname_func_ref,~,errcode,errmsg] = proc_scan(fname_data,fname_info,parms,1,0);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = reg_ref_to_T1(parms)
  errcode = 0; errmsg = [];
  parms.M_T1_to_BOLD = [];
  if ~parms.regT1flag, return; end
  % set reference fname
  if isempty(parms.fnames_data_fmap)
    if parms.revflag
      parms.fname_ref_regT1 = parms.fname_func_rev_ref;
    else
      parms.fname_ref_regT1 = parms.fname_func_for_ref;
    end
  else
    if parms.revflag
      parms.fname_ref_regT1 = parms.fname_fmap_rev_ref;
    else
      parms.fname_ref_regT1 = parms.fname_fmap_for_ref;
    end
  end
  % import processed T1
  tmpdir = cmig_get_fileparts(parms.fname_ref_regT1);
  fname_T1 = parms.fnames_data_t1_proc{1};
  parms.fname_T1 = cmig_import_to_mgz(fname_T1,tmpdir,[],parms.forceflag);
  % import brain mask as {T1stem}_mask so MMIL_Register_T2_to_T1 finds it
  % NOTE: the proc T1w is _desc-preproc_T1w; its mask carries no desc-preproc
  fname_T1_mask = regexprep(fname_T1,'(_desc-preproc)?_([A-Za-z0-9]+)\.nii\.gz$','_space-T1w_desc-brain_mask.nii.gz');
  if exist(fname_T1_mask,'file')
    [~,fstem_T1] = cmig_get_fileparts(parms.fname_T1);
    fname_T1_mask = cmig_import_to_mgz(fname_T1_mask,tmpdir,...
      sprintf('%s_mask',fstem_T1),parms.forceflag);
  end
  % register reference to T1
  tparms = parms;
  tparms.T2_type = ['T2_' parms.fnamestem];
  args = mmil_parms2args(tparms,parms.regT1_tags);
  [parms.M_T1_to_BOLD,parms.RegInfo,errcode] = ...
    MMIL_Register_T2_to_T1(parms.fname_ref_regT1,args{:});
  if errcode, errmsg = 'registration to T1 failed'; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = reg_ref_func_to_fmap(parms)
  errcode = 0; errmsg = [];
  parms.M_GE_to_SE = [];
  % set file names of volumes to be registered
  if parms.revflag
    fname_GE = parms.fname_func_rev_ref;
    fname_SE = parms.fname_fmap_rev_ref;
  else
    fname_GE = parms.fname_func_for_ref;
    fname_SE = parms.fname_fmap_for_ref;
  end
  if isempty(fname_SE)
    errcode = 1;
    errmsg = 'registration of func to fmap failed because SE reference is missing';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  if isempty(fname_GE)
    errcode = 1;
    errmsg = 'registration of func to fmap failed because GE reference is missing';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  [~,fstem_SE] = fileparts(fname_SE);
  [fpath,fstem_GE] = fileparts(fname_GE);
  if ~isempty(regexp(fname_GE,fname_SE,'ONCE'))
    if parms.verbose
      fprintf('%s: skipping registration between %s and %s...\n',mfilename,fstem_SE,fstem_GE);
    end
    return;
  end
  if parms.verbose
    fprintf('%s: registering reference images %s and %s...\n',...
      mfilename,fstem_SE,fstem_GE);
  end
  fname_regref_B0uw = sprintf('%s/%s_regref_B0uw.mat',fpath,fstem_GE);
  if ~exist(fname_regref_B0uw,'file') || parms.forceflag
    tparms = [];
    tparms.outdir = sprintf('%s/tmp_BOLD_corr/reg_GE_%s_SE_%s',fpath,fstem_GE,fstem_SE);
    tparms.cleanup_flag = 0;
    tparms.forceflag = parms.forceflag;
    args = mmil_parms2args(tparms);
    M_GE_to_SE = mmil_reg_SEtoGE(fname_SE,fname_GE,args{:});
    save(fname_regref_B0uw,'M_GE_to_SE','fname_SE','fname_GE');
  else
    tmp = load(fname_regref_B0uw);
    M_GE_to_SE = tmp.M_GE_to_SE;
  end
  parms.M_GE_to_SE = M_GE_to_SE;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname,errcode,errmsg] = resamp_ref(fname,save_reg_flag,regref_flag,parms)
  errcode = 0; errmsg = [];
  if isempty(fname), return; end
  % whether to save updated registration matrix for resampled volume
  if ~exist('save_reg_flag','var') || isempty(save_reg_flag), save_reg_flag = 0; end
  % whether to use GE (gradient echo) to SE (spin echo) registration
  if ~exist('regref_flag','var') || isempty(regref_flag), regref_flag = 0; end
  [fpath,fstem,fext] = fileparts(fname);
  % set fname_out
  outfix = 'resBOLD';
  parms.fname_out = sprintf('%s/%s_%s%s',fpath,fstem,outfix,parms.ext);
  % set fname_reg
  if parms.regT1flag>0
    parms.fname_reg = sprintf('%s/%s_%s_regT1.mat',fpath,fstem,outfix);
  else
    parms.fname_reg = [];
  end
  % set additional parameters
  parms.save_reg_flag = save_reg_flag;
  parms.regT1flag = (parms.regT1flag==2);
  parms.M_T1_to_EPI = parms.M_T1_to_BOLD;
  parms.EPI_type = parms.fnamestem;
  if regref_flag
    parms.M_reg = inv(parms.M_GE_to_SE);
  else
    parms.M_reg = [];
  end
  % resample volume to final resolution
  args = mmil_parms2args(parms,parms.resample_tags);
  try
    epi_resample_data(fname,args{:});
  catch me
    errcode = 1;
    errmsg = sprintf('resampling to final resolution failed for %s:\n%s',...
      fname,me.message);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  fname = parms.fname_out;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = proc_fmap(parms)
  errcode = 0; errmsg = [];
  
  % skip if no field map scans
  if isempty(parms.fnames_data_fmap), return; end

  % skip if output files exist and B0dx files are present
  if ~parms.proc_fmap_flag && ~parms.proc_func_flag && ~parms.proc_b0dx_flag, return; end

  % process reference fmap scans
  [parms.fname_fmap_for_ref,errcode,errmsg] = proc_fmap_ref(parms.s_fmap_for,parms);
  if errcode, return; end
  [parms.fname_fmap_rev_ref,errcode,errmsg] = proc_fmap_ref(parms.s_fmap_rev,parms);
  if errcode, return; end

  % register reference to T1
  [parms,errcode,errmsg] = reg_ref_to_T1(parms);
  if errcode, return; end
  save_reg_flag_for = ~parms.revflag;
  save_reg_flag_rev = parms.revflag;

  % resample reference fmap scans
  parms.fname_fmap_for_ref = resamp_ref(parms.fname_fmap_for_ref,save_reg_flag_for,0,parms);
  parms.fname_fmap_rev_ref = resamp_ref(parms.fname_fmap_rev_ref,save_reg_flag_rev,0,parms);

  % process fmap scans
  nscans = length(parms.fnames_data_fmap);
  errcodes = zeros(nscans,1);
  errmsgs = cell(nscans,1);
  for s=1:nscans
    proc_flag = check_output_files(parms.fnames_data_fmap(s),parms.outdir_fmap,parms.forceflag) || parms.proc_b0dx_flag;
    if proc_flag
      [parms.fname_data_fmap_proc{s},parms.fname_info_fmap_proc{s},errcodes(s),errmsgs{s}] = ...
        proc_scan(parms.fnames_data_fmap{s},parms.fnames_info_fmap{s},parms,0,1);
    end
  end
  [errcode,errmsg] = cmig_summarize_errcodes(errcodes,errmsgs);
  if errcode, return; end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = proc_func(parms)
  errcode = 0; errmsg = [];
  
  % skip if output files exist
  if ~parms.proc_func_flag, return; end

  % initialize M_GE_to_SE; set in reg_ref_func_to_fmap only when fmaps are available
  parms.M_GE_to_SE = [];

  % process reference functional scans
  [parms.fname_func_for_ref,errcode,errmsg] = proc_func_ref(parms.s_func_for,parms.s_fmap_for,parms);
  if errcode, return; end
  [parms.fname_func_rev_ref,errcode,errmsg] = proc_func_ref(parms.s_func_rev,parms.s_fmap_rev,parms);
  if errcode, return; end

  if ~isempty(parms.fnames_data_fmap)
    % register reference func scans to reference fmap scans
    [parms,errcode,errmsg] = reg_ref_func_to_fmap(parms);
    if errcode, return; end;
    save_reg_flag_for = 0;
    save_reg_flag_rev = 0;
    % tie the func reference grid to the resampled fmap reference grid
    %   so func and fmap outputs are voxel-aligned by construction
    if parms.revflag
      fname_fmap_ref = parms.fname_fmap_rev_ref;
    else
      fname_fmap_ref = parms.fname_fmap_for_ref;
    end
    if ~isempty(fname_fmap_ref)
      fpath_fmap_ref = cmig_get_fileparts(fname_fmap_ref);
      [parms.M_ref,parms.volsz_ref] = mmil_load_mgh_info(fname_fmap_ref,[],fpath_fmap_ref);
    end
  else
    % register reference to T1
    [parms,errcode,errmsg] = reg_ref_to_T1(parms);
    if errcode, return; end;
    save_reg_flag_for = ~parms.revflag;
    save_reg_flag_rev = parms.revflag;
  end
  
  % resample reference func scans
  parms.fname_func_for_ref = resamp_ref(parms.fname_func_for_ref,save_reg_flag_for,parms.regref_B0uw_flag,parms);
  parms.fname_func_rev_ref = resamp_ref(parms.fname_func_rev_ref,save_reg_flag_rev,parms.regref_B0uw_flag,parms);

  % process func scans
  nscans = length(parms.fnames_data_func);
  errcodes = zeros(nscans,1);
  errmsgs = cell(nscans,1);
  for s=1:nscans
    proc_flag = check_output_files(parms.fnames_data_func(s),parms.outdir,parms.forceflag);
    if proc_flag
      [parms.fnames_data_func_proc{s},parms.fnames_info_func_proc{s},errcodes(s),errmsgs{s}] = ...
        proc_scan(parms.fnames_data_func{s},parms.fnames_info_func{s},parms,0,0);
    end
  end
  [errcode,errmsg] = cmig_summarize_errcodes(errcodes,errmsgs);
  if errcode, return; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = copy_events(parms)
  errcode = 0; errmsg = [];

  % skip if output files exist
  if ~parms.copy_events_flag, return; end
  % copy events file if found (and output does not exist)
  nfiles = length(parms.fnames_events);
  errcodes = zeros(nfiles,1);
  errmsgs = cell(nfiles,1);
  for s=1:nfiles
    copy_flag = check_events_files(parms.fnames_events(s),parms.outdir,parms.forceflag);
    if copy_flag
      [parms.fnames_events_proc{s},errcodes(s),errmsgs{s}] = copy_events_file(parms.fnames_events{s},parms);
    end
  end
  [errcode,errmsg] = cmig_summarize_errcodes(errcodes,errmsgs);
  if errcode, return; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errmsg] = check_fmap_pairs(parms)
  errmsg = [];
  ees_reltol = 1e-4; % relative tolerance for the echo spacing match
  parms.fmap_info = [];
  parms.s_fmap_for = [];
  parms.s_fmap_rev = [];
  if isempty(parms.fnames_data_fmap), return; end  
  
  % compile info about field maps
  parms.fmap_info = get_scan_info(parms.fnames_info_fmap);

  % identify pairs of forward and reverse field map scans
  fmap_revflags = get_revflags(parms.fmap_info);
  fmap_labels = {parms.fmap_info.B0FieldIdentifier};
  [uniq_fmap_labels,~,idx_all] = unique(fmap_labels,'stable');
  for i=1:length(uniq_fmap_labels)
    % NOTE: idx_all holds each scan's group number, not unique's first-occurrence index
    idx = find(idx_all==i);
    r = find(fmap_revflags(idx));
    f = find(~fmap_revflags(idx));
    if ~isempty(r) && ~isempty(f)
      % NOTE: a pair is one forward and one reverse scan; anything else is skipped,
      %   including unassigned field maps, which share the default label
      if numel(f)>1 || numel(r)>1
        errmsg = sprintf('%d forward and %d reverse field maps for B0FieldIdentifier %s; expected one of each',...
          numel(f),numel(r),uniq_fmap_labels{i});
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      s_for = idx(f);
      s_rev = idx(r);

      % check phase-encode directions
      phasedirs = get_phasedirs(parms.fmap_info([s_for,s_rev]));
      if length(unique(phasedirs)) > 1
        errmsg = sprintf('phase-encode direction differs for %s and %s (%d vs %d)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          parms.fmap_info(s_for).PhaseEncodingDirection,parms.fmap_info(s_rev).PhaseEncodingDirection);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      % check voxel size mismatch
      if ~isnan(parms.fmap_info(s_for).SliceThickness) && ...
         ~isnan(parms.fmap_info(s_rev).SliceThickness) && ...
         parms.fmap_info(s_for).SliceThickness ~= parms.fmap_info(s_rev).SliceThickness
        errmsg = sprintf('SliceThickness differs for %s and %s (%g vs %g)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          parms.fmap_info(s_for).SliceThickness,parms.fmap_info(s_rev).SliceThickness);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      if ~isnan(parms.fmap_info(s_for).ReconMatrixPE) && ...
         ~isnan(parms.fmap_info(s_rev).ReconMatrixPE) && ...
         parms.fmap_info(s_for).ReconMatrixPE ~= parms.fmap_info(s_rev).ReconMatrixPE
        errmsg = sprintf('ReconMatrixPE differs for %s and %s (%d vs %d)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          parms.fmap_info(s_for).ReconMatrixPE,parms.fmap_info(s_rev).ReconMatrixPE);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      if ~isnan(parms.fmap_info(s_for).FrequencyEncodingSteps) && ...
         ~isnan(parms.fmap_info(s_rev).FrequencyEncodingSteps) && ...
         parms.fmap_info(s_for).FrequencyEncodingSteps ~= parms.fmap_info(s_rev).FrequencyEncodingSteps
        errmsg = sprintf('FrequencyEncodingSteps differs for %s and %s (%d vs %d)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          parms.fmap_info(s_for).FrequencyEncodingSteps,parms.fmap_info(s_rev).FrequencyEncodingSteps);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      if ~isnan(parms.fmap_info(s_for).PhaseEncodingSteps) && ...
         ~isnan(parms.fmap_info(s_rev).PhaseEncodingSteps) && ...
         parms.fmap_info(s_for).PhaseEncodingSteps ~= parms.fmap_info(s_rev).PhaseEncodingSteps
        errmsg = sprintf('PhaseEncodingSteps differs for %s and %s (%d vs %d)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          parms.fmap_info(s_for).PhaseEncodingSteps,parms.fmap_info(s_rev).PhaseEncodingSteps);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      % check echo spacing mismatch, within a relative tolerance: dcm2niix
      %   derives Philips echo spacing from the water-fat shift, so the two
      %   directions can differ by a few ppm
      ees_for = parms.fmap_info(s_for).EffectiveEchoSpacing;
      ees_rev = parms.fmap_info(s_rev).EffectiveEchoSpacing;
      if ~isnan(ees_for) && ~isnan(ees_rev) && ...
         abs(ees_for - ees_rev) > ees_reltol * max(abs(ees_for),abs(ees_rev))
        errmsg = sprintf('EffectiveEchoSpacing differs for %s and %s (%g vs %g)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          parms.fmap_info(s_for).EffectiveEchoSpacing,parms.fmap_info(s_rev).EffectiveEchoSpacing);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end      
      % check slice orientation mismatch
      if any(parms.fmap_info(s_for).ImageOrientationPatientDICOM ~= parms.fmap_info(s_rev).ImageOrientationPatientDICOM & ...
             ~(isnan(parms.fmap_info(s_for).ImageOrientationPatientDICOM) & isnan(parms.fmap_info(s_rev).ImageOrientationPatientDICOM)))
        errmsg = sprintf('ImageOrientationPatientDICOM differs for %s and %s (%s vs %s)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          strtrim(sprintf('%d ',parms.fmap_info(s_for).ImageOrientationPatientDICOM(:))),...
          strtrim(sprintf('%d ',parms.fmap_info(s_rev).ImageOrientationPatientDICOM(:))));
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end      
      % check shim setting mismatch
      if any(parms.fmap_info(s_for).ShimSetting ~= parms.fmap_info(s_rev).ShimSetting & ...
             ~(isnan(parms.fmap_info(s_for).ShimSetting) & isnan(parms.fmap_info(s_rev).ShimSetting)))
        errmsg = sprintf('ShimSetting differs for %s and %s (%s vs %s)',...
          parms.fnames_info_fmap{s_for},parms.fnames_info_fmap{s_rev},...
          strtrim(sprintf('%d ',parms.fmap_info(s_for).ShimSetting(:))),...
          strtrim(sprintf('%d ',parms.fmap_info(s_rev).ShimSetting(:))));
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      % check vox2ras/dimension mismatch between the actual images; JSON
      % sidecars can match while FOV/cropping/offset differ (e.g. when a
      % scan operator diverges from protocol). Read headers only (no file
      % is created; this runs on every invocation, before check_output).
      [~,M_for,volsz_for] = cmig_read_nifti(parms.fnames_data_fmap{s_for},1);
      [~,M_rev,volsz_rev] = cmig_read_nifti(parms.fnames_data_fmap{s_rev},1);
      % tight rotation tolerance: any rotation changes the distortion
      % pattern (topup cannot correct it); loose offset tolerance: topup's
      % motion-correction step absorbs a translational offset
      [~,rotmatch,offmatch,maxrot,maxoff] = ...
        fs_vox2ras_match(M_for,M_rev,parms.vox2ras_maxdiff,parms.offset_maxdiff);
      if ~rotmatch || ~offmatch
        errmsg = sprintf('vox2ras mismatch for %s and %s (max rotation %0.4f deg, max offset %0.4f mm)',...
          parms.fnames_data_fmap{s_for},parms.fnames_data_fmap{s_rev},maxrot,maxoff);
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      if any(volsz_for(1:3) ~= volsz_rev(1:3))
        errmsg = sprintf('volume dimensions differ for %s and %s ([%s] vs [%s])',...
          parms.fnames_data_fmap{s_for},parms.fnames_data_fmap{s_rev},...
          strtrim(sprintf('%d ',volsz_for(1:3))),strtrim(sprintf('%d ',volsz_rev(1:3))));
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
        continue;
      end
      parms.s_fmap_for(end+1) = s_for;
      parms.s_fmap_rev(end+1) = s_rev;
    end
  end

  % exclude unpaired fmap scans
  s_fmap = union(parms.s_fmap_for,parms.s_fmap_rev);
  parms.fnames_data_fmap = parms.fnames_data_fmap(s_fmap);
  parms.fnames_info_fmap = parms.fnames_info_fmap(s_fmap);
  parms.fmap_info = parms.fmap_info(s_fmap);
  [~,parms.s_fmap_for] = intersect(s_fmap,parms.s_fmap_for,'stable');
  [~,parms.s_fmap_rev] = intersect(s_fmap,parms.s_fmap_rev,'stable');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_info = get_scan_info(fnames_info)
  scan_info = [];
  for s=1:length(fnames_info)
    tmp_info = load_scan_info(fnames_info{s});
    [scan_info,tmp_info] = reconcile_fields(scan_info,tmp_info);
    scan_info = cat(2,scan_info,tmp_info);    
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_info = load_scan_info(fname_info)
  scan_info = [];
  try
    scan_info = jsondecode(fileread(fname_info));
  catch me
    error('failed to load json file %s:\n%s',fname_info,me.message);
  end
  % ensure key variables have values
  scan_info.AcquisitionDateTime = mmil_getfield(scan_info,'AcquisitionDateTime','UNKNOWN');
  scan_info.AcquisitionMatrixPE = mmil_getfield(scan_info,'AcquisitionMatrixPE',NaN);
  if isfield(scan_info,'B0FieldIdentifier') || isfield(scan_info,'IntendedFor')
    scan_info.B0FieldIdentifier = mmil_getfield(scan_info,'B0FieldIdentifier','UNKNOWN');
  elseif isfield(scan_info,'B0FieldSource')
    scan_info.B0FieldSource = mmil_getfield(scan_info,'B0FieldSource','UNKNOWN');
  end
  % FrequencyEncodingSteps: freq encoding k-space steps; set by manufacturer
  % Siemens: BaseResolution field; Philips: FrequencyEncodingSteps field; GE: AcquisitionMatrixPE
  if mmil_Siemens(scan_info)
    scan_info.FrequencyEncodingSteps = mmil_getfield(scan_info,'BaseResolution',NaN);
  elseif mmil_Philips(scan_info)
    scan_info.FrequencyEncodingSteps = mmil_getfield(scan_info,'FrequencyEncodingSteps',NaN);
  else
    scan_info.FrequencyEncodingSteps = mmil_getfield(scan_info,'FrequencyEncodingSteps',scan_info.AcquisitionMatrixPE);
  end
  scan_info.EchoTrainLength = mmil_getfield(scan_info,'EchoTrainLength',NaN);
  scan_info.EffectiveEchoSpacing = mmil_getfield(scan_info,'EffectiveEchoSpacing',NaN);
  scan_info.ImageOrientationPatientDICOM = mmil_getfield(scan_info,'ImageOrientationPatientDICOM',NaN);
  scan_info.InPlanePhaseEncodingDirectionDICOM = mmil_getfield(scan_info,'InPlanePhaseEncodingDirectionDICOM',NaN);
  scan_info.PhaseEncodingDirection = mmil_getfield(scan_info,'PhaseEncodingDirection','UNKNOWN');
  scan_info.PhaseEncodingSteps = mmil_getfield(scan_info,'PhaseEncodingSteps',NaN);
  scan_info.ReconMatrixPE = mmil_getfield(scan_info,'ReconMatrixPE',NaN);
  scan_info.RepetitionTime = mmil_getfield(scan_info,'RepetitionTime',NaN);
  scan_info.ShimSetting = mmil_getfield(scan_info,'ShimSetting',NaN);
  scan_info.SliceThickness = mmil_getfield(scan_info,'SliceThickness',NaN);
  scan_info.SliceTiming = mmil_getfield(scan_info,'SliceTiming',NaN);
  scan_info.SpacingBetweenSlices = mmil_getfield(scan_info,'SpacingBetweenSlices',NaN);
  scan_info.TotalReadoutTime = mmil_getfield(scan_info,'TotalReadoutTime',NaN);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readout_time = get_readout_time(scan_info)
  % TotalReadoutTime (s) from the scan's sidecar, empty where it has none
  readout_time = [];
  if ~isnan(scan_info.TotalReadoutTime) && scan_info.TotalReadoutTime > 0
    readout_time = double(scan_info.TotalReadoutTime);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revflags = get_revflags(scan_info)
  dirs = {scan_info.PhaseEncodingDirection};
  revflags = contains(dirs,'-');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phasedirs = get_phasedirs(scan_info)
  phasedirs = {scan_info.PhaseEncodingDirection};
  phasedirs = cellfun(@(x) x(1),phasedirs,'UniformOutput',false);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_func_revflags(parms)
  parms.func_info = [];
  parms.s_func_for = [];
  parms.s_func_rev = [];
  parms.revflag = 0;
  % compile info about functional scans
  parms.func_info = get_scan_info(parms.fnames_info_func);
  % identify scans with forward and reverse polarity
  func_revflags = get_revflags(parms.func_info);
  parms.s_func_for = find(~func_revflags);
  parms.s_func_rev = find(func_revflags);
  % set revflag based on number of scans
  if length(parms.s_func_rev) > length(parms.s_func_for), parms.revflag = 1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s_ref = select_ref(s)
  s_ref = [];
  if isempty(s); return; end
  s_ref = s(floor(median(1:length(s))));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s_ref = select_func_ref(s_func,s_fmap,parms)
  % select the session reference func scan for one polarity
  % NOTE: each fmap pair gets the func scan acquired next after it (else nearest
  %   before), closest in head position; the median pair's func scan is used
  s_ref = [];
  if isempty(s_func), return; end
  if ~isempty(s_fmap)
    t_func = {parms.func_info(s_func).AcquisitionDateTime};
    t_fmap = {parms.fmap_info(s_fmap).AcquisitionDateTime};
    if ~any(strcmp('UNKNOWN',cat(2,t_func,t_fmap)))
      % ISO 8601 date-times sort lexicographically, so unique ranks
      %   reproduce acquisition order
      [~,~,ranks] = unique(cat(2,t_func,t_fmap));
      r_func = ranks(1:length(t_func));
      r_fmap = ranks(length(t_func)+1:end);
      npairs = length(s_fmap);
      s_refs = zeros(npairs,1);
      for i=1:npairs
        dr = r_func - r_fmap(i);
        dr(dr<0) = Inf; % prefer func scans after the fmap
        [mindiff,ind] = min(dr);
        if isinf(mindiff) % no func scans after, so look before
          [~,ind] = min(abs(r_func - r_fmap(i)));
        end
        s_refs(i) = s_func(ind);
      end
      % use the func scan associated with the median fmap pair
      s_ref = s_refs(floor(median(1:npairs)));
      return;
    end
  end
  % without fmap pairs or usable acquisition times, use median func scan
  s_ref = select_ref(s_func);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_spatial_reference(parms)
  % identify the scan that defines the session frame and the bids:: URI that
  %   every other scan in that frame points at
  % NOTE: the reference field map defines the frame, since the reference func
  %   scan is resampled into it; with no field maps, the reference func scan
  parms.fname_data_ref = [];
  parms.spatial_reference = [];
  if ~isempty(parms.fnames_data_fmap)
    s_ref = select_ref(parms.s_fmap_for);
    if isempty(s_ref), return; end
    parms.fname_data_ref = parms.fnames_data_fmap{s_ref};
    subdir = 'fmap';
  else
    s_ref = select_ref(parms.s_func_for);
    if isempty(s_ref), return; end
    parms.fname_data_ref = parms.fnames_data_func{s_ref};
    subdir = 'func';
  end
  % the exported derivative, by its proc name
  [~,fstem] = cmig_get_fileparts(parms.fname_data_ref);
  parms.spatial_reference = proc_uri(parms,subdir,fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_for,fname_rev,for_info,rev_info,fallback_flag] = select_fmap_pair(fname_info,parms)
  fname_for = []; for_info = [];
  fname_rev = []; rev_info = [];
  fallback_flag = 0;
  if isempty(parms.fnames_data_fmap), return; end
  [fname_for,fname_rev,fallback_flag] = cmig_select_bids_fmap(fname_info,parms.fnames_data_fmap,parms.fnames_info_fmap);
  if isempty(fname_for) || isempty(fname_rev), return; end
  s_for = find(strcmp(parms.fnames_data_fmap,fname_for));
  s_rev = find(strcmp(parms.fnames_data_fmap,fname_rev));
  if isempty(s_for) || isempty(s_rev), return; end
  for_info = parms.fmap_info(s_for);
  rev_info = parms.fmap_info(s_rev);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function procflag = check_output_files(fnames_data,outdir,forceflag)
  procflag = 0;
  if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end
  if forceflag, procflag = 1; return; end
  % set procflag=1 if any output file does not exist
  for s=1:length(fnames_data)
    [~,fstem] = cmig_get_fileparts(fnames_data{s});
    fstem = set_proc_fstem(fstem);
    fname_out = sprintf('%s/%s.nii.gz',outdir,fstem);
    if ~exist(fname_out,'file'), procflag = 1; return; end
    fname_out_json = sprintf('%s/%s.json',outdir,fstem);
    if ~exist(fname_out_json,'file'), procflag = 1; return; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function procflag = check_events_files(fnames_events,outdir,forceflag)
  procflag = 0;
  if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end
  if forceflag, procflag = 1; return; end
  % set procflag=1 if any output file does not exist, or is older than the
  %   rawdata events file, which is regenerated when its converter changes;
  %   cmig_analyze_bids_func builds its regressors from these copies
  for s=1:length(fnames_events)
    if isempty(fnames_events{s}), continue; end
    [~,fstem,fext] = cmig_get_fileparts(fnames_events{s});
    fname_out = sprintf('%s/%s%s',outdir,fstem,fext);
    if ~exist(fname_out,'file'), procflag = 1; return; end
    d_in = dir(fnames_events{s});
    d_out = dir(fname_out);
    if ~isempty(d_in) && ~isempty(d_out) && d_in(1).datenum > d_out(1).datenum
      procflag = 1; return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function procflag = check_b0dx_file(parms)
  procflag = 0;
  if parms.B0unwarp_flag == 0 || isempty(parms.fnames_data_fmap), return; end
  if parms.forceflag, procflag = 1; return; end
  if isempty(parms.s_fmap_for) || isempty(parms.s_fmap_rev), return; end
  for i=1:length(parms.s_fmap_for)
    [~,fstem_for] = cmig_get_fileparts(parms.fnames_data_fmap{parms.s_fmap_for(i)});
    [~,fstem_rev] = cmig_get_fileparts(parms.fnames_data_fmap{parms.s_fmap_rev(i)});
    fname_B0dx = set_b0dx_fname(fstem_for,fstem_rev,parms);
    fname_B0dx_json = regexprep(fname_B0dx,'\.nii\.gz$','.json');
    if ~exist(fname_B0dx,'file') || ~exist(fname_B0dx_json,'file')
      procflag = 1;
      return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_data_proc,fname_info_proc,errcode,errmsg] = proc_scan(fname_data,fname_info,parms,ref_flag,fmap_flag)
  fname_data_proc = []; fname_info_proc = [];
  errcode = 0; errmsg = [];
  if ~mmil_check_nargs(nargin,3), return; end
  if ~exist('ref_flag','var') || isempty(ref_flag), ref_flag = 0; end
  if ~exist('fmap_flag','var'), fmap_flag = []; end

  % load info file
  scan_info = load_scan_info(fname_info);

  % determine which type of scan this is
  if isempty(fmap_flag)
    if isfield(scan_info,'B0FieldIdentifier')
      fmap_flag = 1;
    else
      fmap_flag = 0;  
    end
  end

  % specify output directory depending on scan type
  if fmap_flag
    outdir = parms.outdir_fmap;
    parent_tmpdir = parms.tmpdir_fmap;
  else
    outdir = parms.outdir;
    parent_tmpdir = parms.tmpdir;
  end
  
  % create temporary parameters for bold_correct_data
  tparms = parms;
  
  % set revflag
  tparms.revflag = get_revflags(scan_info);
  
  % select fmap pair to use for B0 unwarp
  if tparms.B0unwarp_flag
    [fname_for,fname_rev,for_info,rev_info,fmap_fallback_flag] = select_fmap_pair(fname_info,parms);
    % skip if no valid fmap pair
    if isempty(fname_for) || isempty(fname_rev)
      errmsg = sprintf('skipping B0 unwarp due to missing field map scans');
      if tparms.B0unwarp_flag > 1
        errcode = 1;
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      else
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
      end
      tparms.B0unwarp_flag = 0;
    end
  end

  % get gradwarp info
  if tparms.gradunwarp_flag
    [tparms.gruw_type,tparms.gruw_unwarpflag,tparms.gruw_isoctrflag] = cmig_get_bids_gradwarpinfo(scan_info);
    if isnan(tparms.gruw_type)
      errmsg = sprintf('skipping grad unwarp due to missing grad warp type');
      if tparms.gradunwarp_flag > 1 && ~mmil_Philips(scan_info)
        errcode = 1;
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      else
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
      end
      tparms.gradunwarp_flag = 0;
      % NOTE: NaN is outside bold_correct_data's [0:18] range for gruw_type
      tparms.gruw_type = [];
    end
  end

  % create tmp dir nested under parent_tmpdir
  [~,fstem] = cmig_get_fileparts(fname_data);
  tmpdir = sprintf('%s/%s',parent_tmpdir,fstem);
  mmil_mkdir(tmpdir);

  % import data (convert to mgz)
  fname_tmp = cmig_import_to_mgz(fname_data,tmpdir,[],parms.forceflag);

  % import field maps
  if tparms.B0unwarp_flag
    [~,fstem_for] = cmig_get_fileparts(fname_for);
    tmpdir_for = sprintf('%s/%s',parms.tmpdir_fmap,fstem_for);
    mmil_mkdir(tmpdir_for);
    fname_for_tmp = cmig_import_to_mgz(fname_for,tmpdir_for,[],parms.forceflag);
    [~,fstem_rev] = cmig_get_fileparts(fname_rev);
    tmpdir_rev = sprintf('%s/%s',parms.tmpdir_fmap,fstem_rev);
    mmil_mkdir(tmpdir_rev);
    fname_rev_tmp = cmig_import_to_mgz(fname_rev,tmpdir_rev,[],parms.forceflag);
    % set file name of B0dx volume
    tparms.fname_B0dx = sprintf('%s/%sVS%s_B0dx%s',tmpdir_for,fstem_for,fstem_rev,parms.outext);
    % registration of the bold to its field map pair, named here to be read back
    if ~ref_flag && ~fmap_flag && tparms.regref_B0uw_flag
      tparms.fname_regref_B0uw = sprintf('%s/%s_regref_B0uw.mat',tmpdir,fstem);
    end
  else
    fname_for_tmp = [];
    fname_rev_tmp = [];
  end
  tparms.fname_for = fname_for_tmp;
  tparms.fname_rev = fname_rev_tmp;

  % get volume info for input volume
  [M_raw,volsz_raw] = mmil_load_mgh_info(fname_tmp,[],tmpdir);

  % save first frame of reference scan if multi-frame
  if ref_flag && volsz_raw(4)>1
    fname_tmp = copy_first_frame(fname_tmp,parms); 
  end

  if ref_flag
    % no slice timing or motion correction for single-frame reference
    tparms.tshift_flag = 0;
    tparms.mc_flag = 0;
    % self as reference
    tparms.fname_B0uw_ref = fname_tmp;
    tparms.motion_B0uw_iters = 1;
  else
    % set file names for reference scans
    if fmap_flag
      % set reference for motion correction of fmap scans within session
      % NOTE: this should be a processed image in output resolution
      if tparms.revflag
        tparms.fname_ref = parms.fname_fmap_rev_ref;
      else
        tparms.fname_ref = parms.fname_fmap_for_ref;
      end
      tparms.fname_B0uw_ref = [];
    else
      % set reference for motion correction of func scans within session
      % NOTE: this should be a processed image in output resolution
      if tparms.revflag
        tparms.fname_ref = parms.fname_func_rev_ref;
      else
        tparms.fname_ref = parms.fname_func_for_ref;
      end
      % no B0 unwarp reference needed without fmaps
      tparms.fname_B0uw_ref = [];
    end    
  end

  % disable registration to self for SE field map
  if fmap_flag
    tparms.regref_B0uw_flag = 0;
  end

  % set outfix
  % NOTE: distinct for reference passes, so they never share a name with the
  %   scan pass of the same fmap scan; avoids a doubled _resBOLD suffix
  if ref_flag
    tparms.outfix = 'corr_ref';
  else
    tparms.outfix = 'corr_resBOLD';
  end

  % set additional parameters related to B0 unwarping
  if tparms.B0unwarp_flag
    % set echo spacing
    tparms.EchoSpacing = single(scan_info.EffectiveEchoSpacing);  
    tparms.EchoSpacing_for = single(for_info.EffectiveEchoSpacing);  
    tparms.EchoSpacing_rev = single(rev_info.EffectiveEchoSpacing);  
    % set total readout time, for topup and applytopup
    % NOTE: the field map pair's for topup, the scan's own for applytopup;
    %   without the pair's, both are given 1 s
    tparms.TotalReadoutTime_for = get_readout_time(for_info);
    tparms.TotalReadoutTime = get_readout_time(scan_info);
    if isempty(tparms.TotalReadoutTime_for)
      fprintf('%s: WARNING: no TotalReadoutTime for %s, so topup is given 1 s\n',...
        mfilename,fstem_for);
    end
    % set phase-encoding direction and nrows and ncols
    phasedir = get_phasedirs(scan_info);
    if strcmp(phasedir,'j')
      tparms.PhaseDir = 'COL';
      % set number of rows
      tparms.AcquisitionRows = single(scan_info.FrequencyEncodingSteps);
      tparms.AcquisitionRows_for = single(for_info.FrequencyEncodingSteps);
      tparms.AcquisitionRows_rev = single(rev_info.FrequencyEncodingSteps);
      % set number of columns
      tparms.AcquisitionColumns = single(scan_info.AcquisitionMatrixPE);
      tparms.AcquisitionColumns_for = single(for_info.AcquisitionMatrixPE);
      tparms.AcquisitionColumns_rev = single(rev_info.AcquisitionMatrixPE);
    elseif strcmp(phasedir,'i')
      tparms.PhaseDir = 'ROW';
      % set number of rows
      tparms.AcquisitionRows = single(scan_info.AcquisitionMatrixPE);
      tparms.AcquisitionRows_for = single(for_info.AcquisitionMatrixPE);
      tparms.AcquisitionRows_rev = single(rev_info.AcquisitionMatrixPE);
      % set number of columns
      tparms.AcquisitionColumns = single(scan_info.FrequencyEncodingSteps);
      tparms.AcquisitionColumns_for = single(for_info.FrequencyEncodingSteps);
      tparms.AcquisitionColumns_rev = single(rev_info.FrequencyEncodingSteps);
    end
  end

  % call bold_correct_data
  % NOTE: correct_data reports its own errors
  [fname_tmp,fname_B0dx,fname_motion_txt,start_info,errcode,errmsg] = correct_data(fname_tmp,tparms);
  if errcode, return; end

  % export B0dx displacement field to outdir_fmap with BIDS-style name
  % NOTE: on the forward field map scan's pass only; reference passes precede
  %   the output frame, and export_b0dx needs the forward scan's motion matrix
  if ~ref_flag && ~tparms.revflag && tparms.B0unwarp_flag && ...
     ~isempty(fname_B0dx) && exist(fname_B0dx,'file')
    export_b0dx(fname_B0dx,fstem_for,fstem_rev,tmpdir_for,for_info,tparms.outfix,parms);
  end

  % stop here for reference scans
  if ref_flag
    fname_data_proc = fname_tmp;
    return;
  end

  % get volume info for output volume
  [M_proc,volsz_proc] = mmil_load_mgh_info(fname_tmp,[],tmpdir);

  % save output data file
  fname_data_proc = export_data(fname_tmp,set_proc_fstem(fstem),outdir,scan_info.RepetitionTime);

  % write motion tsv
  write_motion_tsv(fname_motion_txt,fname_data_proc);

  % fields the pipeline computes, merged over the copied ones
  added = [];
  added.B0UnwarpApplied = tparms.B0unwarp_flag > 0; % logical, so jsonencode writes true/false
  if tparms.B0unwarp_flag > 0
    % NOTE: the pair actually used, flagged if it substitutes for a rejected one
    [~,fstem_fmap_for] = cmig_get_fileparts(fname_for);
    [~,fstem_fmap_rev] = cmig_get_fileparts(fname_rev);
    added.B0UnwarpFieldMapFor = fstem_fmap_for;
    added.B0UnwarpFieldMapRev = fstem_fmap_rev;
    added.B0UnwarpFieldMapFallback = fmap_fallback_flag > 0;
    if ~fmap_flag
      % B0FieldSource: the identifier of the pair actually used (for a
      %   fallback, the substitute pair's), not the rawdata one
      b0_id = mmil_getfield(for_info,'B0FieldIdentifier',[]);
      if isempty(b0_id), b0_id = mmil_getfield(rev_info,'B0FieldIdentifier',[]); end
      if ~isempty(b0_id), added.B0FieldSource = b0_id; end
      % omitted where the bold has no ShimSetting
      shim_match = check_shim_match(scan_info,for_info,rev_info);
      if ~isempty(shim_match)
        added.B0UnwarpFieldMapShimMatch = shim_match;
      end
      % head displacement (mm) from the field map: translation of the bold's
      %   registration to its pair (bold_unwarpB0 regref, mmil_reg_SEtoGE)
      M_GE_to_SE = load_regref(mmil_getfield(tparms,'fname_regref_B0uw',[]));
      if ~isempty(M_GE_to_SE)
        added.B0UnwarpFieldMapDisplacement = norm(M_GE_to_SE(1:3,4));
      end
      % fraction of the bold slab inside the forward field map's slab
      [M_for,volsz_for] = mmil_load_mgh_info(fname_for_tmp,[],tmpdir_for);
      added.B0UnwarpFieldMapCoverage = calc_fmap_coverage(M_raw,volsz_raw,M_for,volsz_for);
    end
  end
  % start of the between-scan registration kept, and the distance (mm)
  %   between the header and centroid starts' results (bold_estimate_motion)
  if ~fmap_flag && ~isempty(start_info)
    added.BetweenScanStart = start_info.start;
    if isfinite(start_info.distance)
      added.BetweenScanStartDistance = start_info.distance;
    end
  end
  added.RegistrationMatrixT1w = parms.M_T1_to_BOLD; % registration matrix (i.e., T1 to fMRI)
  added.SourceAffine = M_raw; % original vox2ras matrix
  added.Affine = M_proc; % processed vox2ras matrix
  added.SourceVolumeSize = volsz_raw; % raw volume dimensions (x,y,z,t)
  added.VolumeSize = volsz_proc; % processed volume dimensions (x,y,z,t)
  added.Description = set_description(fmap_flag,tparms);
  % Sources: the rawdata scan, then the rawdata field map pair it was unwarped with
  added.Sources = {raw_uri(parms,fname_data)};
  if tparms.B0unwarp_flag > 0
    added.Sources = unique(cat(1,added.Sources,...
      {raw_uri(parms,fname_for);raw_uri(parms,fname_rev)}),'stable');
  end
  added.SkullStripped = false; % no brain extraction is applied to the output
  added.BiasCorrectionApplied = false; % none is applied for func or fmap
  % NOTE: gradunwarp_flag is cleared above if the gradwarp type is missing
  added.GradUnwarpApplied = tparms.gradunwarp_flag > 0;
  % the frame-defining scan is 'orig'; everything else in that frame points at it
  fname_data_ref = mmil_getfield(parms,'fname_data_ref',[]);
  if ~isempty(fname_data_ref) && strcmp(fname_data,fname_data_ref)
    added.SpatialReference = 'orig';
  else
    added.SpatialReference = mmil_getfield(parms,'spatial_reference','orig');
  end

  % assemble the sidecar from the copy-list plus the fields above
  % NOTE: filtered only here, since check_fmap_pairs reads fields that are not
  %   written; see cmig_filter_bids_info for the copy-list
  info = cmig_filter_bids_info(scan_info,added);

  % the encoding directions name a voxel axis and its sign; with std_orient
  %   the output's axes can differ from the raw's (LAS -> LPS reverses j)
  info = reorient_encoding_dirs(info,M_raw,M_proc);

  % write output json
  fname_info_proc = regexprep(fname_data_proc,'\.nii\.gz','.json');

  % a bold not unwarped has no B0FieldSource; one unwarped without a
  %   pair identifier keeps none rather than the rawdata one
  if ~fmap_flag && ~isfield(added,'B0FieldSource') && isfield(info,'B0FieldSource')
    info = rmfield(info,'B0FieldSource');
  end

  % the field maps' IntendedFor is written by set_fmap_intended_for, from
  %   the bold runs actually unwarped with them, not the rawdata intent
  if isfield(info,'IntendedFor')
    info = rmfield(info,'IntendedFor');
  end

  write_json(fname_info_proc,info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = set_description(fmap_flag,tparms)
  % the processing applied to the scan, from the flags proc_scan passes to
  %   bold_correct_data
  if fmap_flag
    desc = 'EPI field map';
  else
    desc = 'BOLD series';
  end
  steps = {};
  % a field map is a single frame: its mc_flag pass is the registration to
  %   the session's reference field map, named at the end instead
  if tparms.mc_flag && ~fmap_flag, steps{end+1} = 'motion correction'; end
  if tparms.B0unwarp_flag > 0, steps{end+1} = 'B0 distortion correction (topup)'; end
  if tparms.gradunwarp_flag > 0, steps{end+1} = 'gradient nonlinearity correction'; end
  if length(steps) > 1
    desc = sprintf('%s with %s and %s',desc,strjoin(steps(1:end-1),', '),steps{end});
  elseif length(steps) == 1
    desc = sprintf('%s with %s',desc,steps{1});
  end
  if fmap_flag
    desc = [desc ', registered and resampled to the session''s reference field map'];
  else
    desc = [desc ', resampled to the session''s output grid'];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function set_fmap_intended_for(parms)
  % IntendedFor of each proc field map: the proc bold runs whose sidecars
  %   name it in B0UnwarpFieldMapFor/Rev (both scans of the pair); omitted for
  %   a field map no processed run used
  % NOTE: read from the sidecars on disk, so a rerun with forceflag 0 that
  %   leaves some bold or field map outputs in place stays consistent
  nfmap = length(parms.fnames_data_fmap);
  if nfmap==0, return; end
  fstems_fmap = cell(nfmap,1);
  for i=1:nfmap
    [~,fstems_fmap{i}] = cmig_get_fileparts(parms.fnames_data_fmap{i});
  end
  uris = cell(nfmap,1);
  for s=1:length(parms.fnames_data_func)
    [~,fstem] = cmig_get_fileparts(parms.fnames_data_func{s});
    fname_out = sprintf('%s/%s',parms.outdir,set_proc_fstem(fstem));
    if ~exist([fname_out '.nii.gz'],'file') || ~exist([fname_out '.json'],'file'), continue; end
    try
      info = jsondecode(fileread([fname_out '.json']));
    catch
      continue;
    end
    if ~isequal(mmil_getfield(info,'B0UnwarpApplied',false),true), continue; end
    pair = {mmil_getfield(info,'B0UnwarpFieldMapFor',''),mmil_getfield(info,'B0UnwarpFieldMapRev','')};
    idx = find(ismember(fstems_fmap,pair(cellfun(@ischar,pair))));
    for i=reshape(idx,1,[])
      uris{i}{end+1,1} = proc_uri(parms,'func',fstem);
    end
  end
  for i=1:nfmap
    fname_json = sprintf('%s/%s.json',parms.outdir_fmap,set_proc_fstem(fstems_fmap{i}));
    if ~exist(fname_json,'file'), continue; end
    replace_intended_for(fname_json,uris{i});
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function replace_intended_for(fname_json,uris)
  % set IntendedFor in a sidecar's text, or remove it when uris is empty
  % NOTE: edited as text, not decoded and re-encoded, since jsondecode turns
  %   a one-element array into a scalar and jsonencode would not restore it
  json = fileread(fname_json);
  json_old = json;
  pat = '"IntendedFor"\s*:\s*("[^"]*"|\[[^\]]*\])';
  entry = regexp(json,pat,'match','once');
  if ~isempty(uris)
    entry_new = ['"IntendedFor": ' jsonencode(uris)];
    if ~isempty(entry)
      json = strrep(json,entry,entry_new);
    else
      k = find(json=='{',1);
      json = [json(1:k) sprintf('\n  %s,',entry_new) json(k+1:end)];
    end
  elseif ~isempty(entry)
    % with the comma that joins it to the next field, or to the previous one
    json_new = regexprep(json,[pat '\s*,\s*'],'','once');
    if strcmp(json_new,json)
      json_new = regexprep(json,['\s*,\s*' pat],'','once');
    end
    json = json_new;
  end
  if strcmp(json,json_old), return; end
  fid = fopen(fname_json,'wt');
  if fid<0, error('failed to open %s for writing',fname_json); end
  fprintf(fid,'%s',json);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shim_match = check_shim_match(scan_info,for_info,rev_info)
  % whether the scan's ShimSetting equals the pair's, compared as
  %   cmig_select_bids_fmap does; empty if the scan has none
  shim_match = [];
  shim = mmil_getfield(scan_info,'ShimSetting',[]);
  if isempty(shim) || all(isnan(shim(:))), return; end
  shim_match = isequal(for_info.ShimSetting,shim) && ...
               isequal(rev_info.ShimSetting,shim);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_GE_to_SE = load_regref(fname_regref)
  % registration matrix saved by bold_unwarpB0, empty if absent
  M_GE_to_SE = [];
  if isempty(fname_regref) || ~exist(fname_regref,'file'), return; end
  vars = whos('-file',fname_regref);
  if ~ismember('M_GE_to_SE',{vars.name}), return; end
  tmp = load(fname_regref,'M_GE_to_SE');
  M_GE_to_SE = tmp.M_GE_to_SE;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coverage = calc_fmap_coverage(M_scan,volsz_scan,M_fmap,volsz_fmap)
  % fraction of the scan's extent along its slice axis (third voxel axis)
  %   that lies inside the field map's extent, in scanner coordinates
  n = M_scan(1:3,3)/norm(M_scan(1:3,3));
  [s0,s1] = calc_slab_extent(M_scan,volsz_scan,n);
  [f0,f1] = calc_slab_extent(M_fmap,volsz_fmap,n);
  coverage = max(0,min(s1,f1)-max(s0,f0))/(s1-s0);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x0,x1] = calc_slab_extent(M,volsz,n)
  % range of the voxel grid's outer corners (1-based vox2ras) projected onto n
  [i,j,k] = ndgrid([0.5 volsz(1)+0.5],[0.5 volsz(2)+0.5],[0.5 volsz(3)+0.5]);
  x = n'*(M(1:3,:)*[i(:)';j(:)';k(:)';ones(1,8)]);
  x0 = min(x); x1 = max(x);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_json(fname,info)
  % NOTE: jsonencode writes the shortest round-trip form of each number
  % NOTE: PrettyPrint, as dcm2niix writes the rawdata sidecars, needs R2021a (9.10)
  if verLessThan('matlab','9.10')
    json = jsonencode(info);
  else
    json = jsonencode(info,'PrettyPrint',true);
  end
  % B0UnwarpFieldMapCoverage as %.6f, so it is a float in every case
  %   (jsonencode writes 1.0 as 1), with the postqc cover_prop's precision
  val = mmil_getfield(info,'B0UnwarpFieldMapCoverage',[]);
  if isnumeric(val) && isscalar(val) && isfinite(val)
    json = regexprep(json,'"B0UnwarpFieldMapCoverage"\s*:\s*[^,\s}]+',...
      ['"B0UnwarpFieldMapCoverage": ' sprintf('%.6f',val)]);
  end
  fid = fopen(fname,'wt');
  if fid<0, error('failed to open %s for writing',fname); end
  fprintf(fid,'%s',json);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode,errmsg] = copy_events_file(fname_events,parms)
  fname_out = []; errcode = 0; errmsg = [];
  if isempty(fname_events), return; end
  [~,fstem,fext] = cmig_get_fileparts(fname_events);
  fname_out = sprintf('%s/%s%s',parms.outdir,fstem,fext);
  [s,r] = mmil_unix(sprintf('cp %s %s',fname_events,fname_out));
  if s
    errcode = 1;
    errmsg = sprintf('failed to copy events file %s to %s:\n%s',...
      fname_events,parms.outdir,r);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  mmil_unix(sprintf('chmod u+w %s',fname_out));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,fname_B0dx,fname_motion_txt,start_info,errcode,errmsg] = correct_data(fname_in,parms)
  fname_out = []; fname_B0dx = []; fname_motion_txt = []; start_info = [];
  errcode = 0; errmsg = [];
  args = mmil_parms2args(parms,parms.correct_tags);
  try
    [fname_out,fname_B0dx,fname_motion_txt,start_info] = bold_correct_data(fname_in,args{:});
  catch me
    fprintf('%s: WARNING: bold_correct_data failed for %s:\nerror in %s at line %d:\n%s\n',...
      mfilename,fname_in,me.stack(1).name,me.stack(1).line,me.message);
    % delete tmp dir
    tpath = fileparts(fname_in);
    tmpdir = [tpath '/tmp_BOLD_corr'];
    if exist(tmpdir,'dir')
      fprintf('%s: removing tmp dir %s...\n',mfilename,tmpdir);
      cmd = sprintf('rm -r %s',tmpdir);
      [status,result,err] = jsystem(cmd);
      if status
        fprintf('%s: WARNING: deletion of temporary directory %s failed:\n%s\n',mfilename,tmpdir,err);
      end
    end
    % try again
    fprintf('%s: re-attempting bold_correct_data...\n',mfilename);
    try
      [fname_out,fname_B0dx,fname_motion_txt,start_info] = bold_correct_data(fname_in,args{:});
    catch me
      errcode = 1;
      errmsg = sprintf('bold_correct_data failed for %s:\nerror in %s at line %d:\n%s',...
        fname_in,me.stack(1).name,me.stack(1).line,me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_motion_tsv(fname_motion_txt,fname_data_proc)
  if isempty(fname_motion_txt) || ~exist(fname_motion_txt,'file'), return; end
  fname_out = regexprep(fname_data_proc,'_desc-preproc_(bold|epi)\.nii\.gz$','_desc-motion_timeseries.tsv');
  if strcmp(fname_out,fname_data_proc)
    fprintf('%s: WARNING: could not determine motion TSV path from %s\n',mfilename,fname_data_proc);
    return;
  end
  col_names = {'t_indx','rot_z','rot_x','rot_y','trans_z','trans_x','trans_y'};
  fid_in = fopen(fname_motion_txt,'rt');
  if fid_in < 0
    fprintf('%s: WARNING: failed to open motion file %s\n',mfilename,fname_motion_txt);
    return;
  end
  data = textscan(fid_in,'%f %f %f %f %f %f %f %f %f');
  fclose(fid_in);
  fid_out = fopen(fname_out,'wt');
  if fid_out < 0
    fprintf('%s: WARNING: failed to open %s for writing\n',mfilename,fname_out);
    return;
  end
  fprintf(fid_out,'%s\n',strjoin(col_names,'\t'));
  for r = 1:length(data{1})
    fprintf(fid_out,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',...
      data{1}(r),data{2}(r),data{3}(r),data{4}(r),data{5}(r),data{6}(r),data{7}(r));
  end
  fclose(fid_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = export_data(fname_in,fstem,outdir,TR)
  % set new name
  fname_out = sprintf('%s/%s.nii.gz',outdir,fstem);
  % convert from mgz to nifti
  % NOTE: the mgz carries TR 0; -tr takes msec and writes seconds to pixdim[4]
  options = [];
  if ~isempty(TR) && ~isnan(TR) && TR > 0
    options = sprintf('-tr %g',1000*TR);
  end
  fs_mri_convert(fname_in,fname_out,'options',options);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_b0dx(fname_B0dx,fstem_for,fstem_rev,tmpdir_for,for_info,outfix,parms)
  % export the B0 displacement field for one forward/reverse pair, in mm and
  %   in the reference field map's frame at output resolution
  % NOTE: B0dx (topup's field times the readout time) is in voxels of its
  %   own grid, so it is converted to mm before resampling, since the grids
  %   differ by vendor
  % NOTE: not smoothed; topup's warp resolution already band-limits it in mm

  % set output file names
  fname_out = set_b0dx_fname(fstem_for,fstem_rev,parms);
  fname_out_json = regexprep(fname_out,'\.nii\.gz$','.json');
  if exist(fname_out,'file') && ~parms.forceflag
    if ~exist(fname_out_json,'file')
      write_b0dx_json(fname_out_json,fstem_for,fstem_rev,parms);
    end
    return;
  end

  % determine phase-encode axis of the field's own grid
  pe_axis = get_pe_axis(for_info);
  if isempty(pe_axis)
    fprintf('%s: WARNING: skipping B0dx export for %s because the phase-encode axis could not be determined\n',...
      mfilename,fstem_for);
    return;
  end

  % between-scan motion correction matrix for the forward field map scan
  fname_motion = sprintf('%s/%s_%s_motion_inter.mat',tmpdir_for,fstem_for,outfix);
  if ~exist(fname_motion,'file')
    fprintf('%s: WARNING: skipping B0dx export for %s because %s was not found\n',...
      mfilename,fstem_for,fname_motion);
    return;
  end
  tmp = load(fname_motion,'M_motion');
  if ~isfield(tmp,'M_motion') || isempty(tmp.M_motion)
    fprintf('%s: WARNING: skipping B0dx export for %s because %s contains no M_motion\n',...
      mfilename,fstem_for,fname_motion);
    return;
  end
  M_motion = tmp.M_motion{1};

  % reference field map scan defines the output frame
  fname_ref = mmil_getfield(parms,'fname_fmap_for_ref',[]);
  if isempty(fname_ref) || ~exist(fname_ref,'file')
    fprintf('%s: WARNING: skipping B0dx export for %s because the reference field map is missing\n',...
      mfilename,fstem_for);
    return;
  end
  [M_ref,volsz_ref] = mmil_load_mgh_info(fname_ref,parms.forceflag,tmpdir_for);

  % convert from voxels of the field's own grid to mm, then rigid resample
  %   into the reference field map's frame and regrid to output resolution
  [~,fstem_B0dx] = fileparts(fname_B0dx);
  fname_tmp = sprintf('%s/%s_mm_%s%s',tmpdir_for,fstem_B0dx,outfix,parms.outext);
  if ~exist(fname_tmp,'file') || parms.forceflag
    [vol_dx,M_dx] = fs_load_mgh(fname_B0dx);
    if ndims(vol_dx)>3, vol_dx = vol_dx(:,:,:,1); end
    vol_dx = vol_dx * norm(M_dx(1:3,pe_axis));
    % signed along the output frame's phase-encode axis: with std_orient,
    %   that axis can point the other way from the field's own (LAS -> LPS)
    u_pe = M_dx(1:3,pe_axis)/norm(M_dx(1:3,pe_axis));
    proj = M_ref(1:3,1:3)' * u_pe;
    [~,k_pe] = max(abs(proj));
    vol_dx = sign(proj(k_pe)) * vol_dx;
    vol_dx = mmil_resample_vol(vol_dx,M_dx,...
      'M_reg',M_motion,'M_ref',M_ref,'nvox_ref',volsz_ref(1:3),'bclamp',0);
    fs_save_mgh(vol_dx,fname_tmp,M_ref);
  end

  % write output displacement field
  fs_mri_convert(fname_tmp,fname_out);
  write_b0dx_json(fname_out_json,fstem_for,fstem_rev,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_b0dx_json(fname,fstem_for,fstem_rev,parms)
  % Sources: the processed field map pair the field was estimated from;
  %   SpatialReference: the reference field map, whose frame it is in
  info = [];
  info.Description = ['B0 displacement field estimated by topup from a',...
    ' forward and reverse field map pair, in mm along the phase-encode axis'];
  info.Sources = {proc_uri(parms,'fmap',fstem_for);proc_uri(parms,'fmap',fstem_rev)};
  info.SpatialReference = mmil_getfield(parms,'spatial_reference','orig');
  if isempty(info.SpatialReference), info.SpatialReference = 'orig'; end
  write_json(fname,info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = reorient_encoding_dirs(info,M_in,M_out)
  % restate PhaseEncodingDirection and SliceEncodingDirection ('i','j','k',
  %   with '-' for decreasing) from the voxel axes of M_in to those of M_out:
  %   the same direction in the world, as the output axis nearest to it
  letters = 'ijk';
  fields = {'PhaseEncodingDirection','SliceEncodingDirection'};
  for f = 1:length(fields)
    val = mmil_getfield(info,fields{f},'');
    if isempty(val) || ~ischar(val), continue; end
    a = strfind(letters,val(1));
    if isempty(a), continue; end
    w = M_in(1:3,a);
    if any(val=='-'), w = -w; end
    proj = M_out(1:3,1:3)' * w;
    [~,b] = max(abs(proj));
    if proj(b) < 0
      info.(fields{f}) = [letters(b) '-'];
    else
      info.(fields{f}) = letters(b);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pe_axis = get_pe_axis(scan_info)
  pe_axis = [];
  pe = mmil_getfield(scan_info,'PhaseEncodingDirection','');
  if isempty(pe) || ~ischar(pe), return; end
  switch pe(1)
    case 'i'
      pe_axis = 1;
    case 'j'
      pe_axis = 2;
    case 'k'
      pe_axis = 3;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = copy_first_frame(fname_in,parms)
  [fpath,fstem,fext] = fileparts(fname_in);
  fname_out = sprintf('%s/%s_f0%s',fpath,fstem,fext);
  if ~exist(fname_out,'file') || parms.forceflag
    [vol,M]=fs_load_mgh(fname_in,[],1);
    fs_save_mgh(vol,fname_out,M);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B] = reconcile_fields(A,B)
  if isempty(A) || isempty(B), return; end
  fnamesA = fieldnames(A);
  fnamesB = fieldnames(B);
  fnames = setdiff(fnamesB,fnamesA);
  for f=1:length(fnames)
    A(1).(fnames{f}) = [];
  end
  fnames = setdiff(fnamesA,fnamesB);
  for f=1:length(fnames)
    B(1).(fnames{f}) = [];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_B0dx = set_b0dx_fname(fstem_for,fstem_rev,parms)
  % name the exported B0dx by the pair it came from, keeping AP and PA as fixed
  %   slots filled from each scan's own dir entity rather than from for/rev
  %   order; the forward scan is dir-PA for ABCD, so filling by order would put
  %   the PA run in the AP slot
  tok = regexp(fstem_for,'_acq-([^_]+)','tokens','once');
  if ~isempty(tok), acq_label = tok{1}; else acq_label = 'unknown'; end
  run_ap = ''; run_pa = '';
  fstems = {fstem_for,fstem_rev};
  for i=1:length(fstems)
    tok = regexp(fstems{i},'_dir-([^_]+)','tokens','once');
    if isempty(tok), continue; end
    dir_label = upper(tok{1});
    tok = regexp(fstems{i},'_run-([^_]+)','tokens','once');
    if ~isempty(tok), run_label = tok{1}; else run_label = ''; end
    switch dir_label
      case 'AP', run_ap = run_label;
      case 'PA', run_pa = run_label;
    end
  end
  fname_B0dx = sprintf('%s/%s_%s_acq-%s_desc-AP%sPA%s_B0dx.nii.gz',...
    parms.outdir_fmap,parms.participant_id,parms.session_id,...
    acq_label,run_ap,run_pa);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = set_proc_fstem(fstem)
  % the proc output of a rawdata scan: desc-preproc before the suffix
  %   (task-mid_run-01_bold -> task-mid_run-01_desc-preproc_bold)
  fstem = regexprep(fstem,'_([A-Za-z0-9]+)$','_desc-preproc_$1');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = raw_uri(parms,fname)
  % bids:rawdata: URI of a rawdata scan; rawdata is the input root's name in
  %   DatasetLinks (RAW_LABEL in cmig_write_bids_dataset_description.py)
  [fpath,fstem,fext] = fileparts(fname);
  [~,subdir] = fileparts(fpath);
  uri = sprintf('bids:rawdata:%s/%s/%s/%s%s',...
    parms.participant_id,parms.session_id,subdir,fstem,fext);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = proc_uri(parms,subdir,fstem)
  % bids:: URI of the proc output of the rawdata scan with file stem fstem
  uri = sprintf('bids::%s/%s/%s/%s.nii.gz',...
    parms.participant_id,parms.session_id,subdir,set_proc_fstem(fstem));
return;
