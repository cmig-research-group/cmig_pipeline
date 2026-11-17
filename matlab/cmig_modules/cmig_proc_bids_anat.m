function [outdir,errcode,errmsg] = cmig_proc_bids_anat(participant_id,session_id,input_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_proc_bids_anat(participant_id,session_id,input_rootdir,output_rootdir,[options])
%
% Purpose: Process anatomical MRI data from BIDS format
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
%   'select_flag': [0|1] whether to select one scan to process
%     based on shortest QC notes in scans.tsv file or last scan acquired
%     {default = 1}
%   'atlas_flag': [0|1] whether to register to atlas
%     {default = 1}
%   'native_flag': [0|1] whether to keep original orientation and resolution
%     {default = 0}
%   'gradunwarp_flag': [0|1|2] whether to correct for gradient warping
%     if 0, do not perform gradwarp
%     if 1, perform gradwarp if gradwarp type available and supported
%     if 2, quit if gradwarp type unavailable or unsupported
%     {default = 2}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     will have form of {output_rootdir}/{participant_id}/{session_id}/anat
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  08/04/2025 by Don Hagler
% Prev Mod: 08/06/2025 by Don Hagler
% Prev Mod: 08/07/2025 by Don Hagler
% Prev Mod: 08/12/2025 by Don Hagler
% Prev Mod: 08/14/2025 by Don Hagler
% Prev Mod: 09/04/2025 by Don Hagler
% Prev Mod: 03/09/2026 by Don Hagler
% Prev Mod: 03/12/2026 by Don Hagler
% Prev Mod: 04/24/2026 by Don Hagler
% Prev Mod: 04/29/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 06/04/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/10/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/19/2026 by Don Hagler
% Prev Mod: 08/26/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/12/2026 by Don Hagler
% Last Mod: 11/13/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; errcode = 0; errmsg = [];
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

% write status file
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-proc',errcode,errmsg,[],startdate,procflag);

outdir = parms.outdir;
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
    'select_flag',true,[false true],...
    'atlas_flag',true,[false true],...
    'native_flag',false,[false true],...
    'gradunwarp_flag',2,[0 1 2],...
    'cleanup_flag',true,[false true],...
    'forceflag',false,[false true],...
  });
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];

  % set output directory
  parms.outdir = sprintf('%s/%s/%s/anat',parms.output_rootdir,parms.participant_id,parms.session_id);

  % get input file names for T1w scans
  [parms.fnames_data_t1,parms.fnames_info_t1] = get_input_data(parms.participant_id,parms.session_id,...
    parms.input_rootdir,'T1w',parms.select_flag);

  % get input file names for T2w scans
  [parms.fnames_data_t2,parms.fnames_info_t2] = get_input_data(parms.participant_id,parms.session_id,...
    parms.input_rootdir,'T2w',parms.select_flag);

  % quit if no input files to process
  if isempty(parms.fnames_data_t1) && isempty(parms.fnames_data_t2)
    errcode = 1;
    errmsg = 'not processed: no usable T1w or T2w scans to process';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fnames_data,fnames_info] = get_input_data(participant_id,session_id,input_rootdir,series_type,select_flag)
  % get input file names
  [fnames_data,fnames_info,scan_info] = cmig_get_bids_data(participant_id,session_id,input_rootdir,'anat',series_type);
  % check that there are input files to process
  if isempty(fnames_data), return; end
  if select_flag
    % choose which scan to process based on scans_info
    idx = cmig_select_scan(scan_info);
    fnames_data = fnames_data(idx);
    fnames_info = fnames_info(idx);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = proc_data(parms)
  errcode = 0; errmsg = [];

  % create output directory
  mmil_mkdir(parms.outdir);

  % process T1w scans
  fname_ref = [];
  nscans_t1 = length(parms.fnames_data_t1);
  errcodes_t1 = zeros(nscans_t1,1);
  errmsgs_t1 = cell(nscans_t1,1);
  for s=1:nscans_t1
    [fname_data_t1_proc,fname_info_t1_proc,errcodes_t1(s),errmsgs_t1{s}] = proc_scan(...
      parms.fnames_data_t1{s},parms.fnames_info_t1{s},...
      parms,fname_ref,'T1w');
    if isempty(fname_ref), fname_ref = fname_data_t1_proc; end
  end
  
  % process T2w scans
  nscans_t2 = length(parms.fnames_data_t2);
  errcodes_t2 = zeros(nscans_t2,1);
  errmsgs_t2 = cell(nscans_t2,1);
  for s=1:nscans_t2
    [fname_data_t2_proc,fname_info_t2_proc,errcodes_t2(s),errmsgs_t2{s}] = proc_scan(...
      parms.fnames_data_t2{s},parms.fnames_info_t2{s},...
      parms,fname_ref,'T2w');
  end

  % set errcode and errmsg based on errcodes and errmsgs
  [errcode,errmsg] = cmig_summarize_errcodes(...
    cat(1,errcodes_t1,errcodes_t2),...
    cat(1,errmsgs_t1,errmsgs_t2));
  if errcode, return; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_data_proc,fname_info_proc,errcode,errmsg] = proc_scan(fname_data,fname_info,parms,fname_ref,scantype)
  fname_data_proc = []; fname_info_proc = [];
  errcode = 0; errmsg = [];

  if ~mmil_check_nargs(nargin,4), return; end
  if ~exist('fname_ref','var'), fname_ref = []; end
  if ~exist('scantype','var') || isempty(scantype), scantype = 'T1w'; end

  % quit if no input
  if isempty(fname_data)
    fprintf('%s: no usable data for scan type %s\n',mfilename,scantype);
    return;
  end

  % load info file
  try
    scan_info = jsondecode(fileread(fname_info));
  catch me
    error('failed to load json file %s:\n%s',fname_info,me.message);
  end

  % create tmp dir nested under parms.tmpdir
  [~,fstem] = cmig_get_fileparts(fname_data);
  tmpdir = sprintf('%s/%s',parms.tmpdir,fstem);
  mmil_mkdir(tmpdir)
    
  % import data (convert to mgz)
  fname_tmp = cmig_import_to_mgz(fname_data,tmpdir,[],parms.forceflag);

  % get volume info for input volume
  [M_raw,volsz_raw] = mmil_load_mgh_info(fname_tmp,[],tmpdir);

  % gradient nonlinearity distortion correction
  % NOTE: grad_unwarp returns fname_in unchanged when it skips -- Philips, a
  %   scanner that already did 3D, or an unknown gradwarp type -- so a changed
  %   file name is what says the correction ran
  gradunwarp_applied = false;
  if parms.gradunwarp_flag
    fname_pre = fname_tmp;
    [fname_tmp,errcode,errmsg] = grad_unwarp(fname_tmp,scan_info,parms.gradunwarp_flag,parms.forceflag);
    if errcode, return; end
    gradunwarp_applied = ~strcmp(fname_tmp,fname_pre);
  end

  atlas_reg_failed = NaN; atlas_reg_corr = NaN; atlas_reg_fvox_miss = NaN;
  atlas_reg_corr_fail = NaN; atlas_reg_fvox_miss_fail = NaN;
  % NOTE: bias_corr also returns fname_in unchanged when it skips, which it
  %   does for T2w without a T1w reference
  biascorr_applied = false;
  switch scantype
    case 'T1w'
      % apply white matter bias correction
      fname_pre = fname_tmp;
      [fname_tmp,errcode,errmsg] = bias_corr(fname_tmp,[],scantype,parms.tmpdir,parms.forceflag);
      if errcode, return; end
      biascorr_applied = ~strcmp(fname_tmp,fname_pre);
      if isempty(fname_ref)
        % register to atlas and resample, create mask
        % clip intensities only if wm bias correction normalized them
        if biascorr_applied
          image_range = get_image_range(scantype);
        else
          image_range = [];
        end
        [fname_tmp,fname_tmp_mask,M_reg,atlas_reg_failed,atlas_reg_corr,atlas_reg_fvox_miss,...
          atlas_reg_corr_fail,atlas_reg_fvox_miss_fail] = cmig_reg_T1w_to_atlas(fname_tmp,...
          'atlas_flag',parms.atlas_flag,'native_flag',parms.native_flag,'cleanupflag',false,...
          'image_range',image_range,'forceflag',parms.forceflag);
      else
        % register to reference and resample
        [fname_tmp,M_reg,errcode,errmsg] = reg_to_ref(fname_tmp,fname_ref,scantype,parms.tmpdir,biascorr_applied,parms.forceflag);
        if errcode, return; end
      end
    case 'T2w'
      % apply white matter bias correction
      fname_pre = fname_tmp;
      [fname_tmp,errcode,errmsg] = bias_corr(fname_tmp,fname_ref,scantype,parms.tmpdir,parms.forceflag);
      if errcode, return; end
      biascorr_applied = ~strcmp(fname_tmp,fname_pre);
      % register to reference and resample
      [fname_tmp,M_reg,errcode,errmsg] = reg_to_ref(fname_tmp,fname_ref,scantype,parms.tmpdir,biascorr_applied,parms.forceflag);
      if errcode, return; end
    otherwise
      fprintf('%s: WARNING: no registration for scantype %s\n',mfilename,scantype);
      return;
  end

  % get volume info for output volume
  [M_proc,volsz_proc] = mmil_load_mgh_info(fname_tmp,[],tmpdir);

  % save output data file
  fname_data_proc = export_data(fname_tmp,get_proc_stem(fstem),parms.outdir);

  % save output reference mask file with BIDS-compliant name
  if isempty(fname_ref)
    [fpath_tmp,fstem_tmp,fext_tmp] = cmig_get_fileparts(fname_tmp);
    fname_tmp_mask = sprintf('%s/%s_mask%s',fpath_tmp,fstem_tmp,fext_tmp);
    fname_mask_proc = export_data(fname_tmp_mask,get_mask_stem(fstem),parms.outdir);
    write_mask_json(fname_mask_proc,fname_data_proc,atlas_reg_failed,parms);
  end

  % fields the pipeline computes, merged over the copied ones
  added = [];
  added.Description = set_description(scantype,gradunwarp_applied,biascorr_applied,...
    fname_ref,M_reg,atlas_reg_failed,parms);
  % Sources: the rawdata scan, then the processed T1w it was corrected and
  %   registered with, if any
  [fpath,fname_raw,fext] = fileparts(fname_data);
  [~,subdir] = fileparts(fpath);
  % rawdata is the input root's name in DatasetLinks (RAW_LABEL in
  %   cmig_write_bids_dataset_description.py)
  added.Sources = {sprintf('bids:rawdata:%s/%s/%s/%s%s',...
    parms.participant_id,parms.session_id,subdir,fname_raw,fext)};
  if ~isempty(fname_ref)
    [~,fstem_ref] = cmig_get_fileparts(fname_ref);
    added.Sources{end+1,1} = sprintf('bids::%s/%s/anat/%s.nii.gz',...
      parms.participant_id,parms.session_id,fstem_ref);
  end
  added.RegistrationMatrixT1w = M_reg; % registration matrix (e.g., atlas to subject, reference to raw)
  added.SourceAffine = M_raw; % original vox2ras matrix
  added.Affine = M_proc; % processed vox2ras matrix
  added.SourceVolumeSize = volsz_raw; % raw volume dimensions (x,y,z,t)
  added.VolumeSize = volsz_proc; % processed volume dimensions (x,y,z,t)
  added.AtlasRegistrationFailed = atlas_reg_failed; % [0|1|NaN] whether registration to atlas failed (NaN if reg to atlas not attempted)
  added.AtlasRegistrationCorrelation = atlas_reg_corr; % correlation between atlas and subject brain, for the registration applied
  added.AtlasRegistrationFracVoxMissing = atlas_reg_fvox_miss; % fraction of voxels missing, for the registration applied
  added.AtlasRegistrationFailedCorrelation = atlas_reg_corr_fail; % correlation of the failed registration (NaN if not failed)
  added.AtlasRegistrationFailedFracVoxMissing = atlas_reg_fvox_miss_fail; % fraction of voxels missing of the failed registration (NaN if not failed)
  added.GradUnwarpApplied = gradunwarp_applied; % whether gradient unwarp ran
  added.BiasCorrectionApplied = biascorr_applied; % whether wm bias correction ran
  added.SkullStripped = false; % the brain mask is a separate file, not applied
  % NOTE: the atlas defines the space, so a scan resampled into it points at
  %   the copy shipped with the derivative dataset
  % NOTE: set from the processing outcome, not the requested flags: the
  %   atlas-reg-failure fallback (resampled to 256 LIA without atlas
  %   registration) is a grid change, not a space change, and a scan whose
  %   registration was soft-skipped stays in its own space
  if isempty(fname_ref)
    if parms.atlas_flag && ~parms.native_flag && atlas_reg_failed==0
      added.SpatialReference = 'bids::atlases/atlas-ABCD_T1w.nii.gz';
    else
      added.SpatialReference = 'orig';
    end
  elseif isempty(M_reg)
    % registration to reference soft-skipped
    added.SpatialReference = 'orig';
  else
    % resampled into register with reference: inherit its atlas space, or
    %   point at the reference itself when it is in its own space
    added.SpatialReference = get_ref_space(fname_ref,parms);
  end

  % assemble the sidecar from the copy-list plus the fields above
  % NOTE: fields not on the copy-list are not written; see
  %   cmig_filter_bids_info
  info = cmig_filter_bids_info(scan_info,added);

  % write output json
  fname_info_proc = regexprep(fname_data_proc,'\.nii\.gz','.json');
  write_json(fname_info_proc,info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_mask_json(fname_mask,fname_data_proc,atlas_reg_failed,parms)
  % the mask lies on the grid of the processed T1w it was made with
  [~,fstem_data] = cmig_get_fileparts(fname_data_proc);
  uri_data = sprintf('bids::%s/%s/anat/%s.nii.gz',...
    parms.participant_id,parms.session_id,fstem_data);
  info = [];
  % NOTE: cmig_reg_T1w_to_atlas falls back to its mmil_brainmask mask when
  %   the atlas registration fails or is not attempted
  if atlas_reg_failed==0
    info.Description = 'Brain mask from the ABCD atlas brain mask, transformed by the registration of the T1w to the atlas';
    info.Sources = {uri_data,'bids::atlases/atlas-ABCD_desc-brain_mask.nii.gz'};
  else
    info.Description = 'Brain mask from mmil_brainmask on the T1w';
    info.Sources = {uri_data};
  end
  info.Type = 'Brain';
  info.SpatialReference = uri_data;
  info.SkullStripped = false;
  write_json(regexprep(fname_mask,'\.nii\.gz$','.json'),info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = set_description(scantype,gradunwarp_applied,biascorr_applied,...
  fname_ref,M_reg,atlas_reg_failed,parms)
  % the processing applied to the scan, from the outcomes proc_scan records
  steps = {};
  if gradunwarp_applied, steps{end+1} = 'gradient nonlinearity correction'; end
  if biascorr_applied, steps{end+1} = 'white matter bias correction'; end
  if isempty(fname_ref)
    if parms.atlas_flag && ~parms.native_flag && atlas_reg_failed==0
      steps{end+1} = 'registration to the ABCD atlas';
    end
  elseif ~isempty(M_reg)
    steps{end+1} = 'registration to the processed T1w';
  end
  desc = sprintf('%s image',scantype);
  if length(steps) > 1
    desc = sprintf('%s with %s and %s',desc,strjoin(steps(1:end-1),', '),steps{end});
  elseif length(steps) == 1
    desc = sprintf('%s with %s',desc,steps{1});
  else
    desc = [desc ', with no correction or registration applied'];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_json(fname_json,info)
  % NOTE: jsonencode has no FloatFormat equivalent, using a shortest
  %   round-trip representation instead of the former '%0.15g'
  % NOTE: PrettyPrint was added in R2021a (MATLAB 9.10) and errors as an
  %   unrecognized parameter before that. Guarded so an older MATLAB writes a
  %   single line rather than failing; dcm2niix pretty-prints the rawdata
  %   sidecars, so the indented form is what matches the inputs
  if verLessThan('matlab','9.10')
    json = jsonencode(info);
  else
    json = jsonencode(info,'PrettyPrint',true);
  end
  fid = fopen(fname_json,'wt');
  if fid<0, error('failed to open %s for writing',fname_json); end;
  fprintf(fid,'%s',json);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem_proc = get_proc_stem(fstem)
  % rawdata stem to processed output stem: desc-preproc before the suffix
  fstem_proc = regexprep(fstem,'_([^_]+)$','_desc-preproc_$1');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = get_raw_stem(fstem_proc)
  % processed output stem to rawdata stem, which names the per-scan tmpdir
  fstem = regexprep(fstem_proc,'_desc-preproc(_[^_]+)$','$1');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mask_stem = get_mask_stem(fstem)
  % rawdata stem to brain mask stem, e.g. {s}_run-01_space-T1w_desc-brain_mask
  tok = regexp(fstem,'^(.+)_([^_]+)$','tokens','once');
  mask_stem = sprintf('%s_space-%s_desc-brain_mask',tok{1},tok{2});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode,errmsg] = grad_unwarp(fname_in,scan_info,gradunwarp_flag,forceflag)
  fname_out = fname_in;
  errcode = 0; errmsg = [];
  if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end
  % skip Philips
  if ~isempty(regexpi(scan_info.Manufacturer,'Philips'))
    fprintf('%s: NOTE: grad unwarp not required for manufacturer %s\n',...
      mfilename,scan_info.Manufacturer);
    return;
  end
  % set output name
  % NOTE: fname_out keeps the input name until past the skips, so skip
  %   branches return fname_in and a changed name means the correction ran
  [fpath,fstem,fext] = fileparts(fname_in);
  fname_uw = sprintf('%s/%s_uw.mgz',fpath,fstem);
  if ~exist(fname_uw,'file') || forceflag
    % get gradwarp info from scan info
    [gwtype,unwarpflag,isoctrflag] = cmig_get_bids_gradwarpinfo(scan_info);
    % skip if 3D grad unwarp done on scanner
    % NOTE: mmil_getfield because unwarpflag is also NaN when the sidecar has no
    %   NonlinearGradientCorrection field at all, as Philips sidecars do not
    if isnan(unwarpflag)
      fprintf('%s: NOTE: grad unwarp not required for NonlinearGradientCorrectionType %s\n',...
        mfilename,mmil_getfield(scan_info,'NonlinearGradientCorrectionType','(absent)'));
      return;
    end
    % skip if gwtype unknown
    if isnan(gwtype)
      errmsg = 'missing gradwarp type necessary for grad unwarp';
      if gradunwarp_flag>1
        errcode = 1;
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      else
        fprintf('%s: WARNING: %s\n',mfilename,errmsg);
      end
      return;
    end
    % apply gradienty nonlinearity distortion correction
    fprintf('%s: unwarping %s...\n',mfilename,fname_in);
    [ctx_vol,mr_parms] = ctx_load_mgh(fname_in);
    ctx_voluw = ctx_unwarp_grad(ctx_vol,gwtype,unwarpflag,isoctrflag);
    ctx_save_mgh(ctx_voluw,fname_uw,mr_parms);
  end
  fname_out = fname_uw;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function image_range = get_image_range(scantype)
  if strcmp(scantype,'T2w')
    image_range = [0,2500];
  else
    image_range = [0,255];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,errcode,errmsg] = bias_corr(fname_in,fname_ref,scantype,parent_tmpdir,forceflag)
  fname_out = fname_in;
  errcode = 0; errmsg = [];

  if ~mmil_check_nargs(nargin,1), return; end
  if ~exist('fname_ref','var'), fname_ref = []; end
  if ~exist('scantype','var') || isempty(scantype), scantype = 'T1w'; end
  if ~exist('parent_tmpdir','var'), parent_tmpdir = []; end
  if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end

  % skip if scantype is T2 and fname_ref is empty or does not exist
  if strcmp(scantype,'T2w') && (isempty(fname_ref) || ~exist(fname_ref,'file'))
    errmsg = sprintf('skipping white matter bias correction for %s because no T1w reference',fname_in);
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end
   
  [fpath,fstem,fext] = fileparts(fname_in);

  % set parameters
  p = [];
  p.erode_wm_flag = false;
  p.dilate_bm_flag = false;
  p.n3flag = false;
  p.vol_target = 110;
  p.ratio_range = [];
  p.image_range = get_image_range(scantype);
  p.outlier_thresh = 0.1;
  p.outlier_edge_flag = true;
  p.lambda1 = 20;
  p.lambda2 = 1e3;
  p.postsmooth_flag = true;
  p.atlasdir = [];
  p.atlasname = 'T1_Atlas/T1_atlas_ABCD';
  p.atlasdir_T1T2 = 'T1T2';
  p.tmpdir = sprintf('%s/wmbc',fpath);
  p.cleanup_flag = false;
  p.verbose = false;
  p.ext = '.mgz';
  p.forceflag = forceflag;

  % set output file names
  p.fname_out = sprintf('%s/%s_wmbc%s',fpath,fstem,fext);
  p.fname_wm = sprintf('%s/%s_wm%s',fpath,fstem,fext);
  p.fname_bm = sprintf('%s/%s_bm%s',fpath,fstem,fext);
  p.fname_gm = sprintf('%s/%s_gm%s',fpath,fstem,fext);

  if strcmp(scantype,'T2w')
    [fpath_ref,fstem_ref,fext_ref] = cmig_get_fileparts(fname_ref);
    fstem_ref = get_raw_stem(fstem_ref);
    if ~isempty(parent_tmpdir)
      tmpdir_ref = sprintf('%s/%s',parent_tmpdir,fstem_ref);
    else
      tmpdir_ref = strcat(fpath_ref,'/',fstem_ref);
    end

    % identify reference input and masks
    infix_ref = '_uw';
    fname_ref_tmp = sprintf('%s/%s%s%s',tmpdir_ref,fstem_ref,infix_ref,fext);
    if ~exist(fname_ref_tmp,'file'), infix_ref = ''; end
    p.fname_T1 = sprintf('%s/%s%s%s',tmpdir_ref,fstem_ref,infix_ref,fext);
    p.fname_wm_T1 = sprintf('%s/%s%s_wm%s',tmpdir_ref,fstem_ref,infix_ref,fext);
    p.fname_bm_T1 = sprintf('%s/%s%s_bm%s',tmpdir_ref,fstem_ref,infix_ref,fext);
    p.fname_gm_T1 = sprintf('%s/%s%s_gm%s',tmpdir_ref,fstem_ref,infix_ref,fext);
    
    % skip if reference masks are missing
    if isempty(p.fname_wm_T1) || isempty(p.fname_bm_T1) || isempty(p.fname_gm_T1) ||...
       ~exist(p.fname_wm_T1,'file') || ~exist(p.fname_bm_T1,'file') || ~exist(p.fname_gm_T1,'file')
       fprintf('%s: missing reference masks for %s, skipping...\n',...
         mfilename,fname_in);
      return;
    end
  end

  if ~exist(p.fname_out,'file') || p.forceflag
    fprintf('%s: white matter bias correction %s...\n',mfilename,fname_in);
  end
  a = mmil_parms2args(p);
  fname_out = mmil_wm_corr(fname_in,a{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_ref_tmp,fname_ref_mask_tmp] = get_ref(fname_ref,parent_tmpdir)
  if ~exist('parent_tmpdir','var'), parent_tmpdir = []; end
  % identify reference input and mask
  [fpath_ref,fstem_ref,fext_ref] = cmig_get_fileparts(fname_ref);
  fstem_ref = get_raw_stem(fstem_ref);
  fname_ref_mask = sprintf('%s/%s.nii.gz',fpath_ref,get_mask_stem(fstem_ref));
  if ~isempty(parent_tmpdir)
    tmpdir_ref = sprintf('%s/%s',parent_tmpdir,fstem_ref);
  else
    tmpdir_ref = sprintf('%s/%s',fpath_ref,fstem_ref);
  end
  infix_ref = '_uw_wmbc_res'; fext = '.mgz';
  fname_ref_tmp = sprintf('%s/%s%s%s',tmpdir_ref,fstem_ref,infix_ref,fext);
  if ~exist(fname_ref_tmp,'file'), infix_ref = '_wmbc_res'; end
  fname_ref_tmp = sprintf('%s/%s%s%s',tmpdir_ref,fstem_ref,infix_ref,fext);
  fname_ref_mask_tmp = sprintf('%s/%s%s_mask%s',tmpdir_ref,fstem_ref,infix_ref,fext);
  % regenerate tmp files if necessary
  if ~exist(fname_ref_tmp,'file')
    fprintf('%s: regenerating temporary reference...\n',mfilename);
    fs_mri_convert(fname_ref,fname_ref_tmp);
  end
  if ~exist(fname_ref_mask_tmp,'file')
    fprintf('%s: regenerating temporary reference mask...\n',mfilename);
    fs_mri_convert(fname_ref_mask,fname_ref_mask_tmp);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_out,M_ref_to_raw,errcode,errmsg] = reg_to_ref(fname_in,fname_ref,scantype,parent_tmpdir,biascorr_applied,forceflag)
  fname_out = fname_in; M_ref_to_raw = [];
  errcode = 0; errmsg = [];
  if ~exist('parent_tmpdir','var'), parent_tmpdir = []; end
  if ~exist('biascorr_applied','var') || isempty(biascorr_applied), biascorr_applied = false; end
  if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end
  if (isempty(fname_ref) || ~exist(fname_ref,'file'))
    errmsg = sprintf('skipping registration for %s because no T1w reference',fname_in);
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end
  % create temporary output dir
  [fpath_in,fstem_in,fext] = fileparts(fname_in);
  tmpdir = sprintf('%s/reg_to_ref',fpath_in);
  mmil_mkdir(tmpdir);
  % identify reference input and mask
  [fname_ref,fname_ref_mask] = get_ref(fname_ref,parent_tmpdir);
  % register volumes
  fname_reg = sprintf('%s/registration_matrix.mat',tmpdir);
  if ~exist(fname_reg,'file') || forceflag
    fprintf('%s: registering %s to %s...\n',mfilename,fname_in,fname_ref);
    if strcmp(scantype,'T2w')
      M_ref_to_raw = mmil_reg_T2toT1(fname_in,fname_ref,...
        'outdir',tmpdir,'fname_T1_mask',fname_ref_mask,...
        'smoothmask_flag',true,...
        'smooth1',20,'thresh1',0.5,'smooth2',40,'thresh2',0.2,'smooth3',10,...
        'cleanup_flag',false,'forceflag',forceflag);
    else
      % smooth/dilate the reference mask before registration, as production
      %   did for additional same-type scans; same parameters as the
      %   smoothed mask in the T2w branch above
      fname_ref_mask_sm = sprintf('%s/ref_mask_sm.mgz',tmpdir);
      if ~exist(fname_ref_mask_sm,'file') || forceflag
        mmil_dilate_mask([],'fname_in',fname_ref_mask,...
          'fname_out',fname_ref_mask_sm,...
          'smooth1',20,'thresh1',0.5,'smooth2',40,'thresh2',0.2,'smooth3',10,...
          'forceflag',forceflag);
      end
      M_ref_to_raw = mmil_reg(fname_ref,fname_in,...
        'fname_maskA',fname_ref_mask_sm,'outdir',tmpdir,'forceflag',forceflag);
    end
    % save output file
    save(fname_reg,'M_ref_to_raw');
  else
    load(fname_reg);
  end

  [fpath_in,fstem_in,fext] = fileparts(fname_in);
  fname_out = sprintf('%s/%s_res%s',fpath_in,fstem_in,fext);
  if ~exist(fname_out,'file') || forceflag
    % resample into register with reference
    fprintf('%s: resampling %s to %s...\n',mfilename,fname_in,fname_out);
    [vol_in,M_in] = fs_load_mgh(fname_in);
    [M_ref,nvox_ref] = mmil_load_mgh_info(fname_ref,[],tmpdir);
    [vol_res,M_res] = mmil_resample_vol(vol_in,M_in,...
      'M_ref',M_ref,'nvox_ref',nvox_ref(1:3),'M_reg',M_ref_to_raw);
    % clip intensities only if wm bias correction normalized them
    if biascorr_applied
      image_range = get_image_range(scantype);
    else
      image_range = [];
    end
    if ~isempty(image_range)
      vol_res(vol_res<image_range(1)) = image_range(1);
      vol_res(vol_res>image_range(2)) = image_range(2);
    end
    % save output file
    fs_save_mgh(vol_res,fname_out,M_res);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function space_ref = get_ref_space(fname_ref,parms)
  % space of the reference volume, from its sidecar json
  % NOTE: a reference in its own space ('orig', or a sidecar that is missing
  %   or unreadable) is itself the space of a volume resampled into it, so
  %   the URI of the processed reference is returned instead of 'orig'
  [~,fstem_ref] = cmig_get_fileparts(fname_ref);
  space_ref = sprintf('bids::%s/%s/anat/%s.nii.gz',...
    parms.participant_id,parms.session_id,fstem_ref);
  fname_json = regexprep(fname_ref,'\.nii\.gz$','.json');
  if ~exist(fname_json,'file'), return; end
  try
    ref_info = jsondecode(fileread(fname_json));
  catch
    return;
  end
  if isfield(ref_info,'SpatialReference') && ~isempty(ref_info.SpatialReference) &&...
     ~strcmp(ref_info.SpatialReference,'orig')
    space_ref = ref_info.SpatialReference;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_out = export_data(fname_in,fstem,outdir)
  % set new name
  fname_out = sprintf('%s/%s.nii.gz',outdir,fstem);
  % convert from nifti to mgz
  fs_mri_convert(fname_in,fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;
  % set parent tmpdir so per-scan tmpdirs nest under it, enabling
  % selective cleanup (cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir))
  parms.tmpdir = sprintf('%s/proc_anat',parms.outdir);
  if parms.forceflag, procflag = 1; return; end
  % set procflag=1 if any output file does not exist
  for s=1:length(parms.fnames_data_t1)
    [~,fstem] = cmig_get_fileparts(parms.fnames_data_t1{s});
    fstem_proc = get_proc_stem(fstem);
    fname_out = sprintf('%s/%s.nii.gz',parms.outdir,fstem_proc);
    if ~exist(fname_out,'file'), procflag = 1; return; end
    fname_out_json = sprintf('%s/%s.json',parms.outdir,fstem_proc);
    if ~exist(fname_out_json,'file'), procflag = 1; return; end
    if s==1
      mask_stem = get_mask_stem(fstem);
      fname_out_mask = sprintf('%s/%s.nii.gz',parms.outdir,mask_stem);
      if ~exist(fname_out_mask,'file'), procflag = 1; return; end
      fname_out_mask_json = sprintf('%s/%s.json',parms.outdir,mask_stem);
      if ~exist(fname_out_mask_json,'file'), procflag = 1; return; end
    end
  end
  for s=1:length(parms.fnames_data_t2)
    [~,fstem] = cmig_get_fileparts(parms.fnames_data_t2{s});
    fstem_proc = get_proc_stem(fstem);
    fname_out = sprintf('%s/%s.nii.gz',parms.outdir,fstem_proc);
    if ~exist(fname_out,'file'), procflag = 1; return; end
    fname_out_json = sprintf('%s/%s.json',parms.outdir,fstem_proc);
    if ~exist(fname_out_json,'file'), procflag = 1; return; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



