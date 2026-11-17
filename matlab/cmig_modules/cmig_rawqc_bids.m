function [outdir,errcode,errmsg] = cmig_rawqc_bids(participant_id,session_id,input_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_rawqc_bids(participant_id,session_id,input_rootdir,output_rootdir,[options])
%
% Purpose: Perform automated quality control procedures on BIDS rawdata
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
%   'modality': modality subdirectory string
%     allowed values: {'anat','dwi','func','fmap'}
%     {default = 'anat'}
%   'cleanup_flag': [0|1] remove temporary files
%     {default = 1}
%   'export_flag': [0|1] create summary images for manual review
%     {default = 1}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Modality-specific Input ('key',value):
%   'motion_radius': head radius (mm) for rotation-to-distance conversion
%     for dwi and func data
%     {default = 50}
%   'method': motion correction method for func data
%     allowed: 'fslmean','fsl','afni'
%     {default = 'fslmean'}
%   'fMRI_minTRs': minimum frames of a func scan for its multi-frame metrics
%     (brain_tSNR, motion, fwhm) and motion plots; these skip the first
%     skipTRs (8) frames, so keep it at least skipTRs + 12
%     {default = 20}
%   'afni_flag': [0|1] for dwi data, use AFNI's 3dvolreg for within-scan
%     registration instead of rbreg_vol2vol_icorr_mask
%     {default = 1}
%   'censor_reg_flag': [0|1] for dwi motion correction, exclude censored
%     slices from the per-frame registration and replace them after
%     alignment, instead of replacing them beforehand
%     {default = 1}
%
% NOTE: afni_flag and censor_reg_flag default ON as of 08/21/2026, so
%   dti_correct_data must be one that accepts them. mmil_args2parms is strict,
%   so they are still named only when enabled, but with both on by default that
%   no longer buys compatibility with an MMPS dti_correct_data that predates
%   them; set both to 0 to run against one.
%   insensitive to the registration mask, and censor_reg_flag writes its own
%   per-frame weight for censored frames regardless.
%
% Output:
%   outdir: output directory for the specified modality
%     form: {output_rootdir}/{participant_id}/{session_id}/{modality}
%   errcode: error code (0 = success, 1 = error)
%     a per-scan failure does not set errcode if at least one scan succeeds
%   errmsg: error message string (empty if no error)
%     accumulates per-scan failure messages, which are also recorded in
%     {participant_id}_{session_id}_desc-rawqc_errors.tsv so that reruns
%     retry the failed scans and keep the failure record in the status file
%
% Created:  02/24/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 06/04/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 06/07/2026 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 06/26/2026 by Don Hagler
% Prev Mod: 06/28/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/22/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/12/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/20/2026 by Don Hagler
% Prev Mod: 08/21/2026 by Don Hagler
% Prev Mod: 08/25/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/22/2026 by Don Hagler
% Prev Mod: 09/24/2026 by Don Hagler
% Prev Mod: 09/26/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/20/2026 by Don Hagler
% Prev Mod: 10/23/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Last Mod: 11/15/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on MMIL_autoQC_Exam
%   Created:  02/20/16 by Don Hagler
%   Last Mod: 12/05/25 by Don Hagler
% and cmig_proc_bids_anat
%   Created:  08/04/2025 by Don Hagler
%   Last Mod: 09/04/2025 by Don Hagler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,4), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');
parms = check_parms(participant_id,session_id,input_rootdir,output_rootdir,varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[parms,errcode,errmsg] = check_input(parms);
if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-rawqc',errcode,errmsg,[],startdate); return; end

procflag = check_output(parms);
if ~procflag
  outdir = parms.outdir;
  % keep per-scan failure record from a previous run in the status file
  errmsg = read_scan_errors(parms);
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-rawqc',errcode,errmsg,[],startdate,procflag);
  return;
end

if parms.verbose
  fprintf('%s: processing modality %s\n',mfilename,parms.modality);
  for f=1:length(parms.fnames_data)
    fprintf('%s: input file: %s\n',mfilename,parms.fnames_data{f});
  end
  fprintf('%s: output directory: %s\n',mfilename,parms.outdir);
end

% get unique BIDS series types within this modality
stypes = get_bids_stypes(parms.fnames_data);

all_qcstruct = {};
failed_scans = {};
for s=1:length(stypes)
  stype = stypes{s};

  % get file list and sidecar list for this series type
  idx = get_stype_indices(parms.fnames_data,stype);
  fnamelist = parms.fnames_data(idx);

  % set parms for this modality and series type
  tparms = parms;
  tparms.fnames_info = parms.fnames_info(idx);

  if tparms.verbose
    fprintf('%s: processing series type %s (%d files)\n',...
      mfilename,stype,length(fnamelist));
  end

  stype_qcstruct = {};
  for f=1:length(fnamelist)
    fname = fnamelist{f};
    fstem = get_fstem(fname);
    auto_qcinfo = [];
    scan_errmsg = [];
    fname_vol_mat = sprintf('%s/%s_vol.mat',tparms.tmpdir,fstem);
    fname_qcinfo_mat = sprintf('%s/%s_qcinfo.mat',tparms.tmpdir,fstem);
    if ~exist(fname_vol_mat,'file') || ~exist(fname_qcinfo_mat,'file') || tparms.forceflag
      try
        mmil_mkdir(tparms.tmpdir);
        auto_qcinfo = check_scan(fname,stype,tparms);
        if ~isempty(auto_qcinfo)
          save(fname_qcinfo_mat,'auto_qcinfo');
        end
      catch me
        scan_errmsg = sprintf('%s: %s',fstem,me.message);
        fprintf('%s: ERROR: %s\n',mfilename,scan_errmsg);
      end
    else
      tmp = load(fname_qcinfo_mat,'auto_qcinfo');
      auto_qcinfo = tmp.auto_qcinfo;
    end
    if ~isempty(auto_qcinfo) && tparms.export_flag
      export_image(fname,stype,tparms,auto_qcinfo);
    end
    if ~isempty(auto_qcinfo)
      qcs = convert_qcinfo(auto_qcinfo,tparms.thresholds,parms.participant_id,parms.session_id);
      all_qcstruct{end+1} = qcs;
      stype_qcstruct{end+1} = qcs;
    else
      fprintf('%s: WARNING: %s auto_qcinfo is empty\n',mfilename,fstem);
      if isempty(scan_errmsg)
        scan_errmsg = sprintf('%s: auto_qcinfo is empty',fstem);
      end
      failed_scans{end+1} = struct('fstem',fstem,'errmsg',scan_errmsg);
    end
  end

  % write series-level TSV for this stype
  if ~isempty(stype_qcstruct)
    if strcmp(parms.modality,'fmap')
      acq_vals = unique(cellfun(@(r) r.acq, stype_qcstruct, 'UniformOutput', false));
      for a=1:length(acq_vals)
        acq = acq_vals{a};
        acq_struct = stype_qcstruct(cellfun(@(r) strcmp(r.acq,acq), stype_qcstruct));
        fname_acq_tsv = sprintf('%s/%s_%s_acq-%s_desc-epi_rawqc.tsv',...
          parms.outdir,parms.participant_id,parms.session_id,acq);
        write_series_rawqc(acq_struct,stype,fname_acq_tsv,parms);
      end
    else
      desc = stype_to_desc(stype);
      fname_stype_tsv = sprintf('%s/%s_%s_desc-%s_rawqc.tsv',...
        parms.outdir,parms.participant_id,parms.session_id,desc);
      write_series_rawqc(stype_qcstruct,stype,fname_stype_tsv,parms);
    end
  end
end

% write session-level TSV
if ~isempty(all_qcstruct)
  fname_session_tsv = sprintf('%s/%s_%s_desc-session_rawqc.tsv',...
    parms.outdir,parms.participant_id,parms.session_id);
  write_session_rawqc(all_qcstruct,fname_session_tsv,parms);
end

% record per-scan failures; job fails outright only if no scan succeeded
[errcode,errmsg] = record_scan_errors(failed_scans,all_qcstruct,parms);

if parms.cleanup_flag
  cleanup_errmsg = cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir);
  if isempty(errmsg), errmsg = cleanup_errmsg; end
end

cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-rawqc',errcode,errmsg,[],startdate,procflag);

outdir = parms.outdir;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_parms(participant_id,session_id,input_rootdir,output_rootdir,options)
  options = cmig_parse_deployed_args(options);
  parms = mmil_args2parms(options,{...
    'modality','anat',{'anat','dwi','func','fmap'},...
    'cleanup_flag',true,[false true],...
  ...
    'export_flag',true,[false true],...
    'mgz_flag',false,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  ...
    'motion_radius',50,[1,100],...
    'method','fslmean',{'fslmean','fsl','afni'},...
    'afni_flag',true,[false true],...
    'censor_reg_flag',true,[false true],...
  ... % hidden
    'fMRI_minTRs',20,[2,1000],...
    'dMRI_minTRs',10,[2,1000],...
    'dti_corr_suffix','ecc_mc',[],...
    'censor_min_fract',0.8,[0,1],...
    'DTfit_nonlin_flag',false,[false true],...
    'DTI_conv_meas_list',{'FA','MD','b0','V0'},[],...
    'DTI_meas_list',{'FA','MD','b0'},[],...
    'DTI_err_list',{'rms_err'},{'err','rms_err','rms'},...
    'BOLD_meas_list',{'mean','std','tsnr'},[],...
    'bb_buffer',20,[0,100],...
    'nii_flag',false,[false true],...
    'thresholds',[0.2 0.3 0.4],[],...
    'dMRI_censor_thresh',3,[],... % fixed, as HBCD raw QC; [] = per-scan percentile
    'dMRI_censor_prctile',99,[0,100],...
  ... % hidden for subplots
    'planestrings',{'HOR','SAG','COR'},[],...
    'multislice_nrows',[4],[],... % may be a vector
    'multislice_min',0.25,[],...
    'multislice_max',0.75,[],...
    'multislice_plim',[1,99],[],...
    'multislice_fig_size',[8 10],[],...
    'multiframe_plim',[5,99.9],[],...
    'multiframe_slice_fract',0.5,[],...
    'multiframe_frame_fracts',linspace(0,1,81),[],...
    'multiframe_fig_size',[10 10],[],...
    'threeviews_plim',[1,99],[],...
    'threeviews_slice_fract',[0.5,0.6,0.5],[],...
    'threeviews_fig_size',[10 4],[],...
    'motion_fig_size',[10 8],[],...
    'tif_dpi',300,[],...
  ... % hidden for motion filtering
    'resp_low',18.6,[],...
    'resp_high',25.7,[],...
  ... % hidden for fwhm and tSNR
    'skipTRs',8,[0 Inf],...
  ...
    'mri_convert_tags',{'out_orient','options','verbose','forceflag'},[],...
    'brainmask_tags',{'outdir','fname_mask','binarize_flag','verbose','forceflag'},[],...
    'tSNR_tags',{'outdir','mcflag','hp_sigma','TR','skipTRs','method','snr_suffix',...
                 'verbose','forceflag','ext'},[],...
    'dti_corr_tags',{'B0unwarp_flag','ecc_flag','censor_flag',...
                     'motion_B0uw_flag','mc_flag','gradunwarp_flag',...
                     'fname_out','fname_qmat','qmat','bvals','fname_censor',...
                     'censor_min_ndirs','censor_thresh','nonlin_flag',...
                     'b0_thresh','interpm','maskoutput','verbose','forceflag',...
                     'require_B0uw_flag','optimize_B0uw_flag',...
                     'inorm_B0uw_flag','fname_B0dx','fname_for','fname_rev',...
                     'fname_B0uw_ref','fname_B0uw_mask','fname_B0uw_reg',...
                     'PhaseDir','revflag','driftcorr',...
                     'motion_B0uw_iters','min_trans',...
                     'min_rot','fname_motion_B0uw',...
                     'gruw_type','gruw_unwarpflag','gruw_isoctrflag',...
                     'fname_ref','fname_mask','fname_reg',...
                     'mstep','scales','ecc_censor_niter','mc_censor_niter',...
                     'gruw_jacobian_flag','kernelWidthMax','lambda2',...
                     'kernelWidthMax_vec','lambda2_vec','multi_opt_flag',...
                     'outfix','outext','tmpext'},[],...
    'dMRI_shells_tags',{'data','infix','outdir','outstem','tmpdir','meas_flag',...
                        'censor_thresh','censor_prctile',...
                        'fname_DTmeas','fname_DTerr',...
                        'cleanupflag','verbose','colorbar_flag','bids_flag','forceflag'},[],...
    'convert_DTmeas_tags',{'regT1flag','zpadflag','M','RegInfo',...
                           'measlist','mgz_flag','dat_flag','dat_measlist',...
                           'nii_flag','nii_measlist','nii_out_orient',...
                           'verbose_flag','forceflag','dat_revsliceflag',...
                           'dat_permvec'},[],...
  });
    
  parms.participant_id = participant_id;
  parms.session_id = session_id;
  parms.input_rootdir = input_rootdir;
  parms.output_rootdir = output_rootdir;

  parms.outdir = sprintf('%s/%s/%s/%s',...
    parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);
  parms.tmpdir = sprintf('%s/autoQC',parms.outdir);
  % NOTE: extension of the intermediate volumes, from mgz_flag. Uncompressed by
  %   default: gzip costs time on every write and every read and saves nothing
  %   on files cleanup_flag deletes. Set mgz_flag with cleanup_flag = 0, where
  %   the intermediates are kept deliberately and the footprint matters
  if parms.mgz_flag
    parms.in_ext = '.mgz';
  else
    parms.in_ext = '.mgh';
  end
  if ~isempty(parms.motion_fig_size) && length(parms.motion_fig_size)==1
    parms.motion_fig_size = [parms.motion_fig_size parms.motion_fig_size];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];
  % NOTE: ignore_scans_flag = 1, the last argument, so that scans marked unusable
  %   in scans.tsv are still returned; raw QC is what determines usability, and
  %   the empty series_type, task and acq_type are only there to reach it
  [parms.fnames_data,parms.fnames_info] = cmig_get_bids_data(...
    parms.participant_id,parms.session_id,parms.input_rootdir,...
    parms.modality,[],[],[],1);
  if isempty(parms.fnames_data)
    errcode = 1;
    errmsg = 'not processed: no data found';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function procflag = check_output(parms)
  procflag = 1;
  mmil_mkdir(parms.outdir);
  if parms.forceflag, return; end
  fname_tsv = sprintf('%s/%s_%s_desc-session_rawqc.tsv',...
    parms.outdir,parms.participant_id,parms.session_id);
  if ~exist(fname_tsv,'file'), return; end
  % reprocess if any expected scan has no row in its series-level TSV
  % (e.g., a scan that failed on a previous run)
  for f=1:length(parms.fnames_data)
    fname = parms.fnames_data{f};
    fstem = get_fstem(fname);
    stype = get_bids_stype(fname);
    if strcmp(parms.modality,'fmap')
      fname_series_tsv = sprintf('%s/%s_%s_acq-%s_desc-epi_rawqc.tsv',...
        parms.outdir,parms.participant_id,parms.session_id,get_acq_entity(fstem));
    else
      fname_series_tsv = sprintf('%s/%s_%s_desc-%s_rawqc.tsv',...
        parms.outdir,parms.participant_id,parms.session_id,stype_to_desc(stype));
    end
    if ~exist(fname_series_tsv,'file'), return; end
    if ~tsv_contains_fstem(fname_series_tsv,fstem), return; end
  end
  procflag = 0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function found = tsv_contains_fstem(fname_tsv,fstem)
  % check whether any row of a series-level TSV refers to this scan
  found = false;
  try
    raw = mmil_readtext(fname_tsv,'\t');
  catch
    return;
  end
  if size(raw,1) < 2, return; end
  j = find(strcmp(raw(1,:),'filename'),1);
  if isempty(j), found = true; return; end % no filename column: cannot check
  for i=2:size(raw,1)
    val = raw{i,j};
    if ischar(val) && ~isempty(strfind(val,fstem))
      found = true;
      return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_errors = get_errors_fname(parms)
  fname_errors = sprintf('%s/%s_%s_desc-rawqc_errors.tsv',...
    parms.outdir,parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = record_scan_errors(failed_scans,all_qcstruct,parms)
  % write per-scan failures to a persistent errors file next to the status
  % file so the record survives reruns; remove the file when all scans succeed
  errcode = 0; errmsg = [];
  fname_errors = get_errors_fname(parms);
  if isempty(failed_scans)
    if exist(fname_errors,'file'), delete(fname_errors); end
    return;
  end
  msgs = cellfun(@(x) x.errmsg,failed_scans,'UniformOutput',false);
  errmsg = strjoin(msgs,'; ');
  if isempty(all_qcstruct), errcode = 1; end
  rows = {};
  for i=1:length(failed_scans)
    row = [];
    row.participant_id = parms.participant_id;
    row.session_id = parms.session_id;
    row.fstem = failed_scans{i}.fstem;
    row.errmsg = failed_scans{i}.errmsg;
    rows{end+1} = row;
  end
  write_rawqc_tsv(rows,fname_errors,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = read_scan_errors(parms)
  % read per-scan failure messages recorded by a previous run
  errmsg = [];
  fname_errors = get_errors_fname(parms);
  if ~exist(fname_errors,'file'), return; end
  try
    raw = mmil_readtext(fname_errors,'\t');
  catch
    return;
  end
  if size(raw,1) < 2, return; end
  j = find(strcmp(raw(1,:),'errmsg'),1);
  if isempty(j), return; end
  msgs = raw(2:end,j);
  msgs = msgs(cellfun(@ischar,msgs));
  errmsg = strjoin(msgs,'; ');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stype = get_bids_stype(fname)
  % extract BIDS suffix (series type) from a BIDS file name
  fname_tmp = regexprep(fname,'\.nii.*','');
  [~,fstem] = fileparts(fname_tmp);
  k = regexp(fstem,'_([^_]+)$','tokens');
  if ~isempty(k)
    stype = k{1}{1};
  else
    stype = fstem;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stypes = get_bids_stypes(fnames_data)
  % get unique BIDS series types from a list of file names
  stypes = {};
  for f=1:length(fnames_data)
    stype = get_bids_stype(fnames_data{f});
    if ~ismember(stype,stypes)
      stypes{end+1} = stype;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idx = get_stype_indices(fnames_data,stype)
  % find indices of files matching a given BIDS series type
  idx = [];
  for f=1:length(fnames_data)
    if strcmp(get_bids_stype(fnames_data{f}),stype)
      idx(end+1) = f;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_tmp,fstem] = import_data(fname_data,tmpdir,in_ext)
  [fpath,fstem] = cmig_get_fileparts(fname_data);
  fname_tmp = sprintf('%s/%s%s',tmpdir,fstem,in_ext);
  if ~exist(fname_tmp,'file')
    fs_mri_convert(fname_data,fname_tmp,'options','-odt float');
  end
  % copy sidecar files so tmpdir MGZ path works with get_scaninfo/get_seriesinfo
  for ext = {'.json','.bvec','.bval'}
    src = sprintf('%s/%s%s',fpath,fstem,ext{1});
    dst = sprintf('%s/%s%s',tmpdir,fstem,ext{1});
    if exist(src,'file') && ~exist(dst,'file')
      copyfile(src,dst);
      fileattrib(dst,'+w','u');
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_scan(fname,stype,parms)
  auto_qcinfo = [];
  % dispatch to check function based on BIDS modality and series type
  if strcmp(parms.modality,'func') || strcmp(stype,'bold')
    auto_qcinfo = check_fMRI(fname,stype,parms);
  elseif strcmp(parms.modality,'dwi') || strcmp(stype,'dwi')
    auto_qcinfo = check_dMRI(fname,stype,parms);
  elseif strcmp(parms.modality,'fmap') || strcmp(stype,'epi')
    auto_qcinfo = check_fmap(fname,stype,parms);
  elseif strcmp(parms.modality,'anat') || ...
         ismember(stype,{'T1w','T2w','T2starw','FLAIR','PDw',...
                         'angio','VFA','MP2RAGE','FLASH',...
                         'TB1map','TB1EPI','TB1AFI','TB1TFL','RB1COR'})
    auto_qcinfo = check_sMRI(fname,stype,parms);
  else
    fprintf('%s: no quality checks implemented for modality %s, type %s\n',...
      mfilename,parms.modality,stype);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_image(fname,stype,parms,auto_qcinfo)
  % dispatch to export function based on BIDS modality and series type
  if strcmp(parms.modality,'func') || strcmp(stype,'bold')
    export_fMRI(fname,stype,parms);
  elseif strcmp(parms.modality,'dwi') || strcmp(stype,'dwi')
    export_dMRI(fname,stype,parms,auto_qcinfo);
  elseif strcmp(parms.modality,'fmap') || strcmp(stype,'epi')
    export_fmap(fname,stype,parms);
  elseif strcmp(parms.modality,'anat') || ...
         ismember(stype,{'T1w','T2w','T2starw','FLAIR','PDw',...
                         'angio','VFA','MP2RAGE','FLASH',...
                         'TB1map','TB1EPI','TB1AFI','TB1TFL','RB1COR'})
    export_sMRI(fname,stype,parms);
  else
    fprintf('%s: no image export implemented for modality %s, type %s\n',...
      mfilename,parms.modality,stype);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = init_auto_qcinfo(fname,fname_mgh,stype,parms)
  auto_qcinfo.stype = stype;
  auto_qcinfo.fstem = get_fstem(fname);
  auto_qcinfo.scaninfo = get_scaninfo(fname,stype,parms);
  auto_qcinfo.seriesinfo = get_seriesinfo(fname,parms);
  auto_qcinfo.imginfo = get_imginfo(fname_mgh,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_sMRI(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking sMRI: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  auto_qcinfo = init_auto_qcinfo(fname,fname_mgh,stype,parms);
  auto_qcinfo.fwhm = [];
  fstem = get_fstem(fname_mgh);
  fname_brainmask = create_brainmask(fname_mgh,parms);
  if isempty(fname_brainmask), return; end
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms,fstem);
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname_mgh,fname_brainmask,parms);
  auto_qcinfo.nonbrain_SNR = calc_nonbrain_SNR(fname_mgh,fname_brainmask,parms);
  % spatial smoothness of the image as imported, no detrending: a structural
  % scan is one frame, so there is nothing to detrend, as for field maps.
  % Expected to carry motion-blur signal, which ABCD has no automated
  % structural metric for, but that association is not yet validated
  % NOTE: fname, the original nii.gz, not the imported copy. cmig_3dFWHM needs
  %   nifti and would otherwise convert the import straight back, a second
  %   mri_convert and a duplicate file per series. Validated as identical to
  %   the four decimals 3dFWHMx prints
  auto_qcinfo.fwhm = calc_fwhm_raw(fname,fname_brainmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_fmap(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking fmap: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  auto_qcinfo = init_auto_qcinfo(fname,fname_mgh,stype,parms);
  auto_qcinfo.fwhm = [];
  fstem = get_fstem(fname_mgh);
  [src_entities,src_suffix] = split_bids_fstem(fstem);
  fname_brainmask = create_brainmask(fname_mgh,parms);
  if isempty(fname_brainmask), return; end
  [auto_qcinfo.line_artifact,auto_qcinfo.line2_artifact] = ...
    detect_lines(fname_mgh,parms,fname_brainmask);
  auto_qcinfo.FOV_cutoff = detect_FOV_cutoff(fname_brainmask,parms,src_entities,src_suffix);
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms,fstem);
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname_mgh,fname_brainmask,parms);
  % NOTE: fname, not fname_mgh, as in check_sMRI and check_fMRI
  auto_qcinfo.fwhm = calc_fwhm_raw(fname,fname_brainmask,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_fMRI(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking fMRI: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  auto_qcinfo = init_auto_qcinfo(fname,fname_mgh,stype,parms);
  fstem = get_fstem(fname_mgh);
  [src_entities,src_suffix] = split_bids_fstem(fstem);
  [fname_brainmask,auto_qcinfo] = calc_fmri_multiframe_metrics(fname_mgh,auto_qcinfo,parms);
  if isempty(fname_brainmask), return; end
  [auto_qcinfo.line_artifact,auto_qcinfo.line2_artifact] = ...
    detect_lines(fname_mgh,parms,fname_brainmask);
  auto_qcinfo.FOV_cutoff = detect_FOV_cutoff(fname_brainmask,parms,src_entities,src_suffix);
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms,fstem);
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname_mgh,fname_brainmask,parms);
  auto_qcinfo = check_pepolar(fname,fname_mgh,fname_brainmask,auto_qcinfo,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_dMRI(fname,stype,parms)
  auto_qcinfo = [];
  if parms.verbose
    fprintf('%s: checking dMRI: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  auto_qcinfo = init_auto_qcinfo(fname,fname_mgh,stype,parms);
  fstem = get_fstem(fname_mgh);
  [src_entities,src_suffix] = split_bids_fstem(fstem);
  fname_brainmask = create_brainmask(fname_mgh,parms);
  if isempty(fname_brainmask), return; end
  [auto_qcinfo.line_artifact,auto_qcinfo.line2_artifact] = ...
    detect_lines(fname_mgh,parms,fname_brainmask);
  auto_qcinfo.FOV_cutoff = detect_FOV_cutoff(fname_brainmask,parms,src_entities,src_suffix);
  auto_qcinfo.brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms,fstem);
  auto_qcinfo.brain_SNR = calc_brain_SNR(fname_mgh,fname_brainmask,parms);
  auto_qcinfo = check_pepolar(fname,fname_mgh,fname_brainmask,auto_qcinfo,parms);
  old_ndirs = size(auto_qcinfo.scaninfo.qmat,1);
  [auto_qcinfo.scaninfo.qmat,errcode] = ...
    check_qmat(auto_qcinfo.scaninfo.qmat,auto_qcinfo.imginfo.nreps);
  if errcode, return; end
  % keep bvals consistent with trimmed qmat
  new_ndirs = size(auto_qcinfo.scaninfo.qmat,1);
  if new_ndirs < old_ndirs && ~isempty(auto_qcinfo.scaninfo.bvals)
    auto_qcinfo.scaninfo.bvals = auto_qcinfo.scaninfo.bvals(1:new_ndirs);
  end
  auto_qcinfo = check_dMRI_fit(fname_mgh,fname_brainmask,auto_qcinfo,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qmat,errcode] = check_qmat(qmat,nframes)
  errcode = 0;
  ndirs = size(qmat,1);
  if nframes == ndirs-1
    fprintf('%s: WARNING: one fewer frame (%d) than directions in qmat (%d)\n',...
      mfilename,nframes,ndirs);
    qmat = qmat(1:nframes,:);
  elseif nframes ~= ndirs
    fprintf('%s: WARNING: mismatch between frames (%d) and directions in qmat (%d)\n',...
      mfilename,nframes,ndirs);
    errcode = 1;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_for,fname_rev] = find_pepolar_fmaps(fname,parms)
  fname_for = []; fname_rev = [];
  % NOTE: ignore_scans_flag = 1 here too, or a scan whose field maps are marked
  %   unusable would lose its pepolar correction during raw QC
  [fnames_fmap_data,fnames_fmap_info] = cmig_get_bids_data(...
    parms.participant_id,parms.session_id,parms.input_rootdir,...
    'fmap','epi',[],[],1);
  if isempty(fnames_fmap_data), return; end
  [fpath,fstem] = cmig_get_fileparts(fname);
  fname_info = sprintf('%s/%s.json',fpath,fstem);
  if ~exist(fname_info,'file')
    fprintf('%s: WARNING: no JSON sidecar for %s\n',mfilename,get_fstem(fname));
    return;
  end
  [fname_for,fname_rev] = cmig_select_bids_fmap(fname_info,fnames_fmap_data,fnames_fmap_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_pepolar(fname,fname_mgh,fname_brainmask,auto_qcinfo,parms)
  auto_qcinfo.pepolar.mi_diff = NaN;
  auto_qcinfo.pepolar.mi_for  = NaN;
  auto_qcinfo.pepolar.mi_rev  = NaN;
  [fname_for,fname_rev] = find_pepolar_fmaps(fname,parms);
  if isempty(fname_for) || isempty(fname_rev)
    fprintf('%s: WARNING: pepolar check skipped — no fmap pair found for %s\n',...
      mfilename,get_fstem(fname));
    return;
  end
  fname_for_mgh = import_data(fname_for,parms.tmpdir,parms.in_ext);
  fname_rev_mgh = import_data(fname_rev,parms.tmpdir,parms.in_ext);
  fstem = get_fstem(fname_mgh);
  fname_cache = sprintf('%s/%s_pepolar.mat',parms.tmpdir,fstem);
  [entities_o,suffix_o] = split_bids_fstem(fstem);
  fname_plot = sprintf('%s/%s_desc-pepolar_%s.jpg',parms.outdir,entities_o,suffix_o);
  if exist(fname_cache,'file') && ~parms.forceflag
    tmp = load(fname_cache,'pepolar');
    auto_qcinfo.pepolar = tmp.pepolar;
    return;
  end
  try
    [mi_diff,mi_for,mi_rev] = mmil_check_pepolar(fname_mgh,fname_for_mgh,fname_rev_mgh,...
      'fname_mask',fname_brainmask,'fname_plot',fname_plot,'forceflag',parms.forceflag);
    auto_qcinfo.pepolar.mi_diff = mi_diff;
    auto_qcinfo.pepolar.mi_for  = mi_for;
    auto_qcinfo.pepolar.mi_rev  = mi_rev;
  catch me
    fprintf('%s: WARNING: pepolar check failed for %s: %s\n',...
      mfilename,get_fstem(fname),me.message);
  end
  pepolar = auto_qcinfo.pepolar;
  save(fname_cache,'pepolar');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_sMRI(fname,stype,parms)
  if parms.verbose
    fprintf('%s: exporting sMRI %s: %s...\n',mfilename,stype,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  plot_multislice(fname_mgh,parms);
  plot_threeviews(fname_mgh,parms);
  plot_threeviews_mask(fname_mgh,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_fmap(fname,stype,parms)
  if parms.verbose
    fprintf('%s: exporting fmap: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  plot_multislice(fname_mgh,parms);
  plot_threeviews(fname_mgh,parms);
  plot_threeviews_mask(fname_mgh,parms);
  imginfo = get_imginfo(fname_mgh,parms);
  if imginfo.nreps > 1
    plot_multiframe(fname_mgh,parms);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_dMRI(fname,stype,parms,auto_qcinfo)
  if parms.verbose
    fprintf('%s: exporting dMRI: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  plot_multislice(fname_mgh,parms);
  plot_multiframe(fname_mgh,parms,stype,0);
  plot_multiframe(fname_mgh,parms,stype,1);
  plot_threeviews(fname_mgh,parms);
  plot_threeviews_mask(fname_mgh,parms);
  plot_multislice_measures(fname_mgh,parms,parms.DTI_meas_list,parms.dti_corr_suffix);
  plot_threeviews_measures(fname_mgh,parms,parms.DTI_meas_list,parms.dti_corr_suffix);
  plot_threeviews_rgb(fname_mgh,parms,parms.dti_corr_suffix);
  plot_multislice_measures(fname_mgh,parms,parms.DTI_err_list,parms.dti_corr_suffix);
  plot_multiframe_measures(fname_mgh,stype,parms,parms.DTI_err_list,parms.dti_corr_suffix);
  plot_threeviews_measures(fname_mgh,parms,parms.DTI_err_list,parms.dti_corr_suffix);
  % reuse the motion estimate check_dMRI_fit already made, rather than repeating
  %   the whole eddy and motion correction. With forceflag = 0 the repeat was a
  %   no-op, since dti_correct_data found its own output; with forceflag = 1 it
  %   ran the correction a second time, which is where the 2.2x cost came from.
  %   mmil_getfield rather than a direct read: the field is empty when the scan
  %   has too few frames to estimate motion from
  motion_data = mmil_getfield(auto_qcinfo,'motion',[]);
  if ~isempty(motion_data)
    fstem = get_fstem(fname_mgh);
    plot_motion(fstem,motion_data,parms);
    plot_motion_fd(fstem,motion_data,parms);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_fMRI(fname,stype,parms)
  if parms.verbose
    fprintf('%s: exporting fMRI: %s...\n',mfilename,fname);
  end
  fname_mgh = import_data(fname,parms.tmpdir,parms.in_ext);
  plot_multislice(fname_mgh,parms);
  plot_multiframe(fname_mgh,parms);
  plot_threeviews(fname_mgh,parms);
  plot_threeviews_mask(fname_mgh,parms);
  plot_multislice_measures(fname_mgh,parms,parms.BOLD_meas_list);
  plot_threeviews_measures(fname_mgh,parms,parms.BOLD_meas_list);
  imginfo = get_imginfo(fname_mgh,parms);
  if imginfo.nreps>=parms.fMRI_minTRs
    fstem = get_fstem(fname_mgh);
    scaninfo = get_scaninfo(fname,stype,parms);
    motion_data = calc_fmri_motion(fname_mgh,scaninfo,parms);
    plot_motion(fstem,motion_data,parms);
    plot_motion_fd(fstem,motion_data,parms);
    plot_motion(fstem,motion_data,parms,1);
    plot_motion_fd(fstem,motion_data,parms,1);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multislice(fname,parms,frame,entities_o,suffix_o)
  if ~exist('frame','var'), frame = []; end
  if ~exist('entities_o','var') || isempty(entities_o)
    [entities_o,suffix_o] = split_bids_fstem(get_fstem(fname));
  end
  tparms = [];
  tparms.frame = frame;
  tparms.outdir = parms.outdir;
  tparms.plim = parms.multislice_plim;
  tparms.fig_size = parms.multislice_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.forceflag = parms.forceflag;
  for j=1:length(parms.multislice_nrows)
    nrows = parms.multislice_nrows(j);
    tparms.slice_fracts = linspace(parms.multislice_min,parms.multislice_max,nrows.^2);
    for i=1:length(parms.planestrings)
      tparms.planestring = parms.planestrings{i};
      plane = lower(tparms.planestring(1:3));  % 'AX','SAG','COR'
      if ~isempty(frame)
        tparms.outstem = sprintf('%s_desc-f%d%sMS_%s',...
          entities_o,frame,plane,suffix_o);
      else
        tparms.outstem = sprintf('%s_desc-%sMS_%s',...
          entities_o,plane,suffix_o);
      end
      fname_out = sprintf('%s/%s.jpg',tparms.outdir,tparms.outstem);
      if ~exist(fname_out,'file') || parms.forceflag
        args = mmil_parms2args(tparms);
        mmil_multislice_subplots(fname,args{:});
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multiframe(fname,parms,stype,groupflag,entities_o,suffix_o)
  if ~exist('stype','var'), stype = []; end
  if ~exist('groupflag','var') || isempty(groupflag), groupflag = 0; end
  if ~exist('entities_o','var') || isempty(entities_o)
    [entities_o,suffix_o] = split_bids_fstem(get_fstem(fname));
  end
  mf_suffix = 'MF';
  tparms = [];
  tparms.outdir = parms.outdir;
  tparms.slice_fract = parms.multiframe_slice_fract;
  tparms.frame_fracts = parms.multiframe_frame_fracts;
  tparms.plim = parms.multiframe_plim;
  tparms.fig_size = parms.multiframe_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.forceflag = parms.forceflag;
  [M,volsz] = mmil_load_mgh_info(fname,0,parms.tmpdir);
  nf = volsz(4);
  if nf==1, return; end
  if strcmp(stype,'dwi')
    tparms.frame_fracts = [];
    if groupflag
      scaninfo = get_scaninfo(fname,stype,parms);
      if isfield(scaninfo,'bvals') && ~isempty(scaninfo.bvals)
        tparms.groups = round(scaninfo.bvals);
      else
        qmat = scaninfo.qmat;
        qlength = max(sqrt(sum(qmat.^2,2)),eps);
        tparms.groups = round(scaninfo.bval.*(qlength.^2));
      end
      mf_suffix = 'MFb';
      if nf~=length(tparms.groups), return; end
    end
  end
  for i=1:length(parms.planestrings)
    tparms.planestring = parms.planestrings{i};
    plane = lower(tparms.planestring(1:3));
    tparms.outstem = sprintf('%s_desc-%s%s_%s',entities_o,plane,mf_suffix,suffix_o);
    fname_out = sprintf('%s/%s.jpg',tparms.outdir,tparms.outstem);
    if ~exist(fname_out,'file') || parms.forceflag
      args = mmil_parms2args(tparms);
      mmil_multiframe_subplots(fname,args{:});
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews(fname,parms,frame,entities_o,suffix_o,space_o)
  if ~exist('frame','var'), frame = []; end
  if ~exist('entities_o','var') || isempty(entities_o)
    [entities_o,suffix_o] = split_bids_fstem(get_fstem(fname));
    space_o = '';
  end
  if ~exist('space_o','var') || isempty(space_o), space_o = ''; end
  if ~isempty(space_o)
    base = sprintf('%s_space-%s',entities_o,space_o);
  else
    base = entities_o;
  end
  tparms = [];
  tparms.frame = frame;
  tparms.outdir = parms.outdir;
  if ~isempty(frame)
    tparms.outstem = sprintf('%s_desc-f%d3view_%s',base,frame,suffix_o);
  else
    tparms.outstem = sprintf('%s_desc-3view_%s',base,suffix_o);
  end
  tparms.plim = parms.threeviews_plim;
  tparms.fig_size = parms.threeviews_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.planestrings = parms.planestrings;
  tparms.slice_fract = parms.threeviews_slice_fract;
  tparms.forceflag = parms.forceflag;
  fname_out = sprintf('%s/%s.jpg',tparms.outdir,tparms.outstem);
  if ~exist(fname_out,'file') || parms.forceflag
    args = mmil_parms2args(tparms);
    mmil_subplots(fname,args{:});
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews_mask(fname,parms,src_entities,src_suffix)
  if ~exist('src_entities','var') || isempty(src_entities)
    [src_entities,src_suffix] = split_bids_fstem(get_fstem(fname));
  end
  switch lower(src_suffix)
    case {'t1w','t2w'}
      space_o = src_suffix;
    otherwise
      space_o = '';
  end
  fname_mask_mgz = sprintf('%s/%s_brainmask%s',parms.tmpdir,get_fstem(fname),parms.in_ext);
  if ~exist(fname_mask_mgz,'file')
    fname_mask_mgz = create_brainmask(fname,parms,src_entities,src_suffix);
  end
  if isempty(fname_mask_mgz) || ~exist(fname_mask_mgz,'file'), return; end
  plot_threeviews(fname_mask_mgz,parms,[],src_entities,'mask',space_o);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multislice_measures(fname,parms,measlist,infix)
  fstem = get_fstem(fname);
  [src_entities,~] = split_bids_fstem(fstem);
  if ~exist('infix','var'), infix = []; end
  if ~isempty(infix), infix = ['_' infix]; end
  for i=1:length(measlist)
    meas = measlist{i};
    fname_meas = sprintf('%s/%s%s_%s%s',parms.tmpdir,fstem,infix,meas,parms.in_ext);
    if ~exist(fname_meas,'file'), continue; end
    plot_multislice(fname_meas,parms,[],src_entities,meas_to_bids_suffix(meas));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multiframe_measures(fname,stype,parms,measlist,infix)
  fstem = get_fstem(fname);
  [src_entities,~] = split_bids_fstem(fstem);
  if ~exist('infix','var'), infix = []; end
  if ~isempty(infix), infix = ['_' infix]; end
  for i=1:length(measlist)
    meas = measlist{i};
    fname_meas = sprintf('%s/%s%s_%s%s',parms.tmpdir,fstem,infix,meas,parms.in_ext);
    if ~exist(fname_meas,'file'), continue; end
    meas_bids = meas_to_bids_suffix(meas);
    plot_multiframe(fname_meas,parms,stype,0,src_entities,meas_bids);
    plot_multiframe(fname_meas,parms,stype,1,src_entities,meas_bids);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews_measures(fname,parms,measlist,infix)
  fstem = get_fstem(fname);
  [src_entities,~] = split_bids_fstem(fstem);
  if ~exist('infix','var'), infix = []; end
  if ~isempty(infix), infix = ['_' infix]; end
  for i=1:length(measlist)
    meas = measlist{i};
    fname_meas = sprintf('%s/%s%s_%s%s',parms.tmpdir,fstem,infix,meas,parms.in_ext);
    if ~exist(fname_meas,'file'), continue; end
    plot_threeviews(fname_meas,parms,[],src_entities,meas_to_bids_suffix(meas),'');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bids_suffix = meas_to_bids_suffix(meas)
  % map internal measure names to BIDS-compliant suffix (no underscores)
  map = struct('rms_err','rmse');
  if isfield(map,meas)
    bids_suffix = map.(meas);
  else
    bids_suffix = meas;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_threeviews_rgb(fname,parms,infix)
  fstem = get_fstem(fname);
  [entities,~] = split_bids_fstem(fstem);
  fname_FA = sprintf('%s/%s_%s_FA%s',parms.tmpdir,fstem,infix,parms.in_ext);
  fname_V0 = sprintf('%s/%s_%s_V0%s',parms.tmpdir,fstem,infix,parms.in_ext);
  if ~exist(fname_FA,'file') || ~exist(fname_V0,'file'), return; end
  tparms = [];
  tparms.outdir = parms.outdir;
  tparms.outstem = sprintf('%s_desc-3view_rgbFA',entities);
  tparms.plim = parms.threeviews_plim;
  tparms.fig_size = parms.threeviews_fig_size;
  tparms.tif_dpi = parms.tif_dpi;
  tparms.planestrings = parms.planestrings;
  tparms.slice_fract = parms.threeviews_slice_fract;
  tparms.forceflag = parms.forceflag;
  args = mmil_parms2args(tparms);
  mmil_dti_rgb_subplots(fname_FA,fname_V0,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_motion(fstem,motion_data,parms,filt_flag)
  if ~exist('filt_flag','var') || isempty(filt_flag), filt_flag = 0; end
  [entities,~] = split_bids_fstem(fstem);
  if filt_flag
    fname_out = sprintf('%s/%s_motionfilt.jpg',parms.outdir,entities);
    motion_tseries = motion_data.motion_tseries_filt;
  else
    fname_out = sprintf('%s/%s_motion.jpg',parms.outdir,entities);
    motion_tseries = motion_data.motion_tseries;
  end
  if ~exist(fname_out,'file') || parms.forceflag
    figure;
    set(gcf,'Visible','off');
    subplot(2,1,1);
    plot(motion_tseries(:,1:3));
    axis tight;
    xlabel('time points');
    ylabel('translation (mm)');
    legend({'L-R (dx)','A-P (dy)','I-S (dz)'},'location','EastOutside');
    subplot(2,1,2);
    plot(motion_tseries(:,4:6));
    axis tight;
    xlabel('time points');
    ylabel('rotation (deg)');
    legend({'pitch (rx)','yaw (ry)','roll (rz)'},'location','EastOutside');
    mmil_save_fig(gcf,fname_out,parms.motion_fig_size,parms.tif_dpi,0,0);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_motion_fd(fstem,motion_data,parms,filt_flag)
  if ~exist('filt_flag','var') || isempty(filt_flag), filt_flag = 0; end
  [entities,~] = split_bids_fstem(fstem);
  if filt_flag
    fname_out = sprintf('%s/%s_motionfdfilt.jpg',parms.outdir,entities);
    motion_fd = motion_data.motion_fd_filt;
  else
    fname_out = sprintf('%s/%s_motionfd.jpg',parms.outdir,entities);
    motion_fd = motion_data.motion_fd;
  end
  if ~exist(fname_out,'file') || parms.forceflag
    figure;
    set(gcf,'Visible','off');
    plot(motion_fd);
    axis tight;
    xlabel('time points');
    ylabel('frame-wise displacement (mm)');
    mmil_save_fig(gcf,fname_out,parms.motion_fig_size,parms.tif_dpi,0,0);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = get_fstem(fname)
  [~,fstem] = cmig_get_fileparts(fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scaninfo = get_scaninfo(fname,stype,parms)
  % read scan parameters from BIDS JSON sidecar and bvec/bval files
  scaninfo = [];
  scaninfo.TR   = NaN;
  scaninfo.qmat = [];
  scaninfo.bvals = [];
  scaninfo.bval  = [];
  if parms.verbose
    fprintf('%s: getting scan info...\n',mfilename);
  end
  % derive base path (works for .nii, .nii.gz, .mgh, .mgz)
  fname_base = regexprep(fname,'\.nii.*','');
  fname_base = regexprep(fname_base,'\.(mgh|mgz)$','');
  % load JSON sidecar
  fname_json = [fname_base '.json'];
  json = [];
  if exist(fname_json,'file')
    try
      json = jsondecode(fileread(fname_json));
    catch me
      fprintf('%s: WARNING: failed to read JSON %s: %s\n',...
        mfilename,fname_json,me.message);
    end
  end
  % TR: JSON RepetitionTime is in seconds; store as ms for backward compat
  if ~isempty(json) && isfield(json,'RepetitionTime')
    scaninfo.TR = json.RepetitionTime * 1000;
  end
  % dMRI: load bvec (3×N) and bval (1×N) files
  fname_bvec = [fname_base '.bvec'];
  fname_bval = [fname_base '.bval'];
  if exist(fname_bvec,'file') && exist(fname_bval,'file')
    try
      bvec = load(fname_bvec);  % 3×N
      bval_vec = load(fname_bval);  % 1×N or N×1
      scaninfo.qmat  = bvec';       % N×3
      scaninfo.bvals = bval_vec(:); % N×1 per-volume b-values
      nonzero = scaninfo.bvals(scaninfo.bvals > 0);
      if ~isempty(nonzero)
        scaninfo.bval = max(nonzero); % scalar max, kept for legacy callers
      else
        scaninfo.bval = 0;
      end
    catch me
      fprintf('%s: WARNING: failed to load bvec/bval for %s: %s\n',...
        mfilename,fname,me.message);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seriesinfo = get_seriesinfo(fname,parms)
  % read series metadata from BIDS JSON sidecar
  seriesinfo = [];
  seriesinfo.StudyDate            = '';
  seriesinfo.StudyInstanceUID     = '';
  seriesinfo.SeriesInstanceUID    = '';
  seriesinfo.Manufacturer         = '';
  seriesinfo.ManufacturersModelName = '';
  seriesinfo.DeviceSerialNumber   = '';
  fname_json = regexprep(fname,'\.nii.*','.json');
  if ~exist(fname_json,'file'), return; end
  try
    json = jsondecode(fileread(fname_json));
  catch me
    fprintf('%s: WARNING: failed to read JSON %s: %s\n',...
      mfilename,fname_json,me.message);
    return;
  end
  seriesinfo.Manufacturer           = mmil_getfield(json,'Manufacturer','');
  seriesinfo.ManufacturersModelName = mmil_getfield(json,'ManufacturersModelName','');
  seriesinfo.DeviceSerialNumber     = mmil_getfield(json,'DeviceSerialNumber','');
  seriesinfo.StudyDate              = mmil_getfield(json,'AcquisitionDateTime','');
  seriesinfo.StudyInstanceUID       = mmil_getfield(json,'StudyInstanceUID','');
  seriesinfo.SeriesInstanceUID      = mmil_getfield(json,'SeriesInstanceUID','');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imginfo = get_imginfo(fname,parms)
  imginfo = [];
  if parms.verbose
    fprintf('%s: getting image info...\n',mfilename);
  end
  [imginfo.M,imginfo.volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.tmpdir);
  imginfo.orient  = fs_read_orient([],imginfo.M);
  imginfo.voxsize = fs_voxsize(imginfo.M);
  imginfo.voxvol  = prod(imginfo.voxsize);
  imginfo.nvox    = prod(imginfo.volsz(1:3));
  imginfo.nreps   = imginfo.volsz(4);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_vol = calc_brain_vol(fname_brainmask,auto_qcinfo,parms,fstem_override)
  brain_vol = [];
  if ~exist('fstem_override','var') || isempty(fstem_override)
    fstem = get_fstem(fname_brainmask);
    fstem = regexprep(fstem,'_space-[^_]+_desc-brain_mask$','');
    fstem = regexprep(fstem,'_desc-brain_mask$','');
    fstem = regexprep(fstem,'_brainmask$','');
  else
    fstem = fstem_override;
  end
  fname_out = sprintf('%s/%s_vol.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating brain volume...\n',mfilename);
    end
    [vol,M] = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol>0);
    nvox_brain = length(ind_brain);
    brain_vol = nvox_brain * auto_qcinfo.imginfo.voxvol;
    save(fname_out,'brain_vol','nvox_brain');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_SNR = calc_brain_SNR(fname,fname_brainmask,parms)
  brain_SNR = [];
  fstem = get_fstem(fname);
  fname_out = sprintf('%s/%s_brain_SNR.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating SNR within brain...\n',mfilename);
    end
    [vol,M] = fs_load_mgh(fname,[],1);
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);
    vals = vol(ind_brain);
    brain_SNR.mean   = mean(vals);
    brain_SNR.std    = std(vals);
    brain_SNR.SNR    = brain_SNR.mean / brain_SNR.std;
    brain_SNR.min    = min(vals);
    brain_SNR.max    = max(vals);
    brain_SNR.median = median(vals);
    brain_SNR.nvox_max = nnz(vals>=brain_SNR.max);
    brain_SNR.fvox_max = brain_SNR.nvox_max / length(vals);
    save(fname_out,'brain_SNR');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nonbrain_SNR = calc_nonbrain_SNR(fname,fname_brainmask,parms)
  nonbrain_SNR = [];
  fstem = get_fstem(fname);
  fname_out = sprintf('%s/%s_nonbrain_SNR.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating SNR outside brain...\n',mfilename);
    end
    [vol,M,tmp,volsz] = fs_load_mgh(fname);
    voxsize = fs_voxsize(M);
    vol_brain = fs_load_mgh(fname_brainmask);
    vol_image = 1.0*(vol>0);
    bb_v = ceil(voxsize * parms.bb_buffer);
    ind_brain = find(vol_brain>0);
    [x,y,z] = ind2sub(volsz,ind_brain);
    bb_x = minmax(x');
    bb_y = minmax(y');
    bb_z = minmax(z');
    vol_box = ones(volsz);
    vol_box(1:max(1,bb_x(1)-bb_v(1)),:,:) = 0;
    vol_box(min(volsz(1),bb_x(2)+bb_v(1)):volsz(1),:,:) = 0;
    vol_box(:,1:max(1,bb_y(1)-bb_v(2)),:) = 0;
    vol_box(:,min(volsz(2),bb_y(2)+bb_v(2)):volsz(2),:) = 0;
    vol_box(:,:,1:max(1,bb_z(1)-bb_v(3))) = 0;
    vol_box(:,:,min(volsz(3),bb_z(2)+bb_v(3)):volsz(3)) = 0;
    vol_nonbrain = 1.0*(vol_box==0) .* vol_image;
    ind_nonbrain = find(vol_nonbrain);
    vals = vol(ind_nonbrain);
    nonbrain_SNR.mean = mean(vals);
    nonbrain_SNR.std  = std(vals);
    nonbrain_SNR.SNR  = nonbrain_SNR.mean / nonbrain_SNR.std;
    save(fname_out,'nonbrain_SNR');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_brainmask,auto_qcinfo] = calc_fmri_multiframe_metrics(fname,auto_qcinfo,parms)
  fname_brainmask = [];
  fstem = get_fstem(fname);
  [src_entities,src_suffix] = split_bids_fstem(fstem);
  auto_qcinfo.brain_tSNR = [];
  auto_qcinfo.motion = [];
  auto_qcinfo.fwhm = [];
  auto_qcinfo.fwhm_raw = [];
  if auto_qcinfo.imginfo.nreps<parms.fMRI_minTRs
    fname_brainmask = create_brainmask(fname,parms,src_entities,src_suffix);
    return;
  end
  fname_out  = sprintf('%s/%s_fmri_multiframe_metrics.mat',parms.tmpdir,fstem);
  fname_mean = sprintf('%s/%s_desc-mean_%s.nii.gz',parms.tmpdir,src_entities,src_suffix);
  % validate cached fname_mean before trusting it as a cache hit; deletes
  % it if corrupt/truncated so it self-heals by regenerating below
  check_for_brainmask(fname_mean,parms);
  if ~exist(fname_out,'file') || ~exist(fname_mean,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating fMRI multiframe metrics...\n',mfilename);
    end
    [fname_snr,fname_mean_tmp,fname_std_tmp,fname_motion_tmp] = ...
      calc_fmri_tSNR(fname,auto_qcinfo,parms);
    fs_mri_convert(fname_mean_tmp,fname_mean);
    % mask from the exported mean, so both cache paths share one mask name
    fname_brainmask = create_brainmask(fname_mean,parms,src_entities,src_suffix);
    if isempty(fname_brainmask), return; end
    brain_tSNR = calc_brain_tSNR(fname_snr,fname_brainmask,parms);
    motion     = calc_fmri_motion(fname,auto_qcinfo.scaninfo,parms);
    fwhm       = calc_fwhm(fname,fname_brainmask,parms);
    % second arm: the imported volume, with no motion correction and no
    % detrending, so func also carries a structure-inclusive value comparable
    % to fmap and dwi. Nearly orthogonal to motion, where the detrended arm is
    % largely motion variance
    fwhm_raw   = calc_fwhm_raw(fname,fname_brainmask,parms);
    save(fname_out,'brain_tSNR','motion','fwhm','fwhm_raw');
  else
    fname_brainmask = create_brainmask(fname_mean,parms,src_entities,src_suffix);
    if isempty(fname_brainmask), return; end
    load(fname_out);
  end
  auto_qcinfo.brain_tSNR = brain_tSNR;
  auto_qcinfo.motion = motion;
  auto_qcinfo.fwhm = fwhm;
  if exist('fwhm_raw','var'), auto_qcinfo.fwhm_raw = fwhm_raw; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_for_brainmask(fname,parms)
  % check that a cached fname (e.g. fname_mean) is still usable before it
  % is trusted as a cache hit; deletes it if corrupt/truncated so the
  % caller's exist() check sees it missing and regenerates it
  if ~exist(fname,'file'), return; end
  fname_mask = sprintf('%s/%s_brainmask%s',parms.tmpdir,get_fstem(fname),parms.in_ext);
  if ~exist(fname_mask,'file')
    try
      fs_load_nifti(fname,1);
    catch me
      fprintf('%s: WARNING: %s\n',mfilename,me.message);
      if exist(fname,'file'), delete(fname); end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_snr,fname_mean,fname_std,fname_motion] = calc_fmri_tSNR(fname,auto_qcinfo,parms)
  fname_snr = []; fname_mean = []; fname_std = []; fname_motion = [];
  tparms = parms;
  tparms.TR = auto_qcinfo.scaninfo.TR/1000;
  if isnan(tparms.TR), tparms.TR = 1; end
  tparms.outdir = parms.tmpdir;
  nreps = max(auto_qcinfo.imginfo.nreps,1);
  tparms.skipTRs = min(parms.skipTRs,nreps-1);
  % the tSNR volume takes the tsnr suffix of BOLD_meas_list
  tparms.snr_suffix = '_tsnr';
  args = mmil_parms2args(tparms,parms.tSNR_tags);
  [fname_snr,fname_mean,fname_std,fname_motion] = mmil_calc_tSNR(fname,args{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_data = calc_fmri_motion(fname,scaninfo,parms)
  motion_data = [];
  fstem = get_fstem(fname);
  tstem = sprintf('%s/%s',parms.tmpdir,fstem);
  if parms.skipTRs, tstem = sprintf('%s_skip%dTRs',tstem,parms.skipTRs); end
  tstem = sprintf('%s_mcf',tstem);
  if strcmp(parms.method,'afni')
    fname_motion_txt = sprintf('%s_motion.1D',tstem);
  else
    fname_motion_txt = sprintf('%s.par',tstem);
  end
  fname_motion = sprintf('%s/%s_motion.mat',parms.tmpdir,fstem);
  if ~exist(fname_motion,'file') || parms.forceflag
    if ~exist(fname_motion_txt,'file'), return; end
    if strcmp(parms.method,'afni')
      motion_tseries = mmil_load_motion_1D(fname_motion_txt);
    else
      motion_tseries = mmil_load_motion_par(fname_motion_txt);
    end
    [motion_stats,motion_fd] = mmil_motion_stats(motion_tseries,parms.motion_radius);
    [motion_stats_nody,motion_fd_nody] = mmil_motion_stats(motion_tseries,...
      parms.motion_radius,[],parms.thresholds,1);
    if isnan(scaninfo.TR)
      motion_tseries_filt = motion_tseries;
      fprintf('%s: WARNING: skipping motion filtering because TR is NaN\n',mfilename);
    else
      motion_tseries_filt = abcd_resp_filter_motion(motion_tseries,...
        scaninfo.TR/1000,parms.resp_low,parms.resp_high);
    end
    [motion_stats_filt,motion_fd_filt] = mmil_motion_stats(...
      motion_tseries_filt,parms.motion_radius);
    [motion_stats_filt_nody,motion_fd_filt_nody] = mmil_motion_stats(...
      motion_tseries_filt,parms.motion_radius,[],parms.thresholds,1);
    save(fname_motion,'motion_tseries',...
                      'motion_fd','motion_stats',...
                      'motion_fd_nody','motion_stats_nody',...
                      'motion_tseries_filt',...
                      'motion_fd_filt','motion_stats_filt',...
                      'motion_fd_filt_nody','motion_stats_filt_nody');
  end
  motion_data = load(fname_motion);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fwhm = calc_fwhm(fname,fname_brainmask,parms)
  % the SynthStrip mask built from the exported mean image, so func matches
  % dwi and fmap; AFNI's -automask would otherwise be used here alone
  fwhm = [];
  fstem = get_fstem(fname);
  if parms.skipTRs
    fname_mcf = sprintf('%s/%s_skip%dTRs_mcf.nii',parms.tmpdir,fstem,parms.skipTRs);
  else
    fname_mcf = sprintf('%s/%s_mcf.nii',parms.tmpdir,fstem);
  end
  fname_out = sprintf('%s/%s_fwhm.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if ~exist(fname_mcf,'file'), return; end
    tparms = [];
    tparms.fname_mask = fname_brainmask;
    tparms.outdir = parms.tmpdir;
    tparms.verbose = parms.verbose;
    tparms.forceflag = parms.forceflag;
    args = mmil_parms2args(tparms);
    fwhm = cmig_3dFWHM(fname_mcf,args{:});
    save(fname_out,'fwhm');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fwhm = calc_fwhm_raw(fname,fname_brainmask,parms)
  % smoothness of the input volume itself, with no motion correction and no
  % detrending, for series that have no usable time series to detrend:
  % AFNI removes nothing from a single frame, and multi-frame field maps are
  % repeat measurements, so detrending only some of them would split one
  % column across two regimes. These values include anatomical structure and
  % are not comparable to the detrended fMRI values
  fwhm = [];
  fstem = get_fstem(fname);
  fname_out = sprintf('%s/%s_fwhm_raw.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if ~exist(fname,'file'), return; end
    tparms = [];
    tparms.fname_mask = fname_brainmask;
    tparms.outdir = parms.tmpdir;
    tparms.detrend_flag = 0;
    tparms.verbose = parms.verbose;
    tparms.forceflag = parms.forceflag;
    args = mmil_parms2args(tparms);
    fwhm = cmig_3dFWHM(fname,args{:});
    save(fname_out,'fwhm');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_tSNR = calc_brain_tSNR(fname_snr,fname_brainmask,parms)
  brain_tSNR = [];
  fstem = get_fstem(fname_snr);
  fname_out = sprintf('%s/%s_brain_tSNR.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating tSNR within brain...\n',mfilename);
    end
    vol_snr   = fs_load_mgh(fname_snr);
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);
    vals = vol_snr(ind_brain);
    brain_tSNR.mean   = mean(vals);
    brain_tSNR.median = median(vals);
    brain_tSNR.std    = std(vals);
    save(fname_out,'brain_tSNR');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = check_dMRI_fit(fname,fname_brainmask,auto_qcinfo,parms)
  % the compute is a subfunction so that whatever it managed to produce is
  %   assigned even when it gives up part way. Its early returns used to leave
  %   auto_qcinfo without a motion field, discarding a good motion estimate
  %   because a later step failed, which also cost the QC row its motion
  %   columns and export_dMRI its plots
  [motion,censor,shells,dmri_diff] = ...
    compute_dMRI_fit(fname,fname_brainmask,auto_qcinfo,parms);
  auto_qcinfo.motion = motion;
  auto_qcinfo.censor = censor;
  auto_qcinfo.shells = shells;
  auto_qcinfo.dmri_diff = dmri_diff;
  auto_qcinfo.brain_DTmeas = [];
  auto_qcinfo.brain_DTerr = [];
  if ~isempty(shells)
    auto_qcinfo.brain_DTmeas = shells.brain_DTmeas;
    % the whole-brain per-frame residual ratio, from the same struct that
    %   carries the shell-wise medians. mmil_getfield rather than a direct read:
    %   a dMRI_fit.mat cached before 08/21/2026 has a brain_DTerr without these
    %   subfields, and the consumer treats an empty value as n/a
    tmp = mmil_getfield(shells,'brain_DTerr',[]);
    if isfield(tmp,'brain_rms_err_ratio_mean')
      auto_qcinfo.brain_DTerr = tmp;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [motion,censor,shells,dmri_diff] = ...
                        compute_dMRI_fit(fname,fname_brainmask,auto_qcinfo,parms)
  motion = []; censor = []; shells = []; dmri_diff = [];
  if auto_qcinfo.imginfo.nreps >= parms.dMRI_minTRs
    [fpath,fstem] = cmig_get_fileparts(fname);
    fname_out = sprintf('%s/%s_dMRI_fit.mat',parms.tmpdir,fstem);
    if ~exist(fname_out,'file') || parms.forceflag
      if parms.verbose
        fprintf('%s: checking dMRI fit for %s...\n',mfilename,fname);
      end
      % perform eddy current and motion correction
      %   and calculate mean frame to frame motion
      try
        [motion,censor,fname_mc,fname_qmat] = ...
          calc_dmri_motion(fname,auto_qcinfo.scaninfo,parms);
      catch me
        fprintf('%s: WARNING: calc_dmri_motion failed: %s\n',mfilename,me.message);
        return;
      end
      % load raw data, for censoring detection
      if parms.verbose
        fprintf('%s: loading %s...\n',mfilename,fname);
      end
      data = [];
      [data.vol_raw,M] = fs_load_mgh(fname);
      data.qmat_raw = auto_qcinfo.scaninfo.qmat;
      % load corrected data, for fits and residual metrics
      if parms.verbose
        fprintf('%s: loading %s...\n',mfilename,fname_mc);
      end
      data.vol = fs_load_mgh(fname_mc);
      data.M = M;
      tmp = load(fname_qmat);
      data.qmat = tmp.qmat;
      data.bvals = tmp.bvals;
      if parms.verbose
        fprintf('%s: loading %s...\n',mfilename,fname_brainmask);
      end
      data.vol_mask = fs_load_mgh(fname_brainmask);
      ndirs   = size(data.qmat,1);
      nframes = size(data.vol,4);
      if nframes ~= ndirs
        fprintf('%s: WARNING: mismatch between number of frames (%d) and directions in qmat (%d)\n',...
          mfilename,nframes,ndirs);
        return;
      end
      % RMS difference between raw and corrected data, relative to raw:
      %   summed severity of eddy distortion, motion, and censoring
      dmri_diff = calc_dmri_diff(data);
      if parms.verbose
        fprintf('%s: running dMRI fit calculations with %d frames...\n',...
          mfilename,nframes);
        tic
      end
      p = parms;
      p.data = data;
      % BIDS names, as in postqc dwi, with the entities of the run
      p.outstem = split_bids_fstem(fstem);
      p.bids_flag = 1;
      p.censor_thresh  = parms.dMRI_censor_thresh;
      p.censor_prctile = parms.dMRI_censor_prctile;
      p.cleanupflag = 0;
      p.meas_flag = 1;
      p.fname_DTmeas = sprintf('%s/%s_%s_DTmeas.mat',parms.tmpdir,fstem,parms.dti_corr_suffix);
      p.fname_DTerr  = sprintf('%s/%s_%s_DTerr.mat',parms.tmpdir,fstem,parms.dti_corr_suffix);
      a = mmil_parms2args(p,parms.dMRI_shells_tags);
      shells = autoqc_dMRI_shells([],a{:});
      if parms.verbose, toc; end
      if isempty(shells)
        fprintf('%s: WARNING: dMRI fit calculations failed for %s\n',mfilename,fname);
        return;
      end
      % export DTI measure and error volumes for image plots
      if exist(p.fname_DTmeas,'file')
        convert_DTmeas(p.fname_DTmeas,parms);
      end
      if exist(p.fname_DTerr,'file')
        convert_DTerr(p.fname_DTerr,parms);
      end
      save(fname_out,'motion','censor','shells','dmri_diff');
    else
      load(fname_out);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dmri_diff = calc_dmri_diff(data)
  % within-brain RMS difference between raw and corrected data per frame,
  %   relative to the within-brain RMS of the raw data
  dmri_diff = [];
  ind_brain = find(data.vol_mask>0);
  nframes = size(data.vol,4);
  ratio = zeros(1,nframes);
  for f=1:nframes
    vol_tmp = data.vol_raw(:,:,:,f);
    rms_raw = sqrt(mean(vol_tmp(ind_brain).^2));
    vol_tmp = vol_tmp - data.vol(:,:,:,f);
    rms_diff = sqrt(mean(vol_tmp(ind_brain).^2));
    ratio(f) = rms_diff / (rms_raw + eps);
  end
  dmri_diff.rms_ratio = ratio;
  dmri_diff.rms_ratio_mean = mean(ratio);
  dmri_diff.rms_ratio_median = median(ratio);
  dmri_diff.rms_ratio_std = std(ratio);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [motion_data,censor_info,fname_mc,fname_qmat] = ...
                                      calc_dmri_motion(fname,scaninfo,parms)
  fname_mc = []; fname_qmat = []; motion_data = []; censor_info = [];
  [fpath,fstem,fext] = fileparts(fname);
  fname_mc = sprintf('%s/%s_%s%s',parms.tmpdir,fstem,parms.dti_corr_suffix,fext);
  fname_qmat   = sprintf('%s/%s_%s_qmat.mat',parms.tmpdir,fstem,parms.dti_corr_suffix);
  fname_censor = sprintf('%s/%s_censor.mat',parms.tmpdir,fstem);
  qmat = scaninfo.qmat;
  % use per-volume bvals from BIDS if available
  if isfield(scaninfo,'bvals') && ~isempty(scaninfo.bvals)
    bvals = scaninfo.bvals;
  else
    bvals = scaninfo.bval*sum(qmat.^2,2);
  end
  parms.B0unwarp_flag = 0;
  parms.ecc_flag = 1;
  parms.censor_flag = 1;
  parms.censor_min_ndirs = parms.censor_min_fract*length(bvals);
  parms.motion_B0uw_flag = 0;
  parms.mc_flag = 1;
  parms.gradunwarp_flag = 0;
  parms.interpm = 1;
  parms.qmat = qmat;
  parms.bvals = bvals;
  parms.fname_out  = fname_mc;
  parms.fname_qmat = fname_qmat;
  if parms.verbose
    if ~exist(fname_qmat,'file') || parms.forceflag
      fprintf('%s: performing eddy current and motion correction...\n',mfilename);
    else
      fprintf('%s: loading motion estimates...\n',mfilename);
    end
    tic;
  end
  % name these dti_correct_data options only when enabled: mmil_args2parms is
  %   strict, so naming them unconditionally would break against an MMPS
  %   dti_correct_data that predates them. Omitting them when off is
  %   behaviorally identical, since they default to off there too. Both default
  %   ON here now, so that fallback only applies if a caller turns them off
  dti_corr_tags = parms.dti_corr_tags;
  if parms.afni_flag
    dti_corr_tags = cat(2,dti_corr_tags,{'afni_flag'});
  end
  if parms.censor_reg_flag
    dti_corr_tags = cat(2,dti_corr_tags,{'censor_reg_flag'});
  end
  args = mmil_parms2args(parms,dti_corr_tags);
  [fname_mc,fname_B0dx,fname_qmat] = dti_correct_data(fname,args{:});
  if parms.verbose, toc; end
  fname_motion = sprintf('%s/%s_motion.mat',parms.tmpdir,fstem);
  if ~exist(fname_motion,'file') || parms.forceflag
    if exist(fname_qmat,'file')
      tmp = load(fname_qmat);
      motion_tseries = tmp.motion_tseries;
      [motion_stats,motion_fd] = mmil_motion_stats(motion_tseries,parms.motion_radius);
      [motion_stats_nody,motion_fd_nody] = mmil_motion_stats(motion_tseries,...
        parms.motion_radius,[],parms.thresholds,1);
      save(fname_motion,'motion_tseries','motion_fd','motion_stats',...
                        'motion_fd_nody','motion_stats_nody');
    else
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_qmat);
    end
  end
  if exist(fname_motion,'file')
    motion_data = load(fname_motion);
  end
  if exist(fname_censor,'file')
    censor_info = load(fname_censor);
    censor_info.nbad_frame_slices = sum(sum(censor_info.censor_mat));
    censor_info.nbad_frames = sum(sum(censor_info.censor_mat,1) > 0);
    censor_info.nbad_slices = sum(sum(censor_info.censor_mat,2) > 0);
    censor_info.nbad_frames_per_slice = sum(censor_info.censor_mat,2);
    censor_info.max_nbad_frames_per_slice = max(censor_info.nbad_frames_per_slice);
    censor_info.nbad_slices_per_frame = sum(censor_info.censor_mat,1);
    censor_info.max_nbad_slices_per_frame = max(censor_info.nbad_slices_per_frame);
  else
    fprintf('%s: WARNING: file %s not found\n',mfilename,fname_censor);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_DTfit = calc_DTfit(fname,fname_qmat,parms)
  fstem = get_fstem(fname);
  fname_DTfit = sprintf('%s/%s_DTfit.mat',parms.tmpdir,fstem);
  if ~exist(fname_DTfit,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: loading %s...\n',mfilename,fname);
    end
    [vol,M] = fs_load_mgh(fname);
    load(fname_qmat);
    nframes = size(qmat,1);
    if parms.verbose
      fprintf('%s: running tensor calculations with %d frames...\n',mfilename,nframes);
      tic
    end
    DTfit = dti_fit_tensor(vol,qmat,...
      'bvals',bvals,'M',M,'nonlin_flag',parms.DTfit_nonlin_flag);
    if parms.verbose, toc; end
    if isempty(DTfit)
      fprintf('%s: WARNING: tensor calculations failed for %s\n',mfilename,fname);
      fname_DTfit = [];
      return;
    end
    save(fname_DTfit,'DTfit','M');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_DTmeas = calc_DTmeas(fname_DTfit,parms)
  fname_DTmeas = regexprep(fname_DTfit,'DTfit','DTmeas');
  if ~exist(fname_DTmeas,'file') || parms.forceflag
    if parms.verbose, fprintf('%s: calculating DT measures...\n',mfilename); tic; end
    load(fname_DTfit);
    DTmeas = dti_calc_DTmeas(DTfit,0);
    save(fname_DTmeas,'DTmeas','M');
    if parms.verbose, toc; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_DTmeas(fname_DTmeas,parms)
  fstem = get_fstem(fname_DTmeas);
  fstem = regexprep(fstem,'_DTmeas','');
  fstem = sprintf('%s/%s',parms.tmpdir,fstem);
  all_exist = 1;
  % NOTE: .nii only when nii_flag asks for it. Requiring it unconditionally
  %   made this guard unsatisfiable, since dti_convert_DTmeas empties
  %   nii_measlist when nii_flag is false and so never writes one -- so
  %   all_exist could never reach 1 and the branch below reloaded the .mat and
  %   reconverted on every run, forceflag or not
  ext_list = {parms.in_ext};
  if parms.nii_flag, ext_list{end+1} = '.nii'; end
  for m=1:length(parms.DTI_conv_meas_list)
    meas = parms.DTI_conv_meas_list{m};
    for k=1:length(ext_list)
      fname = sprintf('%s_%s%s',fstem,meas,ext_list{k});
      if ~exist(fname,'file'), all_exist = 0; break; end
    end
    if ~all_exist, break; end
  end
  if ~all_exist || parms.forceflag
    if parms.verbose, fprintf('%s: converting DT measures...\n',mfilename); tic; end
    load(fname_DTmeas);
    parms.M = M;
    parms.measlist = parms.DTI_conv_meas_list;
    parms.nii_measlist = parms.measlist;
    args = mmil_parms2args(parms,parms.convert_DTmeas_tags);
    dti_convert_DTmeas(DTmeas,fstem,args{:});
    if parms.verbose, toc; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function brain_DTmeas = calc_brain_DTmeas(fname_DTmeas,fname_brainmask,parms)
  brain_DTmeas = [];
  fname_brain_DTmeas = regexprep(fname_DTmeas,'DTmeas','brain_DTmeas');
  if ~exist(fname_brain_DTmeas,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: calculating DTI measures within brain...\n',mfilename);
    end
    load(fname_DTmeas);
    vol_brain = fs_load_mgh(fname_brainmask);
    ind_brain = find(vol_brain>0);
    for m=1:length(parms.DTI_meas_list)
      meas = parms.DTI_meas_list{m};
      if strcmp(meas,'MD')
        vol = 1000*mean(DTmeas.volE,4);
        vals = vol(ind_brain);
      else
        vals = DTmeas.(['vol' meas])(ind_brain);
      end
      brain_DTmeas.(meas).mean   = mean(vals);
      brain_DTmeas.(meas).median = median(vals);
      brain_DTmeas.(meas).std    = std(vals);
    end
    save(fname_brain_DTmeas,'brain_DTmeas');
  else
    load(fname_brain_DTmeas);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_DTerr(fname_DTerr,parms)
  fstem = get_fstem(fname_DTerr);
  fstem = regexprep(fstem,'_DTerr','');
  fstem = sprintf('%s/%s',parms.tmpdir,fstem);
  DTerr = [];
  for m=1:length(parms.DTI_err_list)
    errname  = parms.DTI_err_list{m};
    fname_mgz = sprintf('%s_%s%s',fstem,errname,parms.in_ext);
    fname_nii = sprintf('%s_%s.nii',fstem,errname);
    if ~exist(fname_mgz,'file') || parms.forceflag
      if isempty(DTerr)
        if parms.verbose, fprintf('%s: converting DT err...\n',mfilename); tic; end
        load(fname_DTerr);
      end
      vol = DTerr.(['vol_' errname]);
      fs_save_mgh(vol,fname_mgz,M);
    end
    % NOTE: the nifti copy is for viewing the error map in AFNI and nothing
    %   reads it -- the plot functions build their names from parms.in_ext --
    %   so it is written only when asked for
    if parms.nii_flag
      args = mmil_parms2args(parms,parms.mri_convert_tags);
      fs_mri_convert(fname_mgz,fname_nii,args{:});
    end
    if parms.verbose && ~isempty(DTerr), toc; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [line_artifact,line2_artifact] = detect_lines(fname,parms,fname_brainmask)
  % original line_* detector plus the redesigned fft detector (line2_*), computed
  % in a single mgz load and cached separately so each can be recomputed without
  % disturbing the other. line_* stays byte-identical to the pre-fft behavior.
  line_artifact = [];
  line2_artifact = [];
  if ~exist('fname_brainmask','var'), fname_brainmask = []; end
  [~,fstem] = cmig_get_fileparts(fname);
  fname_out  = sprintf('%s/%s_line_artifact.mat',parms.tmpdir,fstem);
  fname_out2 = sprintf('%s/%s_line2_artifact.mat',parms.tmpdir,fstem);
  need_line  = ~exist(fname_out,'file')  || parms.forceflag;
  need_line2 = ~exist(fname_out2,'file') || parms.forceflag;
  if need_line || need_line2
    if parms.verbose
      fprintf('%s: detecting line artifacts for %s...\n',mfilename,fname);
    end
    % nargout>2 -> mmil_detect_lines also runs the fft detector (line2_*); pass the
    % brain mask so SynthStrip is not recomputed
    [scores,counts,scores_fft,nslices_fft,extents_fft] = ...
      mmil_detect_lines(fname,'outdir',parms.tmpdir,'verbose',parms.verbose,...
        'fname_mask',fname_brainmask);
    entities = regexprep(fstem,'_[^_]+$','');
    d_lines = dir(sprintf('%s/%s_f*_lines_slice*.jpg',parms.tmpdir,fstem));
    names = {d_lines.name};
    % original detector plot (exclude the fft plots)
    idx = find(cellfun(@isempty,regexp(names,'_fft_lines_slice','once')),1);
    if ~isempty(idx)
      dest = sprintf('%s/%s_lines.jpg',parms.outdir,entities);
      copyfile(fullfile(parms.tmpdir,d_lines(idx).name),dest);
    end
    % fft detector plot
    idx2 = find(~cellfun(@isempty,regexp(names,'_fft_lines_slice','once')),1);
    if ~isempty(idx2)
      dest2 = sprintf('%s/%s_lines2.jpg',parms.outdir,entities);
      copyfile(fullfile(parms.tmpdir,d_lines(idx2).name),dest2);
    end
  end
  % line_* (original detector)
  if need_line
    line_artifact.scores = scores;
    line_artifact.counts = counts;
    line_artifact.nframes    = nnz(line_artifact.scores);
    line_artifact.max_score  = max(line_artifact.scores);
    line_artifact.max_count  = max(line_artifact.counts);
    if line_artifact.max_score > 0
      line_artifact.mean_score = mean(line_artifact.scores(line_artifact.scores>0));
      line_artifact.mean_count = mean(line_artifact.counts(line_artifact.scores>0));
    else
      line_artifact.mean_score = 0;
      line_artifact.mean_count = 0;
    end
    if parms.verbose
      fprintf('%s: frames with non-zero line score: %d\n',mfilename,line_artifact.nframes);
    end
    save(fname_out,'line_artifact');
  else
    load(fname_out);
  end
  % line2_* (fft detector) - NaN-aware
  if need_line2
    line2_artifact = summarize_line2(scores_fft,nslices_fft,extents_fft);
    if parms.verbose
      fprintf('%s: frames with fft line artifact: %d\n',mfilename,line2_artifact.nframes);
    end
    save(fname_out2,'line2_artifact');
  else
    S = load(fname_out2);
    % re-reduce from the cached per-frame arrays, so a cache written before a
    % summary field existed still yields that field
    line2_artifact = summarize_line2(S.line2_artifact.scores,...
      S.line2_artifact.nslices,S.line2_artifact.extents);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function line2_artifact = summarize_line2(scores,nslices,extents)
  % NaN-aware reduction of per-frame fft outputs. NaN = a frame the detector could
  % not evaluate (distinct from 0 = evaluated, no line). line2_extent is taken from
  % the max-score frame (extent is only meaningful where a line is present), and
  % line2_best_* from the lowest-scoring frame (the least banded).
  % A frame is "banded" when nslices > 0; score is continuous (clean frames still
  % score ~4-6), so nslices - not score - defines the banded-frame mask that
  % nframes/mean_* are conditional on.
  line2_artifact = [];
  line2_artifact.scores = scores;
  line2_artifact.nslices = nslices;
  line2_artifact.extents = extents;
  valid = ~isnan(scores);
  pos = valid & (nslices > 0);
  line2_artifact.nframes = nnz(pos);
  if any(valid)
    [line2_artifact.max_score,imax] = max(scores);   % max ignores NaN
    line2_artifact.max_nslices = max(nslices);
    line2_artifact.extent = extents(imax);
    % the frame to correct from, where a multi-frame field map offers a choice;
    % one-based, ties taking the first, min ignoring NaN as max does
    [line2_artifact.best_score,line2_artifact.best_frame] = min(scores);
  else
    line2_artifact.max_score = NaN;
    line2_artifact.max_nslices = NaN;
    line2_artifact.extent = NaN;
    line2_artifact.best_score = NaN;
    line2_artifact.best_frame = NaN;
  end
  if line2_artifact.nframes > 0
    line2_artifact.mean_score = mean(scores(pos));
    line2_artifact.mean_nslices = mean(nslices(pos));
  else
    line2_artifact.mean_score = 0;
    line2_artifact.mean_nslices = 0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FOV_cutoff = detect_FOV_cutoff(fname_mask,parms,src_entities,src_suffix)
  FOV_cutoff = [];
  [~,fstem] = cmig_get_fileparts(fname_mask);
  if ~exist('src_entities','var') || isempty(src_entities)
    src_entities = regexprep(fstem,'_brainmask$','');
    src_entities = regexprep(src_entities,'_desc-brain_mask$','');
    src_entities = regexprep(src_entities,'_space-[^_]+_desc-brain_mask$','');
    src_suffix = '';
  end
  fname_out = sprintf('%s/%s_fovcutoff.mat',parms.tmpdir,fstem);
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: detecting FOV cutoff for %s...\n',mfilename,fname_mask);
    end
    fname_img = sprintf('%s/%s_fovcutoff.jpg',parms.outdir,src_entities);
    [FOV_cutoff.cutoff,FOV_cutoff.dorsal,FOV_cutoff.ventral] = ...
      mmil_detect_cutoff(fname_mask,'fname_plot',fname_img,'forceflag',parms.forceflag);
    save(fname_out,'FOV_cutoff');
  else
    load(fname_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_rawqc_tsv(all_qcstruct,fname_tsv,parms)
  % collect all unique field names across all rows (preserving insertion order)
  all_fields = {};
  for i=1:length(all_qcstruct)
    flds = fieldnames(all_qcstruct{i});
    for k=1:length(flds)
      if ~ismember(flds{k},all_fields)
        all_fields{end+1} = flds{k};
      end
    end
  end
  % build data cell matrix (nrows x ncols)
  nrows = length(all_qcstruct);
  ncols = length(all_fields);
  data = cell(nrows,ncols);
  for i=1:nrows
    row = all_qcstruct{i};
    for k=1:ncols
      fn = all_fields{k};
      if isfield(row,fn)
        val = row.(fn);
        if isnumeric(val) && numel(val)==1 && isnan(val)
          data{i,k} = 'n/a';
        else
          data{i,k} = val;
        end
      else
        data{i,k} = 'n/a';
      end
    end
  end
  mmil_write_csv(fname_tsv,cat(1,all_fields,data),'separator','\t','quote_flag',0);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcstruct = convert_qcinfo(qcinfo,thresholds,participant_id,session_id)
  if ~mmil_check_nargs(nargin,4), return; end
  nthresh = length(thresholds);
  for j=1:length(qcinfo)
    qcstruct(j).participant_id = participant_id;
    qcstruct(j).session_id = session_id;
    % stype-independent metrics
    qcstruct(j).fstem              = qcinfo(j).fstem;
    qcstruct(j).acq                = get_acq_entity(qcinfo(j).fstem);
    qcstruct(j).stype              = qcinfo(j).stype;
    qcstruct(j).StudyDate          = qcinfo(j).seriesinfo.StudyDate;
    qcstruct(j).StudyInstanceUID   = qcinfo(j).seriesinfo.StudyInstanceUID;
    qcstruct(j).SeriesInstanceUID  = qcinfo(j).seriesinfo.SeriesInstanceUID;
    qcstruct(j).Manufacturer       = qcinfo(j).seriesinfo.Manufacturer;
    qcstruct(j).ManufacturersModelName = qcinfo(j).seriesinfo.ManufacturersModelName;
    qcstruct(j).DeviceSerialNumber = qcinfo(j).seriesinfo.DeviceSerialNumber;
    qcstruct(j).nvox   = qcinfo(j).imginfo.volsz(1:3);
    qcstruct(j).nreps  = qcinfo(j).imginfo.nreps;
    qcstruct(j).orient = qcinfo(j).imginfo.orient;
    qcstruct(j).voxvol = qcinfo(j).imginfo.voxvol;
    qcstruct(j).brainvol = mmil_getfield(qcinfo(j),'brain_vol',NaN);
    % stype-specific metrics
    switch qcstruct(j).stype
      case {'T1w','T2w','T2starw','FLAIR','PDw','epi','angio',...
            'VFA','MP2RAGE','FLASH','TB1map','TB1EPI','TB1AFI','TB1TFL','RB1COR',...
            'MPR','FLASHhi','XetaT2','QALAS','B1map'} % include legacy names
        % line artifact (fmap/epi only; guarded so anat types skip it)
        if isfield(qcinfo,'line_artifact') && ~isempty(qcinfo(j).line_artifact)
          qcstruct(j).line_nframes     = qcinfo(j).line_artifact.nframes;
          qcstruct(j).line_max_score   = qcinfo(j).line_artifact.max_score;
          qcstruct(j).line_mean_score  = qcinfo(j).line_artifact.mean_score;
          qcstruct(j).line_max_count   = qcinfo(j).line_artifact.max_count;
          qcstruct(j).line_mean_count  = qcinfo(j).line_artifact.mean_count;
        end
        % line2 artifact (fft detector)
        if isfield(qcinfo,'line2_artifact') && ~isempty(qcinfo(j).line2_artifact)
          qcstruct(j).line2_nframes      = qcinfo(j).line2_artifact.nframes;
          qcstruct(j).line2_max_score    = qcinfo(j).line2_artifact.max_score;
          qcstruct(j).line2_mean_score   = qcinfo(j).line2_artifact.mean_score;
          qcstruct(j).line2_max_nslices  = qcinfo(j).line2_artifact.max_nslices;
          qcstruct(j).line2_mean_nslices = qcinfo(j).line2_artifact.mean_nslices;
          qcstruct(j).line2_extent       = qcinfo(j).line2_artifact.extent;
          % which frame of a multi-frame field map to correct from, and its
          % score; defaulted so a qcinfo cached before the fields existed reads
          qcstruct(j).line2_best_frame = ...
            mmil_getfield(qcinfo(j).line2_artifact,'best_frame',NaN);
          qcstruct(j).line2_best_score = ...
            mmil_getfield(qcinfo(j).line2_artifact,'best_score',NaN);
        end
        % fmap/epi only; check_sMRI computes no cutoff, so anat rows skip it
        if isfield(qcinfo,'FOV_cutoff') && ~isempty(qcinfo(j).FOV_cutoff)
          qcstruct(j).cutoff         = qcinfo(j).FOV_cutoff.cutoff;
          qcstruct(j).dorsal_cutoff  = qcinfo(j).FOV_cutoff.dorsal;
          qcstruct(j).ventral_cutoff = qcinfo(j).FOV_cutoff.ventral;
        end
        if isfield(qcinfo,'brain_SNR') && ~isempty(qcinfo(j).brain_SNR)
          qcstruct(j).brain_mean      = qcinfo(j).brain_SNR.mean;
          qcstruct(j).brain_std       = qcinfo(j).brain_SNR.std;
          qcstruct(j).brain_SNR       = qcinfo(j).brain_SNR.SNR;
          qcstruct(j).brain_min       = qcinfo(j).brain_SNR.min;
          qcstruct(j).brain_max       = qcinfo(j).brain_SNR.max;
          qcstruct(j).brain_median    = qcinfo(j).brain_SNR.median;
          qcstruct(j).brain_nvox_max  = qcinfo(j).brain_SNR.nvox_max;
          qcstruct(j).brain_fvox_max  = qcinfo(j).brain_SNR.fvox_max;
        end
        if isfield(qcinfo,'nonbrain_SNR') && ~isempty(qcinfo(j).nonbrain_SNR)
          qcstruct(j).nonbrain_mean = qcinfo(j).nonbrain_SNR.mean;
          qcstruct(j).nonbrain_std  = qcinfo(j).nonbrain_SNR.std;
          qcstruct(j).nonbrain_SNR  = qcinfo(j).nonbrain_SNR.SNR;
        end
        % spatial smoothness, for fmap/epi and the structural types; guarded
        % so any type that does not compute it is skipped
        if isfield(qcinfo,'fwhm') && ~isempty(qcinfo(j).fwhm)
          qcstruct(j).fwhm_x = qcinfo(j).fwhm(1);
          qcstruct(j).fwhm_y = qcinfo(j).fwhm(2);
          qcstruct(j).fwhm_z = qcinfo(j).fwhm(3);
        end
      case {'bold','BOLD'}
        if isfield(qcinfo,'line_artifact') && ~isempty(qcinfo(j).line_artifact)
          qcstruct(j).line_nframes     = qcinfo(j).line_artifact.nframes;
          qcstruct(j).line_max_score   = qcinfo(j).line_artifact.max_score;
          qcstruct(j).line_mean_score  = qcinfo(j).line_artifact.mean_score;
          qcstruct(j).line_max_count   = qcinfo(j).line_artifact.max_count;
          qcstruct(j).line_mean_count  = qcinfo(j).line_artifact.mean_count;
        end
        % line2 artifact (fft detector)
        if isfield(qcinfo,'line2_artifact') && ~isempty(qcinfo(j).line2_artifact)
          qcstruct(j).line2_nframes      = qcinfo(j).line2_artifact.nframes;
          qcstruct(j).line2_max_score    = qcinfo(j).line2_artifact.max_score;
          qcstruct(j).line2_mean_score   = qcinfo(j).line2_artifact.mean_score;
          qcstruct(j).line2_max_nslices  = qcinfo(j).line2_artifact.max_nslices;
          qcstruct(j).line2_mean_nslices = qcinfo(j).line2_artifact.mean_nslices;
          qcstruct(j).line2_extent       = qcinfo(j).line2_artifact.extent;
        end
        if isfield(qcinfo,'FOV_cutoff') && ~isempty(qcinfo(j).FOV_cutoff)
          qcstruct(j).cutoff         = qcinfo(j).FOV_cutoff.cutoff;
          qcstruct(j).dorsal_cutoff  = qcinfo(j).FOV_cutoff.dorsal;
          qcstruct(j).ventral_cutoff = qcinfo(j).FOV_cutoff.ventral;
        end
        if isfield(qcinfo,'brain_SNR') && ~isempty(qcinfo(j).brain_SNR)
          qcstruct(j).brain_mean      = qcinfo(j).brain_SNR.mean;
          qcstruct(j).brain_std       = qcinfo(j).brain_SNR.std;
          qcstruct(j).brain_SNR       = qcinfo(j).brain_SNR.SNR;
          qcstruct(j).brain_min       = qcinfo(j).brain_SNR.min;
          qcstruct(j).brain_max       = qcinfo(j).brain_SNR.max;
          qcstruct(j).brain_median    = qcinfo(j).brain_SNR.median;
          qcstruct(j).brain_nvox_max  = qcinfo(j).brain_SNR.nvox_max;
          qcstruct(j).brain_fvox_max  = qcinfo(j).brain_SNR.fvox_max;
        end
        qcstruct(j).TR = qcinfo(j).scaninfo.TR/1000;
        if isnan(qcstruct(j).TR), qcstruct(j).TR = 1; end
        if isfield(qcinfo,'motion') && ~isempty(qcinfo(j).motion)
          qcstruct(j).mean_motion      = qcinfo(j).motion.motion_stats.mean_motion;
          qcstruct(j).mean_motion_nody = qcinfo(j).motion.motion_stats_nody.mean_motion;
          qcstruct(j).mean_trans       = qcinfo(j).motion.motion_stats.mean_trans;
          qcstruct(j).mean_trans_nody  = qcinfo(j).motion.motion_stats_nody.mean_trans;
          qcstruct(j).mean_rot         = qcinfo(j).motion.motion_stats.mean_rot;
          qcstruct(j).max_motion       = qcinfo(j).motion.motion_stats.max_motion;
          qcstruct(j).max_trans        = qcinfo(j).motion.motion_stats.max_trans;
          qcstruct(j).max_rot          = qcinfo(j).motion.motion_stats.max_rot;
          qcstruct(j).max_dx           = qcinfo(j).motion.motion_stats.max_dx;
          qcstruct(j).max_dy           = qcinfo(j).motion.motion_stats.max_dy;
          qcstruct(j).max_dz           = qcinfo(j).motion.motion_stats.max_dz;
          qcstruct(j).max_rx           = qcinfo(j).motion.motion_stats.max_rx;
          qcstruct(j).max_ry           = qcinfo(j).motion.motion_stats.max_ry;
          qcstruct(j).max_rz           = qcinfo(j).motion.motion_stats.max_rz;
          for i=1:nthresh
            qcstruct(j).(sprintf('subthresh_0%0.0f',10*thresholds(i))) = ...
              qcinfo(j).motion.motion_stats.subthresh_nvols(i)*qcstruct(j).TR;
            qcstruct(j).(sprintf('subthresh_0%0.0f_nody',10*thresholds(i))) = ...
              qcinfo(j).motion.motion_stats_nody.subthresh_nvols(i)*qcstruct(j).TR;
          end
        end
        if isfield(qcinfo,'motion') && ~isempty(qcinfo(j).motion) && ...
           isfield(qcinfo(j).motion,'motion_stats_filt')
          qcstruct(j).mean_motion_filt      = qcinfo(j).motion.motion_stats_filt.mean_motion;
          qcstruct(j).mean_motion_filt_nody = qcinfo(j).motion.motion_stats_filt_nody.mean_motion;
          qcstruct(j).mean_trans_filt       = qcinfo(j).motion.motion_stats_filt.mean_trans;
          qcstruct(j).mean_trans_filt_nody  = qcinfo(j).motion.motion_stats_filt_nody.mean_trans;
          qcstruct(j).mean_rot_filt         = qcinfo(j).motion.motion_stats_filt.mean_rot;
          qcstruct(j).max_dx_filt           = qcinfo(j).motion.motion_stats_filt.max_dx;
          qcstruct(j).max_dy_filt           = qcinfo(j).motion.motion_stats_filt.max_dy;
          qcstruct(j).max_dz_filt           = qcinfo(j).motion.motion_stats_filt.max_dz;
          qcstruct(j).max_rx_filt           = qcinfo(j).motion.motion_stats_filt.max_rx;
          qcstruct(j).max_ry_filt           = qcinfo(j).motion.motion_stats_filt.max_ry;
          qcstruct(j).max_rz_filt           = qcinfo(j).motion.motion_stats_filt.max_rz;
          for i=1:nthresh
            qcstruct(j).(sprintf('subthresh_0%0.0f_filt',10*thresholds(i))) = ...
              qcinfo(j).motion.motion_stats_filt.subthresh_nvols(i)*qcstruct(j).TR;
            qcstruct(j).(sprintf('subthresh_0%0.0f_filt_nody',10*thresholds(i))) = ...
              qcinfo(j).motion.motion_stats_filt_nody.subthresh_nvols(i)*qcstruct(j).TR;
          end
        end
        if isfield(qcinfo,'fwhm') && ~isempty(qcinfo(j).fwhm)
          qcstruct(j).fwhm_x = qcinfo(j).fwhm(1);
          qcstruct(j).fwhm_y = qcinfo(j).fwhm(2);
          qcstruct(j).fwhm_z = qcinfo(j).fwhm(3);
        end
        if isfield(qcinfo,'fwhm_raw') && ~isempty(qcinfo(j).fwhm_raw)
          qcstruct(j).fwhm_raw_x = qcinfo(j).fwhm_raw(1);
          qcstruct(j).fwhm_raw_y = qcinfo(j).fwhm_raw(2);
          qcstruct(j).fwhm_raw_z = qcinfo(j).fwhm_raw(3);
        end
        if isfield(qcinfo,'brain_tSNR') && ~isempty(qcinfo(j).brain_tSNR)
          qcstruct(j).brain_tSNR_mean   = qcinfo(j).brain_tSNR.mean;
          qcstruct(j).brain_tSNR_median = qcinfo(j).brain_tSNR.median;
          qcstruct(j).brain_tSNR_std    = qcinfo(j).brain_tSNR.std;
        end
        if isfield(qcinfo,'pepolar') && ~isempty(qcinfo(j).pepolar)
          qcstruct(j).pepolar_mi_diff = qcinfo(j).pepolar.mi_diff;
          qcstruct(j).pepolar_mi_for  = qcinfo(j).pepolar.mi_for;
          qcstruct(j).pepolar_mi_rev  = qcinfo(j).pepolar.mi_rev;
        end
      case {'dwi','DTI'}
        if isfield(qcinfo,'line_artifact') && ~isempty(qcinfo(j).line_artifact)
          qcstruct(j).line_nframes     = qcinfo(j).line_artifact.nframes;
          qcstruct(j).line_max_score   = qcinfo(j).line_artifact.max_score;
          qcstruct(j).line_mean_score  = qcinfo(j).line_artifact.mean_score;
          qcstruct(j).line_max_count   = qcinfo(j).line_artifact.max_count;
          qcstruct(j).line_mean_count  = qcinfo(j).line_artifact.mean_count;
        end
        % line2 artifact (fft detector)
        if isfield(qcinfo,'line2_artifact') && ~isempty(qcinfo(j).line2_artifact)
          qcstruct(j).line2_nframes      = qcinfo(j).line2_artifact.nframes;
          qcstruct(j).line2_max_score    = qcinfo(j).line2_artifact.max_score;
          qcstruct(j).line2_mean_score   = qcinfo(j).line2_artifact.mean_score;
          qcstruct(j).line2_max_nslices  = qcinfo(j).line2_artifact.max_nslices;
          qcstruct(j).line2_mean_nslices = qcinfo(j).line2_artifact.mean_nslices;
          qcstruct(j).line2_extent       = qcinfo(j).line2_artifact.extent;
        end
        if isfield(qcinfo,'FOV_cutoff') && ~isempty(qcinfo(j).FOV_cutoff)
          qcstruct(j).cutoff         = qcinfo(j).FOV_cutoff.cutoff;
          qcstruct(j).dorsal_cutoff  = qcinfo(j).FOV_cutoff.dorsal;
          qcstruct(j).ventral_cutoff = qcinfo(j).FOV_cutoff.ventral;
        end
        if isfield(qcinfo,'brain_SNR') && ~isempty(qcinfo(j).brain_SNR)
          qcstruct(j).brain_mean      = qcinfo(j).brain_SNR.mean;
          qcstruct(j).brain_std       = qcinfo(j).brain_SNR.std;
          qcstruct(j).brain_SNR       = qcinfo(j).brain_SNR.SNR;
          qcstruct(j).brain_min       = qcinfo(j).brain_SNR.min;
          qcstruct(j).brain_max       = qcinfo(j).brain_SNR.max;
          qcstruct(j).brain_median    = qcinfo(j).brain_SNR.median;
          qcstruct(j).brain_nvox_max  = qcinfo(j).brain_SNR.nvox_max;
          qcstruct(j).brain_fvox_max  = qcinfo(j).brain_SNR.fvox_max;
        end
        qcstruct(j).TR = qcinfo(j).scaninfo.TR/1000;
        if isfield(qcinfo,'motion') && ~isempty(qcinfo(j).motion)
          qcstruct(j).mean_motion      = qcinfo(j).motion.motion_stats.mean_motion;
          qcstruct(j).mean_motion_nody = qcinfo(j).motion.motion_stats_nody.mean_motion;
          qcstruct(j).mean_trans       = qcinfo(j).motion.motion_stats.mean_trans;
          qcstruct(j).mean_trans_nody  = qcinfo(j).motion.motion_stats_nody.mean_trans;
          qcstruct(j).mean_rot         = qcinfo(j).motion.motion_stats.mean_rot;
          qcstruct(j).max_motion       = qcinfo(j).motion.motion_stats.max_motion;
          qcstruct(j).max_trans        = qcinfo(j).motion.motion_stats.max_trans;
          qcstruct(j).max_rot          = qcinfo(j).motion.motion_stats.max_rot;
          qcstruct(j).max_dx           = qcinfo(j).motion.motion_stats.max_dx;
          qcstruct(j).max_dy           = qcinfo(j).motion.motion_stats.max_dy;
          qcstruct(j).max_dz           = qcinfo(j).motion.motion_stats.max_dz;
          qcstruct(j).max_rx           = qcinfo(j).motion.motion_stats.max_rx;
          qcstruct(j).max_ry           = qcinfo(j).motion.motion_stats.max_ry;
          qcstruct(j).max_rz           = qcinfo(j).motion.motion_stats.max_rz;
          for i=1:nthresh
            qcstruct(j).(sprintf('subthresh_0%0.0f',10*thresholds(i))) = ...
              qcinfo(j).motion.motion_stats.subthresh_nvols(i)*qcstruct(j).TR;
            qcstruct(j).(sprintf('subthresh_0%0.0f_nody',10*thresholds(i))) = ...
              qcinfo(j).motion.motion_stats_nody.subthresh_nvols(i)*qcstruct(j).TR;
          end
        end
        if isfield(qcinfo,'censor') && ~isempty(qcinfo(j).censor)
          qcstruct(j).nbad_frame_slices        = qcinfo(j).censor.nbad_frame_slices;
          qcstruct(j).nbad_frames              = qcinfo(j).censor.nbad_frames;
          qcstruct(j).nbad_slices              = qcinfo(j).censor.nbad_slices;
          qcstruct(j).max_nbad_frames_per_slice = qcinfo(j).censor.max_nbad_frames_per_slice;
          qcstruct(j).max_nbad_slices_per_frame = qcinfo(j).censor.max_nbad_slices_per_frame;
          qcstruct(j).ngood_frames = qcstruct(j).nreps - qcstruct(j).nbad_frames;
        end
        if isfield(qcinfo,'brain_DTmeas') && ~isempty(qcinfo(j).brain_DTmeas)
          meas_list = fieldnames(qcinfo(j).brain_DTmeas);
          for m=1:length(meas_list)
            meas = meas_list{m};
            qcstruct(j).([meas '_mean'])   = qcinfo(j).brain_DTmeas.(meas).mean;
            qcstruct(j).([meas '_median']) = qcinfo(j).brain_DTmeas.(meas).median;
            qcstruct(j).([meas '_std'])    = qcinfo(j).brain_DTmeas.(meas).std;
          end
        end
        if isfield(qcinfo,'brain_DTerr') && ~isempty(qcinfo(j).brain_DTerr)
          qcstruct(j).DTerr_mean   = qcinfo(j).brain_DTerr.brain_rms_err_ratio_mean;
          qcstruct(j).DTerr_median = qcinfo(j).brain_DTerr.brain_rms_err_ratio_median;
          qcstruct(j).DTerr_std    = qcinfo(j).brain_DTerr.brain_rms_err_ratio_std;
        end
        if isfield(qcinfo,'dmri_diff') && ~isempty(qcinfo(j).dmri_diff)
          qcstruct(j).diff_rms_mean   = qcinfo(j).dmri_diff.rms_ratio_mean;
          qcstruct(j).diff_rms_median = qcinfo(j).dmri_diff.rms_ratio_median;
          qcstruct(j).diff_rms_std    = qcinfo(j).dmri_diff.rms_ratio_std;
        end
        if isfield(qcinfo,'shells') && ~isempty(qcinfo(j).shells)
          if isfield(qcinfo(j).shells,'RSI_censor_info') && ...
             ~isempty(qcinfo(j).shells.RSI_censor_info)
            censor_info = qcinfo(j).shells.RSI_censor_info;
          elseif isfield(qcinfo(j).shells,'DTI_censor_info') && ...
                 ~isempty(qcinfo(j).shells.DTI_censor_info)
            censor_info = qcinfo(j).shells.DTI_censor_info;
          else
            censor_info = [];
          end
          if ~isempty(censor_info)
            qcstruct(j).nbad_frame_slices        = censor_info.nbad_frame_slices;
            qcstruct(j).nbad_frames              = censor_info.nbad_frames;
            qcstruct(j).nbad_slices              = censor_info.nbad_slices;
            qcstruct(j).max_nbad_frames_per_slice = censor_info.max_nbad_frames_per_slice;
            qcstruct(j).max_nbad_slices_per_frame = censor_info.max_nbad_slices_per_frame;
            qcstruct(j).censor_thresh = censor_info.censor_thresh;
            qcstruct(j).ngood_frames = qcstruct(j).nreps - qcstruct(j).nbad_frames;
            for b=1:qcinfo(j).shells.nbvals
              bval = qcinfo(j).shells.unique_bvals(b);
              qcstruct(j).(sprintf('nframes_b%d',bval))           = censor_info.nframes_per_bval(b);
              qcstruct(j).(sprintf('nbad_frame_slices_b%d',bval)) = censor_info.nbad_frame_slices_per_bval(b);
              qcstruct(j).(sprintf('nbad_frames_b%d',bval))       = censor_info.nbad_frames_per_bval(b);
              qcstruct(j).(sprintf('nbad_slices_b%d',bval))       = censor_info.nbad_slices_per_bval(b);
              qcstruct(j).(sprintf('ngood_frames_b%d',bval))      = censor_info.nframes_per_bval(b) - ...
                                                                     censor_info.nbad_frames_per_bval(b);
            end
          end
          if isfield(qcinfo(j).shells,'brain_FWHM_resid') && ...
             ~isempty(qcinfo(j).shells.brain_FWHM_resid)
            resid = qcinfo(j).shells.brain_FWHM_resid;
            for b=1:length(resid.nframes_bvals)
              bval = qcinfo(j).shells.unique_bvals(b);
              qcstruct(j).(sprintf('FWHMx_resid_b%d',bval)) = resid.mean_brain_FWHMx_bvals(b);
              qcstruct(j).(sprintf('FWHMy_resid_b%d',bval)) = resid.mean_brain_FWHMy_bvals(b);
              qcstruct(j).(sprintf('FWHMz_resid_b%d',bval)) = resid.mean_brain_FWHMz_bvals(b);
            end
            qcstruct(j).FWHMx_resid = resid.mean_brain_FWHMx;
            qcstruct(j).FWHMy_resid = resid.mean_brain_FWHMy;
            qcstruct(j).FWHMz_resid = resid.mean_brain_FWHMz;
          end
          for b=1:qcinfo(j).shells.nbvals
            bval = qcinfo(j).shells.unique_bvals(b);
            qcstruct(j).(sprintf('FWHMx_b%d',bval)) = qcinfo(j).shells.brain_FWHM.mean_brain_FWHMx_bvals(b);
            qcstruct(j).(sprintf('FWHMy_b%d',bval)) = qcinfo(j).shells.brain_FWHM.mean_brain_FWHMy_bvals(b);
            qcstruct(j).(sprintf('FWHMz_b%d',bval)) = qcinfo(j).shells.brain_FWHM.mean_brain_FWHMz_bvals(b);
            qcstruct(j).(sprintf('tSNR_b%d',bval))  = qcinfo(j).shells.brain_tSNR.median_brain_tSNR_bvals(b);
            qcstruct(j).(sprintf('DTerr_norm_b%d',bval)) = qcinfo(j).shells.brain_DTerr.median_brain_rms_err_norm_bvals(b);
            qcstruct(j).(sprintf('DTerr_rel_b%d',bval))  = qcinfo(j).shells.brain_DTerr.median_brain_rms_err_rel_bvals(b);
            qcstruct(j).(sprintf('RSIerr_norm_b%d',bval)) = qcinfo(j).shells.brain_RSIerr.median_brain_rms_err_norm_bvals(b);
            qcstruct(j).(sprintf('RSIerr_rel_b%d',bval))  = qcinfo(j).shells.brain_RSIerr.median_brain_rms_err_rel_bvals(b);
          end
          if isfield(qcinfo(j).shells,'brain_FWHM') && ...
             isfield(qcinfo(j).shells.brain_FWHM,'mean_brain_FWHMx')
            qcstruct(j).FWHMx = qcinfo(j).shells.brain_FWHM.mean_brain_FWHMx;
            qcstruct(j).FWHMy = qcinfo(j).shells.brain_FWHM.mean_brain_FWHMy;
            qcstruct(j).FWHMz = qcinfo(j).shells.brain_FWHM.mean_brain_FWHMz;
          end
          qcstruct(j).DTerr_norm  = qcinfo(j).shells.brain_DTerr.median_brain_rms_err_norm;
          qcstruct(j).DTerr_rel   = qcinfo(j).shells.brain_DTerr.median_brain_rms_err_rel;
          qcstruct(j).RSIerr_norm = qcinfo(j).shells.brain_RSIerr.median_brain_rms_err_norm;
          qcstruct(j).RSIerr_rel  = qcinfo(j).shells.brain_RSIerr.median_brain_rms_err_rel;
        end
        if isfield(qcinfo,'pepolar') && ~isempty(qcinfo(j).pepolar)
          qcstruct(j).pepolar_mi_diff = qcinfo(j).pepolar.mi_diff;
          qcstruct(j).pepolar_mi_for  = qcinfo(j).pepolar.mi_for;
          qcstruct(j).pepolar_mi_rev  = qcinfo(j).pepolar.mi_rev;
        end
    end % switch
  end % for
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_mask = create_brainmask(fname,parms,src_entities,src_suffix)
  % NOTE: the requested name has to carry parms.in_ext, because
  %   mmil_synthstrip_brainmask derives the mask name from the INPUT's
  %   extension; if the two disagree it writes one name and we look for another
  fname_mask_mgz = sprintf('%s/%s_brainmask%s',parms.tmpdir,get_fstem(fname),parms.in_ext);
  fname_mask = cmig_create_brainmask(fname,...
    'warning_flag',true,...
    'outdir',parms.tmpdir,...
    'fname_mask',fname_mask_mgz,...
    'binarize_flag',true,...
    'verbose',parms.verbose,...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [entities,suffix] = split_bids_fstem(fstem)
  tok = regexp(fstem,'^(.+)_([^_]+)$','tokens','once');
  if isempty(tok)
    entities = fstem;
    suffix = '';
  else
    entities = tok{1};
    suffix = tok{2};
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = stype_to_desc(stype)
  switch stype
    case {'bold','BOLD'}
      desc = 'bold';
    case {'dwi','DTI'}
      desc = 'dwi';
    case 'epi'
      desc = 'epi';
    otherwise
      desc = stype;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function task_label = get_task_label(fstem)
  % [^_]+ not \w+: \w includes the underscore, so task-mid_run-01 captured
  % 'mid_run'. That put _run into every func column prefix, left rest
  % unrecognized so its rows took the tfmri branch as tfmri__rest_run instead
  % of rsfmri, and pushed 22 per-run names past the 63-character field-name
  % limit. Matches get_dir_label and get_acq_entity, which already exclude
  % underscores
  tok = regexp(fstem,'(?:^|_)task-([^_]+)','tokens','once');
  if ~isempty(tok)
    task_label = lower(tok{1});
  else
    task_label = 'rest';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dir_label = get_dir_label(fstem)
  tok = regexp(fstem,'(?:^|_)dir-([^_]+)','tokens','once');
  if ~isempty(tok)
    dir_label = lower(tok{1});
  else
    dir_label = '';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = merge_session_fields(s, new_s)
  skip = {'participant_id','session_id'};
  flds = fieldnames(new_s);
  for k=1:length(flds)
    if ~any(strcmp(flds{k},skip))
      s.(flds{k}) = new_s.(flds{k});
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_session_rawqc(all_qcstruct, fname_tsv, parms)
  pid = parms.participant_id;
  sid = parms.session_id;
  session_row.participant_id = pid;
  session_row.session_id = sid;

  n = length(all_qcstruct);
  norm_stypes = cell(1,n);
  for i=1:n
    norm_stypes{i} = lower(mmil_getfield(all_qcstruct{i},'stype',''));
  end

  % sMRI T1w
  idx = find(strcmp(norm_stypes,'t1w'));
  if ~isempty(idx)
    sub = build_smri_session_row(all_qcstruct(idx),'mr_y_qc__raw__smri__t1',pid,sid);
    session_row = merge_session_fields(session_row,sub);
  end

  % sMRI T2w
  idx = find(strcmp(norm_stypes,'t2w'));
  if ~isempty(idx)
    sub = build_smri_session_row(all_qcstruct(idx),'mr_y_qc__raw__smri__t2',pid,sid);
    session_row = merge_session_fields(session_row,sub);
  end

  % bold: group by task
  idx_bold = find(strcmp(norm_stypes,'bold'));
  if ~isempty(idx_bold)
    bold_rows = all_qcstruct(idx_bold);
    task_labels = cellfun(@(r) get_task_label(mmil_getfield(r,'fstem','')), bold_rows,...
      'UniformOutput',false);
    unique_tasks = unique(task_labels);
    for t=1:length(unique_tasks)
      task = unique_tasks{t};
      task_rows = bold_rows(strcmp(task_labels,task));
      sub = build_bold_session_row(task_rows,task,pid,sid);
      session_row = merge_session_fields(session_row,sub);
    end
  end

  % dwi
  idx = find(ismember(norm_stypes,{'dwi','dti'}));
  if ~isempty(idx)
    sub = build_dwi_session_row(all_qcstruct(idx),pid,sid);
    session_row = merge_session_fields(session_row,sub);
  end

  % fmap (epi)
  idx = find(strcmp(norm_stypes,'epi'));
  if ~isempty(idx)
    sub = build_fmap_session_row(all_qcstruct(idx),pid,sid);
    session_row = merge_session_fields(session_row,sub);
  end

  write_rawqc_tsv({session_row},fname_tsv,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_smri_session_row(rows, prefix, pid, sid)
  n = length(rows);
  session_row.participant_id = pid;
  session_row.session_id = sid;
  % the series count is protocol compliance's, written by cmig_pc_bids_abcd

  brain_means = NaN(1,n); brain_stds = NaN(1,n);
  brainvols   = NaN(1,n); voxvols    = NaN(1,n);
  mean_map = [brain_snr_ext_colmap(); nonbrain_colmap(); fwhm_colmap()];
  nfm = size(mean_map,1);
  mean_accum = NaN(nfm,n);

  for i=1:n
    row = rows{i};
    v = mmil_getfield(row,'brain_mean',NaN);
    brain_means(i) = v;
    v = mmil_getfield(row,'brain_std',NaN);
    brain_stds(i)  = v;
    v = mmil_getfield(row,'brainvol',NaN);
    brainvols(i)   = v;
    v = mmil_getfield(row,'voxvol',NaN);
    voxvols(i)     = v;
    for m=1:nfm
      v = mmil_getfield(row,mean_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), mean_accum(m,i) = v; end
    end
  end

  session_row.(sprintf('%s__intensity__brain_mean',prefix)) = mean(brain_means,'omitnan');
  session_row.(sprintf('%s__intensity__brain_std',prefix))  = mean(brain_stds,'omitnan');
  session_row.(sprintf('%s__vol__brain_mean',prefix))       = mean(brainvols,'omitnan');
  session_row.(sprintf('%s__vol__vox_mean',prefix))         = mean(voxvols,'omitnan');
  for m=1:nfm
    session_row.(sprintf('%s__%s',prefix,mean_map{m,2})) = mean(mean_accum(m,:),'omitnan');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_bold_session_row(rows, task, pid, sid)
  switch lower(task)
    case 'rest'
      prefix = 'mr_y_qc__raw__rsfmri';
    otherwise
      prefix = sprintf('mr_y_qc__raw__tfmri__%s',lower(task));
  end
  n = length(rows);
  session_row.participant_id = pid;
  session_row.session_id = sid;
  % the series count is protocol compliance's, written by cmig_pc_bids_abcd

  % fields aggregated by mean across runs
  mean_map = {
    'brain_mean',            'intensity__brain_mean';
    'brain_std',             'intensity__brain_std';
    'brain_tSNR_mean',       'snr_mean';
    'brain_tSNR_median',     'snr_med';
    'brain_tSNR_std',        'snr_std';
    'mean_motion',           'mot_mean';
    'mean_motion_nody',      'mot__nody_mean';
    'mean_motion_filt',      'mot__filt_mean';
    'mean_motion_filt_nody', 'mot__filt__nody_mean';
    'fwhm_x',                'dtr__xdir_fwhm';
    'fwhm_y',                'dtr__ydir_fwhm';
    'fwhm_z',                'dtr__zdir_fwhm';
    'fwhm_raw_x',            'xdir_fwhm';
    'fwhm_raw_y',            'ydir_fwhm';
    'fwhm_raw_z',            'zdir_fwhm';
    'brainvol',              'vol__brain_mean';
    'voxvol',                'vol__vox_mean';
  };
  % fields aggregated by sum across runs (seconds of FD below threshold)
  sum_map = {
    'subthresh_02',           'fdsub__02_sum';    % 0.2 mm
    'subthresh_02_nody',      'fdsub__02__nody_sum';
    'subthresh_02_filt',      'fdsub__02__filt_sum';
    'subthresh_02_filt_nody', 'fdsub__02__filt__nody_sum';
    'subthresh_03',           'fdsub__03_sum';    % 0.3 mm [8.0]
    'subthresh_03_nody',      'fdsub__03__nody_sum';
    'subthresh_03_filt',      'fdsub__03__filt_sum';
    'subthresh_03_filt_nody', 'fdsub__03__filt__nody_sum';
    'subthresh_04',           'fdsub__04_sum';    % 0.4 mm [8.0]
    'subthresh_04_nody',      'fdsub__04__nody_sum';
    'subthresh_04_filt',      'fdsub__04__filt_sum';
    'subthresh_04_filt_nody', 'fdsub__04__filt__nody_sum';
  };
  mean_map = [mean_map; line_artifact_colmap(); brain_snr_ext_colmap(); ...
              motion_comp_colmap(); motion_comp_filt_colmap(); pepolar_colmap()];
  % fields aggregated by MAX across runs (worst run)
  max_map = [fov_cutoff_colmap(); motion_max_colmap(); ...
             motion_max_filt_colmap()];

  nfm = size(mean_map,1);
  nfs = size(sum_map,1);
  nfx = size(max_map,1);
  mean_accum = NaN(nfm,n);
  sum_accum  = NaN(nfs,n);
  max_accum  = NaN(nfx,n);

  for i=1:n
    row = rows{i};
    for m=1:nfm
      v = mmil_getfield(row,mean_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), mean_accum(m,i) = v; end
    end
    for m=1:nfs
      v = mmil_getfield(row,sum_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), sum_accum(m,i) = v; end
    end
    for m=1:nfx
      v = mmil_getfield(row,max_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), max_accum(m,i) = v; end
    end
  end

  % session-level aggregates
  for m=1:nfm
    session_row.(sprintf('%s__%s',prefix,mean_map{m,2})) = mean(mean_accum(m,:),'omitnan');
  end
  for m=1:nfs
    vals = sum_accum(m,~isnan(sum_accum(m,:)));
    if isempty(vals)
      session_row.(sprintf('%s__%s',prefix,sum_map{m,2})) = NaN;
    else
      session_row.(sprintf('%s__%s',prefix,sum_map{m,2})) = sum(vals);
    end
  end
  for m=1:nfx
    vals = max_accum(m,~isnan(max_accum(m,:)));
    if isempty(vals)
      session_row.(sprintf('%s__%s',prefix,max_map{m,2})) = NaN;
    else
      session_row.(sprintf('%s__%s',prefix,max_map{m,2})) = max(vals);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_dwi_session_row(rows, pid, sid)
  prefix = 'mr_y_qc__raw__dmri';
  n = length(rows);
  session_row.participant_id = pid;
  session_row.session_id = sid;
  % the series count is protocol compliance's, written by cmig_pc_bids_abcd

  mean_map = {
    'brain_mean',       'intensity__brain_mean';
    'brain_std',        'intensity__brain_std';
    'mean_motion',      'mot_mean';
    'mean_motion_nody', 'mot__nody_mean';
    'DTerr_mean',       'dterr_mean';
    'DTerr_median',     'dterr_med';
    'DTerr_std',        'dterr_std';
    'FA_mean',          'fa__brain_mean';
    'FA_median',        'fa__brain_med';
    'FA_std',           'fa__brain_std';
    'brainvol',         'vol__brain_mean';
    'voxvol',           'vol__vox_mean';
  };
  % censor counts: frames and frame-slice pairs sum across runs, slices
  %   average
  sum_map = {
    'nbad_frame_slices', 'censor__slice__frame_nbad';
    'nbad_frames',       'censor__frame_nbad';
  };
  censor_mean_map = {
    'nbad_slices',       'censor__slice_nbad';
  };
  mean_map = [mean_map; censor_mean_map; line_artifact_colmap(); ...
              brain_snr_ext_colmap(); motion_comp_colmap(); pepolar_colmap(); ...
              dwi_extra_mean_colmap()];
  sum_map = [sum_map; dwi_censor_sum_colmap()];
  % fields aggregated by MAX across runs (worst run)
  max_map = [fov_cutoff_colmap(); dwi_censor_max_colmap(); motion_max_colmap()];

  nfm = size(mean_map,1);
  nfs = size(sum_map,1);
  nfx = size(max_map,1);
  mean_accum  = NaN(nfm,n);
  sum_accum   = NaN(nfs,n);
  max_accum   = NaN(nfx,n);

  % collect unique b-values across all rows (from tSNR_b{bval} field names)
  bvals = [];
  for i=1:n
    flds = fieldnames(rows{i});
    for k=1:length(flds)
      tok = regexp(flds{k},'^tSNR_b(\d+)$','tokens','once');
      if ~isempty(tok)
        bv = str2double(tok{1});
        if ~any(bvals==bv), bvals(end+1) = bv; end %#ok<AGROW>
      end
    end
  end
  bvals = sort(bvals);
  nb = length(bvals);
  % per-bval columns driven by dwi_pershell_colmap (mixed mean/sum aggregation)
  pershell = dwi_pershell_colmap();
  np = size(pershell,1);
  bval_accum = NaN(nb*np,n);

  for i=1:n
    row = rows{i};
    for m=1:nfm
      v = mmil_getfield(row,mean_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), mean_accum(m,i) = v; end
    end
    for m=1:nfs
      v = mmil_getfield(row,sum_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), sum_accum(m,i) = v; end
    end
    for m=1:nfx
      v = mmil_getfield(row,max_map{m,1},NaN);
      if isnumeric(v) && isscalar(v), max_accum(m,i) = v; end
    end
    for b=1:nb
      bval = bvals(b);
      fi = (b-1)*np;
      for m=1:np
        fn  = sprintf(pershell{m,1},bval);
        sfx = sprintf(pershell{m,2},bval);
        v = mmil_getfield(row,fn,NaN);
        if isnumeric(v)&&isscalar(v), bval_accum(fi+m,i) = v; end
      end
    end
  end

  % session-level aggregates
  for m=1:nfm
    session_row.(sprintf('%s__%s',prefix,mean_map{m,2})) = mean(mean_accum(m,:),'omitnan');
  end
  for m=1:nfs
    vals = sum_accum(m,~isnan(sum_accum(m,:)));
    if isempty(vals)
      session_row.(sprintf('%s__%s',prefix,sum_map{m,2})) = NaN;
    else
      session_row.(sprintf('%s__%s',prefix,sum_map{m,2})) = sum(vals);
    end
  end
  for m=1:nfx
    vals = max_accum(m,~isnan(max_accum(m,:)));
    if isempty(vals)
      session_row.(sprintf('%s__%s',prefix,max_map{m,2})) = NaN;
    else
      session_row.(sprintf('%s__%s',prefix,max_map{m,2})) = max(vals);
    end
  end
  for b=1:nb
    bval = bvals(b);
    fi = (b-1)*np;
    for m=1:np
      sfx = sprintf(pershell{m,2},bval);
      rowvals = bval_accum(fi+m,:);
      if strcmp(pershell{m,3},'sum')
        vals = rowvals(~isnan(rowvals));
        if isempty(vals)
          session_row.(sprintf('%s__%s',prefix,sfx)) = NaN;
        else
          session_row.(sprintf('%s__%s',prefix,sfx)) = sum(vals);
        end
      else
        session_row.(sprintf('%s__%s',prefix,sfx)) = mean(rowvals,'omitnan');
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_fmap_session_row(rows, pid, sid)
  base_prefix = 'mr_y_qc__raw__';
  n = length(rows);
  session_row.participant_id = pid;
  session_row.session_id = sid;

  % group by acquisition and direction entities (e.g. func/dwi x ap/pa);
  % grouping by direction alone pools the fmaps acquired for fMRI with those
  % acquired for dMRI, which differ in resolution and smoothness, and collides
  % their run labels when both are run-01
  grp_labels = cell(1,n);
  grp_acqs   = cell(1,n);
  grp_dirs   = cell(1,n);
  for i=1:n
    fstem = mmil_getfield(rows{i},'fstem','');
    acq = acq_to_modality(get_acq_entity(fstem));
    dir_label = get_dir_label(fstem);
    grp_acqs{i} = acq;
    grp_dirs{i} = dir_label;
    parts = {};
    if ~isempty(acq), parts{end+1} = acq; end
    if ~isempty(dir_label), parts{end+1} = dir_label; end
    grp_labels{i} = strjoin(parts,'__');
  end
  unique_grps = unique(grp_labels);

  for d=1:length(unique_grps)
    grp = unique_grps{d};
    grp_idx = find(strcmp(grp_labels,grp));
    grp_prefix = fmap_group_prefix(base_prefix,grp_acqs{grp_idx(1)},...
                                   grp_dirs{grp_idx(1)});
    nd = length(grp_idx);
    % the field map series count is protocol compliance's, written by
    %   cmig_pc_bids_abcd

    brain_means = NaN(1,nd); brain_stds = NaN(1,nd);
    brainvols   = NaN(1,nd); voxvols    = NaN(1,nd);
    mean_map = [line_artifact_colmap(); brain_snr_ext_colmap(); fwhm_colmap()];
    nfm = size(mean_map,1);
    mean_accum = NaN(nfm,nd);
    max_map = fov_cutoff_colmap();
    nfx = size(max_map,1);
    max_accum = NaN(nfx,nd);

    % accumulate per scan, but write no r{NN}__ per-run columns: session-level
    % variables cmig adds are aggregates only. Every per-scan value is in the
    % series TSV, one row per scan, and the study date and the uids come from
    % cmig_pc_bids_abcd. This also removed 24 column names that exceeded
    % MATLAB's 63-character field-name limit once the acquisition token was
    % added to the prefix, e.g.
    % mr_y_qc__raw__fmri__fm__pa__r04__intensity__brain__maxvox_count
    for ii=1:nd
      row = rows{grp_idx(ii)};
      brain_means(ii) = mmil_getfield(row,'brain_mean',NaN);
      brain_stds(ii)  = mmil_getfield(row,'brain_std',NaN);
      brainvols(ii)   = mmil_getfield(row,'brainvol',NaN);
      voxvols(ii)     = mmil_getfield(row,'voxvol',NaN);
      for m=1:nfm
        v = mmil_getfield(row,mean_map{m,1},NaN);
        if isnumeric(v) && isscalar(v), mean_accum(m,ii) = v; end
      end
      for m=1:nfx
        v = mmil_getfield(row,max_map{m,1},NaN);
        if isnumeric(v) && isscalar(v), max_accum(m,ii) = v; end
      end
    end

    session_row.(sprintf('%s__intensity__brain_mean',grp_prefix)) = mean(brain_means,'omitnan');
    session_row.(sprintf('%s__intensity__brain_std',grp_prefix))  = mean(brain_stds,'omitnan');
    session_row.(sprintf('%s__vol__brain_mean',grp_prefix))       = mean(brainvols,'omitnan');
    session_row.(sprintf('%s__vol__vox_mean',grp_prefix))         = mean(voxvols,'omitnan');
    for m=1:nfm
      session_row.(sprintf('%s__%s',grp_prefix,mean_map{m,2})) = mean(mean_accum(m,:),'omitnan');
    end
    for m=1:nfx
      vals = max_accum(m,~isnan(max_accum(m,:)));
      if isempty(vals)
        session_row.(sprintf('%s__%s',grp_prefix,max_map{m,2})) = NaN;
      else
        session_row.(sprintf('%s__%s',grp_prefix,max_map{m,2})) = max(vals);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_str = get_run_entity(fstem)
  tok = regexp(fstem,'(?:^|_)run-(\d+)','tokens','once');
  if ~isempty(tok)
    run_str = tok{1};
  else
    run_str = 'n/a';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prefix = fmap_group_prefix(base_prefix,acq,dir_label)
  % the modality carries the fm token and the phase encode direction follows
  % it, as cmig_pc_bids_abcd spells the same grouping. A field map with no
  % acq entity has no modality to sit under and keeps the bare fmap stem
  if isempty(acq)
    prefix = sprintf('%sfmap',base_prefix);
  else
    prefix = sprintf('%s%s__fm',base_prefix,acq);
  end
  if ~isempty(dir_label)
    prefix = sprintf('%s__%s',prefix,dir_label);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modality = acq_to_modality(acq_str)
  % release-style modality label for a BIDS acq entity value, for use in
  % column names; file names keep the BIDS value
  switch lower(acq_str)
    case 'func'
      modality = 'fmri';
    case {'dwi','dti'}
      modality = 'dmri';
    case {'n/a',''}
      modality = '';
    otherwise
      modality = lower(acq_str);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acq_str = get_acq_entity(fstem)
  tok = regexp(fstem,'(?:^|_)acq-([^_]+)','tokens','once');
  if ~isempty(tok)
    acq_str = tok{1};
  else
    acq_str = 'n/a';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rec_str = get_rec_entity(fstem)
  tok = regexp(fstem,'(?:^|_)rec-([^_]+)','tokens','once');
  if ~isempty(tok)
    rec_str = tok{1};
  else
    rec_str = 'n/a';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_row = build_smri_series_row(row, prefix)
  series_row.participant_id = mmil_getfield(row,'participant_id','');
  series_row.session_id     = mmil_getfield(row,'session_id','');
  series_row.modality       = mmil_getfield(row,'stype','');
  fstem = mmil_getfield(row,'fstem','');
  series_row.acq_time = mmil_getfield(row,'StudyDate','');
  series_row.filename = sprintf('anat/%s.nii.gz', fstem);
  series_row.run            = get_run_entity(fstem);
  series_row.acq            = get_acq_entity(fstem);
  series_row.rec            = get_rec_entity(fstem);
  series_row.(sprintf('%s__intensity__brain_mean',prefix)) = mmil_getfield(row,'brain_mean',NaN);
  series_row.(sprintf('%s__intensity__brain_std',prefix))  = mmil_getfield(row,'brain_std',NaN);
  series_row.(sprintf('%s__vol__brain_mean',prefix))       = mmil_getfield(row,'brainvol',NaN);
  series_row.(sprintf('%s__vol__vox_mean',prefix))         = mmil_getfield(row,'voxvol',NaN);
  extra_map = [brain_snr_ext_colmap(); nonbrain_colmap(); fwhm_colmap()];
  for m=1:size(extra_map,1)
    series_row.(sprintf('%s__%s',prefix,extra_map{m,2})) = mmil_getfield(row,extra_map{m,1},NaN);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = line_artifact_colmap()
  % qcstruct field -> 7.0-style column suffix for the automated line-artifact
  % detectors; line__ = original detector, line__fft__ = fft detector. Aggregated
  % by mean across runs at the session level (like the other mean_map fields).
  map = {
    'line_max_score',     'auto__line__score_max';
    'line_mean_score',    'auto__line__score_mean';
    'line_max_count',     'auto__line__count_max';
    'line_mean_count',    'auto__line__count_mean';
    'line_nframes',       'auto__line__nframe_count';
    'line2_max_score',    'auto__line__fft__score_max';
    'line2_mean_score',   'auto__line__fft__score_mean';
    'line2_max_nslices',  'auto__line__fft__nslice_max';
    'line2_mean_nslices', 'auto__line__fft__nslice_mean';
    'line2_extent',       'auto__line__fft__extent_prop';
    'line2_nframes',      'auto__line__fft__nframe_count';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = line2_best_colmap()
  % the frame a multi-frame field map should be corrected from, and that
  % frame's fft line score. Series level only: the frame index is a property of
  % one scan and has no meaning averaged over the scans of a session
  map = {
    'line2_best_frame',   'auto__line__fft__frame_best';
    'line2_best_score',   'auto__line__fft__score_best';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = fwhm_colmap()
  % spatial smoothness per axis of the as-imported data ({axis}dir_fwhm, with
  % the average across frames left implicit); dtr__ marks the func values
  % measured after motion correction and detrending. Aggregated by mean across
  % runs at the session level
  map = {
    'fwhm_x', 'xdir_fwhm';
    'fwhm_y', 'ydir_fwhm';
    'fwhm_z', 'zdir_fwhm';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = brain_snr_ext_colmap()
  % extended within-brain intensity stats (brain_mean/std already emitted
  % elsewhere as intensity__brain_mean/std). Aggregated by mean across runs.
  map = {
    'brain_SNR',      'intensity__brain_snr';
    'brain_min',      'intensity__brain_min';
    'brain_max',      'intensity__brain_max';
    'brain_median',   'intensity__brain_med';
    'brain_nvox_max', 'intensity__brain__maxvox_count';
    'brain_fvox_max', 'intensity__brain__maxvox_prop';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = nonbrain_colmap()
  % intensity stats in the non-brain (background) region. Mean across runs.
  map = {
    'nonbrain_mean', 'intensity__nonbrain_mean';
    'nonbrain_std',  'intensity__nonbrain_std';
    'nonbrain_SNR',  'intensity__nonbrain_snr';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = fov_cutoff_colmap()
  % FOV cutoff scores (matches datadict/mr_y_qc__post__aut). Aggregated by MAX
  % across runs at the session level (worst run).
  map = {
    'cutoff',         'fov__cutoff_max';
    'dorsal_cutoff',  'fov__cutoff__dorsal_max';
    'ventral_cutoff', 'fov__cutoff__ventral_max';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = motion_comp_colmap()
  % translation/rotation motion components (mot_mean/nody already emitted).
  % Mean across runs. The maxima moved to motion_max_colmap
  map = {
    'mean_trans',      'mot__trans_mean';
    'mean_trans_nody', 'mot__trans__nody_mean';
    'mean_rot',        'mot__rot_mean';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = motion_max_colmap()
  % motion maxima, aggregated by MAX across runs, so the session value is the
  % worst excursion anywhere in the session rather than a mean of per-run
  % maxima, which the _max name would not describe. Note that the session name
  % is spelled the same as the per-run one, so nothing else signals which
  % aggregation was used
  map = {
    'max_motion',      'mot_max';
    'max_trans',       'mot__trans_max';
    'max_rot',         'mot__rot_max';
    'max_dx',          'mot__dx_max';
    'max_dy',          'mot__dy_max';
    'max_dz',          'mot__dz_max';
    'max_rx',          'mot__rx_max';
    'max_ry',          'mot__ry_max';
    'max_rz',          'mot__rz_max';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = motion_comp_filt_colmap()
  % respiratory-notch-filtered motion components (bold only). Mean across runs.
  % The maxima moved to motion_max_filt_colmap
  map = {
    'mean_trans_filt',      'mot__trans__filt_mean';
    'mean_trans_filt_nody', 'mot__trans__filt__nody_mean';
    'mean_rot_filt',        'mot__rot__filt_mean';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = motion_max_filt_colmap()
  % respiratory-notch-filtered motion maxima (bold only), aggregated by MAX
  % across runs, as the unfiltered maxima in motion_max_colmap are
  map = {
    'max_dx_filt',          'mot__dx__filt_max';
    'max_dy_filt',          'mot__dy__filt_max';
    'max_dz_filt',          'mot__dz__filt_max';
    'max_rx_filt',          'mot__rx__filt_max';
    'max_ry_filt',          'mot__ry__filt_max';
    'max_rz_filt',          'mot__rz__filt_max';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = pepolar_colmap()
  % pepolar (phase-encode reversal) mutual-information values. Mean across runs.
  map = {
    'pepolar_mi_for',  'pepolar__for_mi';
    'pepolar_mi_rev',  'pepolar__rev_mi';
    'pepolar_mi_diff', 'pepolar__diff_mi';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = dwi_extra_mean_colmap()
  % dwi non-FA tensor measures (FA already emitted) + whole-brain DTI/RSI fit
  % errors + raw-vs-corrected RMS diff + censoring threshold.
  % Aggregated by mean across runs.
  map = {
    'MD_mean',      'md__brain_mean';
    'MD_median',    'md__brain_med';
    'MD_std',       'md__brain_std';
    'b0_mean',      'intensity__b0_mean';
    'b0_median',    'intensity__b0_med';
    'b0_std',       'intensity__b0_std';
    'FWHMx',        'xdir_fwhm';
    'FWHMy',        'ydir_fwhm';
    'FWHMz',        'zdir_fwhm';
    'FWHMx_resid',  'dterr__xdir_fwhm';
    'FWHMy_resid',  'dterr__ydir_fwhm';
    'FWHMz_resid',  'dterr__zdir_fwhm';
    'DTerr_norm',   'dterr__norm_med';
    'DTerr_rel',    'dterr__rel_med';
    'RSIerr_norm',  'rsierr__norm_med';
    'RSIerr_rel',   'rsierr__rel_med';
    'diff_rms_mean',   'diff__rms__brain_mean';
    'diff_rms_median', 'diff__rms__brain_med';
    'diff_rms_std',    'diff__rms__brain_std';
    'censor_thresh','censor_thresh';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = dwi_censor_max_colmap()
  % dwi per-slice/per-frame worst-case censor counts.
  % Aggregated by MAX across runs.
  map = {
    'max_nbad_frames_per_slice',    'censor__frame__slice_max';
    'max_nbad_slices_per_frame',    'censor__slice__frame_max';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = dwi_censor_sum_colmap()
  % dwi censor totals and good-frame count. Summed across runs.
  map = {
    'ngood_frames',        'censor__frame_ngood';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function map = dwi_pershell_colmap()
  % per-b-value dwi columns. Columns: {field_fmt,
  % suffix_fmt, agg}; %d is the b-value. agg = session aggregation across runs.
  map = {
    'tSNR_b%d',              'snr__b%d_med',                  'mean';
    'FWHMx_b%d',             'b%d__xdir_fwhm',                'mean';
    'FWHMy_b%d',             'b%d__ydir_fwhm',                'mean';
    'FWHMz_b%d',             'b%d__zdir_fwhm',                'mean';
    'FWHMx_resid_b%d',       'b%d__dterr__xdir_fwhm',         'mean';
    'FWHMy_resid_b%d',       'b%d__dterr__ydir_fwhm',         'mean';
    'FWHMz_resid_b%d',       'b%d__dterr__zdir_fwhm',         'mean';
    'DTerr_norm_b%d',        'dterr__norm__b%d_med',          'mean';
    'DTerr_rel_b%d',         'dterr__rel__b%d_med',           'mean';
    'RSIerr_norm_b%d',       'rsierr__norm__b%d_med',         'mean';
    'RSIerr_rel_b%d',        'rsierr__rel__b%d_med',          'mean';
    'nframes_b%d',           'frame__b%d_count',              'sum';
    'nbad_frame_slices_b%d', 'censor__slice__frame__b%d_nbad','sum';
    'nbad_frames_b%d',       'censor__frame__b%d_nbad',       'sum';
    'nbad_slices_b%d',       'censor__slice__b%d_nbad',       'mean';
    'ngood_frames_b%d',      'censor__frame__b%d_ngood',      'sum';
  };
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_row = build_bold_series_row(row, prefix)
  series_row.participant_id = mmil_getfield(row,'participant_id','');
  series_row.session_id     = mmil_getfield(row,'session_id','');
  series_row.modality       = mmil_getfield(row,'stype','');
  fstem = mmil_getfield(row,'fstem','');
  series_row.acq_time = mmil_getfield(row,'StudyDate','');
  series_row.filename = sprintf('func/%s.nii.gz', fstem);
  series_row.run  = get_run_entity(fstem);
  series_row.task = get_task_label(fstem);
  series_row.acq  = get_acq_entity(fstem);
  series_row.rec  = get_rec_entity(fstem);
  mean_map = {
    'brain_mean',            'intensity__brain_mean';
    'brain_std',             'intensity__brain_std';
    'brain_tSNR_mean',       'snr_mean';
    'brain_tSNR_median',     'snr_med';
    'brain_tSNR_std',        'snr_std';
    'mean_motion',           'mot_mean';
    'mean_motion_nody',      'mot__nody_mean';
    'mean_motion_filt',      'mot__filt_mean';
    'mean_motion_filt_nody', 'mot__filt__nody_mean';
    'fwhm_x',                'dtr__xdir_fwhm';
    'fwhm_y',                'dtr__ydir_fwhm';
    'fwhm_z',                'dtr__zdir_fwhm';
    'fwhm_raw_x',            'xdir_fwhm';
    'fwhm_raw_y',            'ydir_fwhm';
    'fwhm_raw_z',            'zdir_fwhm';
    'brainvol',              'vol__brain_mean';
    'voxvol',                'vol__vox_mean';
  };
  sum_map = {
    'subthresh_02',           'fdsub__02_sum';
    'subthresh_02_nody',      'fdsub__02__nody_sum';
    'subthresh_02_filt',      'fdsub__02__filt_sum';
    'subthresh_02_filt_nody', 'fdsub__02__filt__nody_sum';
    'subthresh_03',           'fdsub__03_sum';
    'subthresh_03_nody',      'fdsub__03__nody_sum';
    'subthresh_03_filt',      'fdsub__03__filt_sum';
    'subthresh_03_filt_nody', 'fdsub__03__filt__nody_sum';
    'subthresh_04',           'fdsub__04_sum';
    'subthresh_04_nody',      'fdsub__04__nody_sum';
    'subthresh_04_filt',      'fdsub__04__filt_sum';
    'subthresh_04_filt_nody', 'fdsub__04__filt__nody_sum';
  };
  % series-level: aggregation type is irrelevant (one value per series), so all
  % extra metric groups can share the mean_map emission loop
  mean_map = [mean_map; line_artifact_colmap(); brain_snr_ext_colmap(); ...
              fov_cutoff_colmap(); motion_comp_colmap(); motion_max_colmap(); ...
              motion_comp_filt_colmap(); motion_max_filt_colmap(); ...
              pepolar_colmap()];
  for m=1:size(mean_map,1)
    series_row.(sprintf('%s__%s',prefix,mean_map{m,2})) = mmil_getfield(row,mean_map{m,1},NaN);
  end
  for m=1:size(sum_map,1)
    series_row.(sprintf('%s__%s',prefix,sum_map{m,2})) = mmil_getfield(row,sum_map{m,1},NaN);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_row = build_dwi_series_row(row, prefix)
  series_row.participant_id = mmil_getfield(row,'participant_id','');
  series_row.session_id     = mmil_getfield(row,'session_id','');
  series_row.modality       = mmil_getfield(row,'stype','');
  fstem_dwi = mmil_getfield(row,'fstem','');
  series_row.acq_time = mmil_getfield(row,'StudyDate','');
  series_row.filename = sprintf('dwi/%s.nii.gz', fstem_dwi);
  series_row.run            = get_run_entity(fstem_dwi);
  series_row.acq            = get_acq_entity(fstem_dwi);
  series_row.rec            = get_rec_entity(fstem_dwi);
  mean_map = {
    'brain_mean',       'intensity__brain_mean';
    'brain_std',        'intensity__brain_std';
    'mean_motion',      'mot_mean';
    'mean_motion_nody', 'mot__nody_mean';
    'DTerr_mean',       'dterr_mean';
    'DTerr_median',     'dterr_med';
    'DTerr_std',        'dterr_std';
    'FA_mean',          'fa__brain_mean';
    'FA_median',        'fa__brain_med';
    'FA_std',           'fa__brain_std';
    'brainvol',         'vol__brain_mean';
    'voxvol',           'vol__vox_mean';
  };
  censor_map = {
    'nbad_frame_slices', 'censor__slice__frame_nbad';
    'nbad_frames',       'censor__frame_nbad';
    'nbad_slices',       'censor__slice_nbad';
  };
  % series-level: aggregation type is irrelevant, so all static extra groups
  % share the mean_map emission loop
  mean_map = [mean_map; line_artifact_colmap(); brain_snr_ext_colmap(); ...
              fov_cutoff_colmap(); motion_comp_colmap(); motion_max_colmap(); ...
              pepolar_colmap(); dwi_extra_mean_colmap(); ...
              dwi_censor_max_colmap(); dwi_censor_sum_colmap()];
  for m=1:size(mean_map,1)
    series_row.(sprintf('%s__%s',prefix,mean_map{m,2})) = mmil_getfield(row,mean_map{m,1},NaN);
  end
  for m=1:size(censor_map,1)
    series_row.(sprintf('%s__%s',prefix,censor_map{m,2})) = mmil_getfield(row,censor_map{m,1},NaN);
  end
  % per-bval fields (dynamic; b-values vary by acquisition)
  pershell = dwi_pershell_colmap();
  flds = fieldnames(row);
  for k=1:length(flds)
    tok = regexp(flds{k},'^tSNR_b(\d+)$','tokens','once');
    if ~isempty(tok)
      bval = str2double(tok{1});
      for m=1:size(pershell,1)
        fn  = sprintf(pershell{m,1},bval);
        sfx = sprintf(pershell{m,2},bval);
        series_row.(sprintf('%s__%s',prefix,sfx)) = mmil_getfield(row,fn,NaN);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_row = build_fmap_series_row(row, prefix)
  series_row.participant_id = mmil_getfield(row,'participant_id','');
  series_row.session_id     = mmil_getfield(row,'session_id','');
  series_row.modality       = mmil_getfield(row,'stype','');
  fstem = mmil_getfield(row,'fstem','');
  series_row.acq_time = mmil_getfield(row,'StudyDate','');
  series_row.filename = sprintf('fmap/%s.nii.gz', fstem);
  series_row.run = get_run_entity(fstem);
  series_row.acq = get_acq_entity(fstem);
  series_row.rec = get_rec_entity(fstem);
  series_row.dir = get_dir_label(fstem);
  series_row.(sprintf('%s__intensity__brain_mean',prefix)) = mmil_getfield(row,'brain_mean',NaN);
  series_row.(sprintf('%s__intensity__brain_std',prefix))  = mmil_getfield(row,'brain_std',NaN);
  series_row.(sprintf('%s__vol__brain_mean',prefix))       = mmil_getfield(row,'brainvol',NaN);
  series_row.(sprintf('%s__vol__vox_mean',prefix))         = mmil_getfield(row,'voxvol',NaN);
  extra_map = [line_artifact_colmap(); line2_best_colmap(); ...
               brain_snr_ext_colmap(); fwhm_colmap(); fov_cutoff_colmap()];
  for m=1:size(extra_map,1)
    series_row.(sprintf('%s__%s',prefix,extra_map{m,2})) = mmil_getfield(row,extra_map{m,1},NaN);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_row = build_series_row(row, stype)
  % the series data is one table across the modalities, so every column takes
  % that table's name as its prefix. What the per-modality stem used to say is
  % in the row, which carries modality and task
  prefix = 'mr_y_qc__raw__series';
  switch lower(stype)
    case {'t1w','t2w'}
      series_row = build_smri_series_row(row,prefix);
    case 'bold'
      series_row = build_bold_series_row(row,prefix);
    case {'dwi','dti'}
      series_row = build_dwi_series_row(row,prefix);
    case 'epi'
      series_row = build_fmap_series_row(row,prefix);
    otherwise
      % other anat stypes (T2starw, FLAIR, etc.)
      series_row = build_smri_series_row(row,prefix);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_series_rawqc(stype_qcstruct, stype, fname_tsv, parms)
  series_rows = {};
  for i=1:length(stype_qcstruct)
    series_rows{end+1} = build_series_row(stype_qcstruct{i},stype);
  end
  write_rawqc_tsv(series_rows,fname_tsv,parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
