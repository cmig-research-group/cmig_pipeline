function [outdir,errcode,errmsg] = cmig_postqc_bids_dwi(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_postqc_bids_dwi(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Automated post-processing QC for BIDS DWI derivatives
%   Checks per-scan registration to T1w; computes session-level
%   brain cutoff and dMRI shell-specific metrics; copies the outlier slice
%   counts of the processed data's sidecar, per scan and per session
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%     processed data at {proc_rootdir}/{participant_id}/{session_id}/dwi
%     output from cmig_proc_bids_func (DWI processing path)
%   fsurf_rootdir: BIDS FreeSurfer derivatives root directory
%     FreeSurfer recons at {fsurf_rootdir}/{participant_id}_{session_id}
%     may be empty to skip FreeSurfer-dependent checks
%   output_rootdir: BIDS output root directory for QC results
%
% Optional Input ('key',value):
%   'min_ndirs': minimum number of diffusion directions (excluding b=0
%     frames) required for a DWI file to be included in dMRI
%     shell-specific metrics
%     {default = 10}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'verbose': [0|1] display status messages and warnings
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     form: {output_rootdir}/{participant_id}/{session_id}/dwi
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 07/14/2026 by Don Hagler
% Prev Mod: 08/04/2026 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/29/2026 by Don Hagler
% Prev Mod: 10/01/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/08/2026 by Don Hagler
% Prev Mod: 10/09/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/15/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,5), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin);
parms.status_extra = [];

[parms,errcode,errmsg] = check_input(parms);
if errcode
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-postqc',errcode,errmsg,parms.status_extra,startdate);
  return;
end

[parms,procflag] = check_output(parms);

if procflag
  [series_qcstruct,session_row,errcode,tmp_errmsg] = qc_dwi(parms);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  if ~errcode
    if ~isempty(series_qcstruct)
      fname_series = sprintf('%s/%s_desc-dwi_postqc.tsv',parms.outdir,parms.outstem);
      cmig_write_tabdat_postqc(series_qcstruct,fname_series,parms.forceflag);
    end
    if ~isempty(session_row)
      fname_session = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);
      cmig_write_tabdat_postqc({session_row},fname_session,parms.forceflag);
    end
  end
end

if parms.cleanup_flag
  cleanup_errmsg = cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir);
  errmsg = cmig_append_errmsg(errmsg,cleanup_errmsg);
end
% a tolerated failure from a resource error (e.g. out of memory) fails the
%   module, so the job is rerun
if ~errcode && cmig_is_resource_error(errmsg), errcode = 1; end
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-postqc',errcode,errmsg,parms.status_extra,startdate,procflag);

outdir = parms.outdir;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,options)
  options = cmig_parse_deployed_args(options);
  parms = mmil_args2parms(options,{...
    'participant_id',participant_id,[],...
    'session_id',session_id,[],...
    'proc_rootdir',proc_rootdir,[],...
    'fsurf_rootdir',fsurf_rootdir,[],...
    'output_rootdir',output_rootdir,[],...
  ...
    'min_ndirs',10,[2,1000],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });

  parms.outdir  = sprintf('%s/%s/%s/dwi',...
    parms.output_rootdir,parms.participant_id,parms.session_id);
  parms.tmpdir  = sprintf('%s/autoQC',parms.outdir);
  parms.outstem = sprintf('%s_%s',parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];
  pid = parms.participant_id;
  sid = parms.session_id;

  fprintf('%s: checking input for %s_%s (modality=dwi)\n',mfilename,pid,sid);

  parms.fnames_dwi = []; parms.fnames_info_dwi = [];
  [parms.fnames_dwi,parms.fnames_info_dwi] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'dwi','dwi');

  parms.fnames_t1w = []; parms.fnames_info_t1w = [];
  [parms.fnames_t1w,parms.fnames_info_t1w] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'anat','T1w');

  if isempty(parms.fsurf_rootdir)
    parms.fspath = [];
  else
    [parms.fspath,parms.fs_version,parms.fs_build,parms.fs_complete_flag,errcode,errmsg] = ...
      cmig_check_bids_fsurf(parms.participant_id,parms.session_id,parms.fsurf_rootdir);
    if errcode
      fprintf('%s: WARNING: %s\n',mfilename,errmsg);
      parms.fspath = [];
      errcode = 0;
    else
      parms.status_extra = struct('fs_version',parms.fs_version,...
        'fs_build',parms.fs_build,'fs_complete_flag',parms.fs_complete_flag);
    end
  end

  for f=1:length(parms.fnames_dwi)
    fprintf('%s: DWI processed: %s\n',mfilename,parms.fnames_dwi{f});
  end

  if isempty(parms.fnames_dwi)
    errcode = 1;
    errmsg = 'no processed DWI data found';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;
  output_files = {};

  output_files{end+1} = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);
  if ~isempty(parms.fnames_dwi)
    output_files{end+1} = sprintf('%s/%s_desc-dwi_postqc.tsv',parms.outdir,parms.outstem);
  end

  if ~exist(parms.outdir,'dir')
    procflag = 1;
  elseif parms.forceflag
    procflag = 1;
  else
    for f=1:length(output_files)
      if ~exist(output_files{f},'file')
        [~,name,ext] = fileparts(output_files{f});
        fprintf('%s: missing output file %s\n',mfilename,[name ext]);
        procflag = 1;
        break;
      end
    end
  end

  if procflag
    mmil_mkdir(parms.tmpdir);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [series_qcstruct,session_row,errcode,errmsg] = qc_dwi(parms)
  series_qcstruct = {}; session_row = []; errcode = 0; errmsg = [];

  n = length(parms.fnames_dwi);
  errcodes = zeros(n,1);
  errmsgs  = cell(n,1);
  raw_qcinfos = {};

  for f=1:n
    fname = parms.fnames_dwi{f};
    fname_info = parms.fnames_info_dwi{f};
    [qcinfo,errcodes(f),errmsgs{f}] = qc_dwi_scan(fname,fname_info,parms);
    if ~isempty(qcinfo)
      series_qcstruct{end+1} = convert_dwi_row(qcinfo,parms.participant_id,parms.session_id);
      raw_qcinfos{end+1} = qcinfo;
    end
  end

  for f=1:n
    errmsg = cmig_append_errmsg(errmsg,errmsgs{f});
  end

  if n > 0
    [base_session,errcode,tmp_errmsg] = calc_dMRI_session_qc(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end
    if ~isempty(base_session)
      % NOTE: b0_rms is a property of the field map pair, not of the dMRI
      %   data, so it is computed by cmig_postqc_bids_fmap
      b0_unwarp_indicator = cmig_get_b0_unwarp_indicator(parms.fnames_dwi,...
        parms.fnames_info_dwi,parms.proc_rootdir,parms.participant_id,...
        parms.session_id,'dwi');
      session_row = build_dwi_session_row(base_session,raw_qcinfos,...
        parms.participant_id,parms.session_id,b0_unwarp_indicator);
    end
  end

  if all(errcodes), errcode = 1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auto_qcinfo,errcode,errmsg] = qc_dwi_scan(fname,fname_info,parms)
  auto_qcinfo = []; errcode = 0; errmsg = [];

  if parms.verbose
    fprintf('%s: running QC for DWI scan: %s\n',mfilename,fname);
  end

  auto_qcinfo = cmig_init_epi_qcinfo(fname,'dwi',parms.participant_id,parms.session_id);
  [auto_qcinfo,errcode,errmsg] = cmig_read_proc_json_epi(auto_qcinfo,fname,fname_info);
  if errcode, return; end
  auto_qcinfo.outlier = read_outlier_info(fname_info);

  if ~isempty(parms.fnames_t1w) && ~isempty(auto_qcinfo.registration_matrix_T1)
    [~,fstem] = cmig_get_fileparts(fname);
    entities = regexprep(rm_desc_preproc(fstem),'_[^_]+$','');
    fname_plot = sprintf('%s/%s_regT1.jpg',parms.outdir,entities);
    [auto_qcinfo,errcode,errmsg] = cmig_check_regT1_epi(auto_qcinfo,fname,...
      'fnames_t1w',parms.fnames_t1w,'tmpdir',parms.tmpdir,'outdir',parms.outdir,...
      'fname_plot',fname_plot,...
      'fspath',parms.fspath,'verbose',parms.verbose,'forceflag',parms.forceflag);
    if errcode, return; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [session_row,errcode,errmsg] = calc_dMRI_session_qc(parms)
  session_row = []; errcode = 0; errmsg = [];
  n = length(parms.fnames_dwi);
  if n == 0, return; end

  if parms.verbose
    fprintf('%s: computing session-level DWI QC...\n',mfilename);
  end

  nscans_for = 0; nscans_rev = 0;
  for f=1:length(parms.fnames_dwi)
    fname_info = [];
    if f<=length(parms.fnames_info_dwi), fname_info = parms.fnames_info_dwi{f}; end
    direction = get_dwi_direction(parms.fnames_dwi{f},fname_info);
    switch direction
      case 'AP', nscans_for = nscans_for + 1;
      case 'PA', nscans_rev = nscans_rev + 1;
    end
  end

  qcinfo = [];
  qcinfo.VisitID                = sprintf('%s_%s',parms.participant_id,parms.session_id);
  qcinfo.stype                  = 'dwi_session';
  qcinfo.nscans                 = n;
  qcinfo.nscans_for              = nscans_for;
  qcinfo.nscans_rev              = nscans_rev;
  qcinfo.fstem                  = sprintf('%s_dwi_session',parms.outstem);
  qcinfo.brain_cutoff           = NaN;
  qcinfo.brain_cutoff_dorsal    = NaN;
  qcinfo.brain_cutoff_ventral   = NaN;
  qcinfo.FWHMx                  = NaN;
  qcinfo.FWHMy                  = NaN;
  qcinfo.FWHMz                  = NaN;
  qcinfo.FWHMx_resid            = NaN;
  qcinfo.FWHMy_resid            = NaN;
  qcinfo.FWHMz_resid            = NaN;
  qcinfo.DTerr_norm             = NaN;
  qcinfo.DTerr_rel              = NaN;
  qcinfo.RSIerr_norm            = NaN;
  qcinfo.RSIerr_rel             = NaN;

  [qcinfo,errcode,errmsg] = cmig_detect_brain_cutoff(qcinfo,...
    parms.fnames_dwi,parms.fnames_info_dwi,'DTI',...
    'tmpdir',parms.tmpdir,'outdir',parms.outdir,'outstem',parms.outstem,...
    'fname_plot',sprintf('%s/%s_braincutoff.jpg',parms.outdir,parms.outstem),...
    'fspath',parms.fspath,'verbose',parms.verbose,'forceflag',parms.forceflag);
  if errcode, return; end

  [qcinfo,errcode,tmp_errmsg] = check_dMRI_shells(qcinfo,parms);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  if errcode, return; end

  session_row = convert_dwi_row(qcinfo,parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qcinfo,errcode,errmsg] = check_dMRI_shells(qcinfo,parms)
  errcode = 0; errmsg = [];
  checkname = 'dMRI shell-specific metrics';
  n = length(parms.fnames_dwi);
  if n == 0, return; end

  fname_cache = sprintf('%s/DTI_dMRI_shells.mat',parms.tmpdir);
  if ~exist(fname_cache,'file') || parms.forceflag
    shells_tmpdir = sprintf('%s/DTI_dMRI_shells',parms.tmpdir);
    mmil_mkdir(shells_tmpdir);

    vol_all = []; qmat_all = []; bvals_all = []; M_all = [];
    for f=1:n
      fname = parms.fnames_dwi{f};

      fname_bvec = regexprep(fname,'\.nii.*','.bvec');
      fname_bval = regexprep(fname,'\.nii.*','.bval');
      if ~exist(fname_bvec,'file') || ~exist(fname_bval,'file')
        tmp_errmsg = sprintf('bvec/bval not found for %s, skipping',fname);
        fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
        errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
        continue;
      end
      bvec = load(fname_bvec);  % 3xN
      bval = load(fname_bval);  % 1xN or Nx1
      % count diffusion-weighted frames only, excluding b=0
      if sum(bval>0) < parms.min_ndirs, continue; end
      if size(bvec,1) ~= 3, bvec = bvec'; end
      bval = bval(:);
      qmat = bvec';  % Nx3

      fname_mgz = cmig_import_to_mgz(fname,shells_tmpdir,[],parms.forceflag);
      if isempty(fname_mgz)
        errmsg = sprintf('failed to import DWI file: %s',fname);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        errcode = 1;
        return;
      end
      [vol,M] = fs_load_mgh(fname_mgz);
      if isempty(M_all), M_all = M; end

      vol_all   = cat(4,vol_all,vol);
      qmat_all  = cat(1,qmat_all,qmat);
      bvals_all = cat(1,bvals_all,bval);
    end

    if isempty(vol_all)
      tmp_errmsg = sprintf('no DWI volumes loaded, skipping %s',checkname);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      return;
    end

    data = [];
    data.vol   = vol_all;
    data.M     = M_all;
    data.qmat  = qmat_all;
    data.bvals = bvals_all;

    if parms.verbose
      fprintf('%s: checking dMRI shell-specific metrics...\n',mfilename);
      tic;
    end
    dMRI_shells_metrics = autoqc_dMRI_shells([],...
      'data',data,...
      'outdir',parms.outdir,...
      'outstem',parms.outstem,...
      'tmpdir',shells_tmpdir,...
      'bids_flag',true,...
      'min_ndirs',parms.min_ndirs,...
      'verbose',parms.verbose,...
      'forceflag',parms.forceflag,...
      'cleanupflag',false);
    if parms.verbose, toc; end

    if isempty(dMRI_shells_metrics)
      tmp_errmsg = sprintf('autoqc_dMRI_shells returned empty, skipping %s',checkname);
      fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      return;
    end
    save(fname_cache,'dMRI_shells_metrics');
  else
    output = mmil_load_matfile(fname_cache,{'dMRI_shells_metrics'});
    dMRI_shells_metrics = output.dMRI_shells_metrics;
  end

  if isfield(dMRI_shells_metrics,'brain_FWHM') && ~isempty(dMRI_shells_metrics.brain_FWHM)
    qcinfo.FWHMx = mmil_getfield(dMRI_shells_metrics.brain_FWHM,'mean_brain_FWHMx',NaN);
    qcinfo.FWHMy = mmil_getfield(dMRI_shells_metrics.brain_FWHM,'mean_brain_FWHMy',NaN);
    qcinfo.FWHMz = mmil_getfield(dMRI_shells_metrics.brain_FWHM,'mean_brain_FWHMz',NaN);
  end
  % the fit-residual arm: autoqc_dMRI_shells measures the residual of the data
  % it was given, which for post-processing QC is the processed data
  if isfield(dMRI_shells_metrics,'brain_FWHM_resid') && ...
     ~isempty(dMRI_shells_metrics.brain_FWHM_resid)
    qcinfo.FWHMx_resid = mmil_getfield(dMRI_shells_metrics.brain_FWHM_resid,'mean_brain_FWHMx',NaN);
    qcinfo.FWHMy_resid = mmil_getfield(dMRI_shells_metrics.brain_FWHM_resid,'mean_brain_FWHMy',NaN);
    qcinfo.FWHMz_resid = mmil_getfield(dMRI_shells_metrics.brain_FWHM_resid,'mean_brain_FWHMz',NaN);
  end
  if isfield(dMRI_shells_metrics,'brain_DTerr') && ~isempty(dMRI_shells_metrics.brain_DTerr)
    qcinfo.DTerr_norm = mmil_getfield(dMRI_shells_metrics.brain_DTerr,'median_brain_rms_err_norm',NaN);
    qcinfo.DTerr_rel  = mmil_getfield(dMRI_shells_metrics.brain_DTerr,'median_brain_rms_err_rel',NaN);
  end
  if isfield(dMRI_shells_metrics,'brain_RSIerr') && ~isempty(dMRI_shells_metrics.brain_RSIerr)
    if isfield(dMRI_shells_metrics.brain_RSIerr,'median_brain_rms_err_norm')
      qcinfo.RSIerr_norm = mmil_getfield(dMRI_shells_metrics.brain_RSIerr,'median_brain_rms_err_norm',NaN);
      qcinfo.RSIerr_rel  = mmil_getfield(dMRI_shells_metrics.brain_RSIerr,'median_brain_rms_err_rel',NaN);
    end
  end

  for b=1:dMRI_shells_metrics.nbvals
    bval = dMRI_shells_metrics.unique_bvals(b);
    if isfield(dMRI_shells_metrics,'brain_FWHM') && ~isempty(dMRI_shells_metrics.brain_FWHM)
      qcinfo.(sprintf('FWHMx_b%d',bval)) = dMRI_shells_metrics.brain_FWHM.mean_brain_FWHMx_bvals(b);
      qcinfo.(sprintf('FWHMy_b%d',bval)) = dMRI_shells_metrics.brain_FWHM.mean_brain_FWHMy_bvals(b);
      qcinfo.(sprintf('FWHMz_b%d',bval)) = dMRI_shells_metrics.brain_FWHM.mean_brain_FWHMz_bvals(b);
    end
    if isfield(dMRI_shells_metrics,'brain_FWHM_resid') && ...
       ~isempty(dMRI_shells_metrics.brain_FWHM_resid)
      resid = dMRI_shells_metrics.brain_FWHM_resid;
      qcinfo.(sprintf('FWHMx_resid_b%d',bval)) = resid.mean_brain_FWHMx_bvals(b);
      qcinfo.(sprintf('FWHMy_resid_b%d',bval)) = resid.mean_brain_FWHMy_bvals(b);
      qcinfo.(sprintf('FWHMz_resid_b%d',bval)) = resid.mean_brain_FWHMz_bvals(b);
    end
    if isfield(dMRI_shells_metrics,'brain_tSNR') && ~isempty(dMRI_shells_metrics.brain_tSNR)
      qcinfo.(sprintf('tSNR_b%d',bval)) = dMRI_shells_metrics.brain_tSNR.median_brain_tSNR_bvals(b);
    end
    if isfield(dMRI_shells_metrics,'brain_DTerr') && ~isempty(dMRI_shells_metrics.brain_DTerr)
      qcinfo.(sprintf('DTerr_norm_b%d',bval)) = dMRI_shells_metrics.brain_DTerr.median_brain_rms_err_norm_bvals(b);
      qcinfo.(sprintf('DTerr_rel_b%d',bval))  = dMRI_shells_metrics.brain_DTerr.median_brain_rms_err_rel_bvals(b);
    end
    if isfield(dMRI_shells_metrics,'brain_RSIerr') && ~isempty(dMRI_shells_metrics.brain_RSIerr)
      if isfield(dMRI_shells_metrics.brain_RSIerr,'median_brain_rms_err_norm_bvals')
        qcinfo.(sprintf('RSIerr_norm_b%d',bval)) = dMRI_shells_metrics.brain_RSIerr.median_brain_rms_err_norm_bvals(b);
        qcinfo.(sprintf('RSIerr_rel_b%d',bval))  = dMRI_shells_metrics.brain_RSIerr.median_brain_rms_err_rel_bvals(b);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function direction = get_dwi_direction(fname,fname_info)
  direction = '';  % empty char (not []) so 'switch direction' is safe when unknown
  [~,fstem] = cmig_get_fileparts(fname);
  tok = regexp(fstem,'(?:^|_)dir-([^_]+)','tokens','once');
  if ~isempty(tok)
    val = upper(tok{1});
    if ismember(val,{'AP','PA'})
      direction = val;
      return;
    end
  end
  if isempty(fname_info) || ~exist(fname_info,'file')
    fname_info = regexprep(fname,'\.nii.*','.json');
  end
  if ~exist(fname_info,'file'), return; end
  try
    info = jsondecode(fileread(fname_info));
    post = pe_toward_posterior(info);
    if isempty(post), return; end
    if post
      direction = 'AP';
    else
      direction = 'PA';
    end
  catch
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post = pe_toward_posterior(info)
  % true when the phase encoding runs toward posterior (AP): the
  %   PhaseEncodingDirection axis and sign through the image's vox2ras
  %   (Affine), so right in any output orientation; without an Affine, or for
  %   an axis with no A/P component, the sign as in dcm2niix's LAS, where j
  %   points anterior
  post = [];
  pe = mmil_getfield(info,'PhaseEncodingDirection','');
  if isempty(pe) || ~ischar(pe), return; end
  a = strfind('ijk',pe(1));
  if isempty(a), return; end
  s = 1;
  if any(pe=='-'), s = -1; end
  M = mmil_getfield(info,'Affine',[]);
  if isnumeric(M) && isequal(size(M),[4 4]) && ...
      abs(M(2,a)) > 0.5*norm(M(1:3,a))
    post = s*M(2,a) < 0;
  else
    post = s < 0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = convert_dwi_row(qcinfo, participant_id, session_id)
  row.participant_id = participant_id;
  row.session_id     = session_id;
  switch qcinfo.stype
    case 'dwi'
      % series-identifying entity columns
      fstem = mmil_getfield(qcinfo,'fstem','');
      row.modality = 'dwi';
      row.acq_time = qcinfo.StudyDate;
      fstem = rm_desc_preproc(fstem);
      row.filename = sprintf('dwi/%s.nii.gz', fstem);
      tok = regexp(fstem,'(?:^|_)run-([^_]+)','tokens','once');
      if ~isempty(tok), row.run = tok{1}; else row.run = 'n/a'; end
      tok = regexp(fstem,'(?:^|_)acq-([^_]+)','tokens','once');
      if ~isempty(tok), row.acq = tok{1}; else row.acq = 'n/a'; end
      tok = regexp(fstem,'(?:^|_)rec-([^_]+)','tokens','once');
      if ~isempty(tok), row.rec = tok{1}; else row.rec = 'n/a'; end
      % series columns are mr_y_qc__post__series__{measure}, with no modality
      row.mr_y_qc__post__series__regt1__imgcorr_score  = qcinfo.regT1_corr;
      row.mr_y_qc__post__series__rigid_score = qcinfo.regT1_norm;
      row.mr_y_qc__post__series__warp_score  = qcinfo.regT1_warp;
      row.mr_y_qc__post__series__regt1__imgcorr1_score = qcinfo.regT1_imgcorr1;
      row.mr_y_qc__post__series__regt1__imgcorr2_score = qcinfo.regT1_imgcorr2;
      row.mr_y_qc__post__series__regt1__imgdist1_score = qcinfo.regT1_imgdist3;
      row.mr_y_qc__post__series__regt1__imgdist2_score = qcinfo.regT1_imgdist4;
      row.mr_y_qc__post__series__regt1__segreg__cost_score = qcinfo.regT1_cost0;
      row = add_outlier_cols(row,'mr_y_qc__post__series__outlier__',...
        mmil_getfield(qcinfo,'outlier',[]));
    case 'dwi_session'
      row.mr_y_qc__post__aut__dmri__nscans_count              = qcinfo.nscans;
      row.mr_y_qc__post__aut__dmri__nscans__for_count          = qcinfo.nscans_for;
      row.mr_y_qc__post__aut__dmri__nscans__rev_count          = qcinfo.nscans_rev;
      row.mr_y_qc__post__aut__dmri__brain__cutoff_score        = qcinfo.brain_cutoff;
      row.mr_y_qc__post__aut__dmri__brain__cutoff__dorsal_score = qcinfo.brain_cutoff_dorsal;
      row.mr_y_qc__post__aut__dmri__brain__cutoff__ventral_score = qcinfo.brain_cutoff_ventral;
      row.mr_y_qc__post__aut__dmri__xdir_fwhm            = qcinfo.FWHMx;
      row.mr_y_qc__post__aut__dmri__ydir_fwhm            = qcinfo.FWHMy;
      row.mr_y_qc__post__aut__dmri__zdir_fwhm            = qcinfo.FWHMz;
      row.mr_y_qc__post__aut__dmri__dterr__xdir_fwhm     = qcinfo.FWHMx_resid;
      row.mr_y_qc__post__aut__dmri__dterr__ydir_fwhm     = qcinfo.FWHMy_resid;
      row.mr_y_qc__post__aut__dmri__dterr__zdir_fwhm     = qcinfo.FWHMz_resid;
      row.mr_y_qc__post__aut__dmri__dterr__norm_med            = qcinfo.DTerr_norm;
      row.mr_y_qc__post__aut__dmri__dterr__rel_med             = qcinfo.DTerr_rel;
      row.mr_y_qc__post__aut__dmri__rsierr__norm_med           = qcinfo.RSIerr_norm;
      row.mr_y_qc__post__aut__dmri__rsierr__rel_med            = qcinfo.RSIerr_rel;
      % per-bval fields: rename to 7.0-style
      pfx = 'mr_y_qc__post__aut__dmri__';
      flds = fieldnames(qcinfo);
      bval_flds = flds(~cellfun(@isempty, regexp(flds,'^(FWHM[xyz](_resid)?|tSNR|DTerr_norm|DTerr_rel|RSIerr_norm|RSIerr_rel)_b\d+')));
      for k = 1:length(bval_flds)
        fld = bval_flds{k};
        new_fld = regexprep(fld,'^FWHM([xyz])_resid_b(\d+)',[pfx 'b$2__dterr__$1dir_fwhm']);
        new_fld = regexprep(new_fld,'^FWHM([xyz])_b(\d+)',[pfx 'b$2__$1dir_fwhm']);
        new_fld = regexprep(new_fld,'^tSNR_b(\d+)',[pfx 'temp__snr__b$1_med']);
        new_fld = regexprep(new_fld,'^DTerr_(norm|rel)_b(\d+)',[pfx 'dterr__$1__b$2_med']);
        new_fld = regexprep(new_fld,'^RSIerr_(norm|rel)_b(\d+)',[pfx 'rsierr__$1__b$2_med']);
        row.(new_fld) = qcinfo.(fld);
      end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_dwi_session_row(existing_row,raw_qcinfos,...
  participant_id,session_id,b0_unwarp_indicator)
  session_row = struct();
  session_row.participant_id = participant_id;
  session_row.session_id     = session_id;
  skip = {'participant_id','session_id'};
  flds = fieldnames(existing_row);
  for k=1:length(flds)
    if ~ismember(flds{k},skip)
      session_row.(flds{k}) = existing_row.(flds{k});
    end
  end
  session_row.mr_y_qc__post__aut__dmri__b0__unwarp_indicator = b0_unwarp_indicator;
  regT1_src  = {'regT1_corr','regT1_norm','regT1_warp',...
                'regT1_imgcorr1','regT1_imgcorr2','regT1_imgdist3','regT1_imgdist4',...
                'regT1_cost0'};
  regT1_dest = {'mr_y_qc__post__aut__dmri__regt1__imgcorr_score',...
                'mr_y_qc__post__aut__dmri__rigid_score',...
                'mr_y_qc__post__aut__dmri__warp_score',...
                'mr_y_qc__post__aut__dmri__regt1__imgcorr1_score',...
                'mr_y_qc__post__aut__dmri__regt1__imgcorr2_score',...
                'mr_y_qc__post__aut__dmri__regt1__imgdist1_score',...
                'mr_y_qc__post__aut__dmri__regt1__imgdist2_score',...
                'mr_y_qc__post__aut__dmri__regt1__segreg__cost_score'};
  for k=1:length(regT1_src)
    vals = cellfun(@(q) mmil_getfield(q,regT1_src{k},NaN), raw_qcinfos);
    vals = vals(~isnan(vals));
    if ~isempty(vals)
      session_row.(regT1_dest{k}) = mean(vals);
    else
      session_row.(regT1_dest{k}) = NaN;
    end
  end
  % outlier slices across runs: counts summed, maxima the highest
  outliers = cellfun(@(q) mmil_getfield(q,'outlier',[]),raw_qcinfos,...
    'UniformOutput',false);
  session_row = add_outlier_cols(session_row,...
    'mr_y_qc__post__aut__dmri__outlier__',combine_outlier_info(outliers));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flds = outlier_fields()
  % the outlier slice detection fields of the cmig_proc_bids_dwi sidecar:
  %   {sidecar field, column tail, 1 for a maximum or 0 for a count}
  flds = {'OutlierVolumeCount','frame_nbad',0;...
          'OutlierSliceCount','slice_nbad',0;...
          'OutlierSliceVolumeCount','slice__frame_nbad',0;...
          'MaxOutlierVolumesPerSlice','frame__slice_max',1;...
          'MaxOutlierSlicesPerVolume','slice__frame_max',1};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outlier = read_outlier_info(fname_info)
  % the run totals and the per-shell values (with Shells, the b-values) of
  %   the outlier fields, or empty where the sidecar has none
  outlier = [];
  if isempty(fname_info) || ~exist(fname_info,'file'), return; end
  info = jsondecode(fileread(fname_info));
  flds = outlier_fields();
  if ~all(isfield(info,flds(:,1))), return; end
  outlier.shells = [];
  if isfield(info,'Shells')
    outlier.shells = round(info.Shells(:)');
  end
  for k=1:size(flds,1)
    outlier.(flds{k,1}) = double(info.(flds{k,1}));
    fld_shell = [flds{k,1} 'PerShell'];
    vals = [];
    if isfield(info,fld_shell), vals = double(info.(fld_shell)(:)'); end
    if length(vals) ~= length(outlier.shells), vals = []; end
    outlier.(fld_shell) = vals;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outlier = combine_outlier_info(outliers)
  % one session's outlier values from its runs': counts summed, maxima the
  %   highest, per shell by b-value
  outlier = [];
  outliers = outliers(~cellfun(@isempty,outliers));
  if isempty(outliers), return; end
  flds = outlier_fields();
  outlier.shells = unique(cell2mat(cellfun(@(o) o.shells,outliers,...
    'UniformOutput',false)));
  for k=1:size(flds,1)
    vals = cellfun(@(o) o.(flds{k,1}),outliers);
    fld_shell = [flds{k,1} 'PerShell'];
    shell_vals = nan(length(outliers),length(outlier.shells));
    for r=1:length(outliers)
      [~,i_out,i_run] = intersect(outlier.shells,outliers{r}.shells);
      if ~isempty(outliers{r}.(fld_shell))
        shell_vals(r,i_out) = outliers{r}.(fld_shell)(i_run);
      end
    end
    if flds{k,3}
      outlier.(flds{k,1}) = max(vals);
      outlier.(fld_shell) = max(shell_vals,[],1);
    else
      outlier.(flds{k,1}) = sum(vals);
      outlier.(fld_shell) = sum(shell_vals,1,'omitnan');
      outlier.(fld_shell)(all(isnan(shell_vals),1)) = NaN;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = add_outlier_cols(row,pfx,outlier)
  % {pfx}{tail} for the run totals and {pfx}{tail} with __b{b} before its
  %   last component for each shell, as cmig_rawqc_bids names its censoring
  if isempty(outlier), return; end
  flds = outlier_fields();
  for k=1:size(flds,1)
    tail = flds{k,2};
    row.([pfx tail]) = outlier.(flds{k,1});
    vals = outlier.([flds{k,1} 'PerShell']);
    stat = regexprep(tail,'^.*_','');
    stem = tail(1:end-length(stat)-1);
    for b=1:length(vals)
      row.(sprintf('%s%s__b%d_%s',pfx,stem,outlier.shells(b),stat)) = vals(b);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = rm_desc_preproc(fstem)
  % the rawdata name of a proc file, which adds desc-preproc: the filename
  %   column joins with raw QC on it, and the jpg names carry the run entities
  fstem = regexprep(fstem,'_desc-preproc(?=_[^_]+$)','');
return;
