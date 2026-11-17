function [outdir,errcode,errmsg] = cmig_postqc_bids_func(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_postqc_bids_func(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Automated post-processing QC for BIDS func derivatives
%   Computes session-level fMRI stdev/tSNR and brain cutoff, and each
%   run's correlation and centre of mass distance to the other runs; the
%   registration to T1w is checked on the fMRI field maps, by
%   cmig_postqc_bids_fmap
%   Each run's row also carries, from its proc sidecar, the field map pair
%   used for B0 unwarping (fallback, shim match, displacement, coverage)
%   and the start of its between-scan registration (cmig_proc_bids_func)
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%     processed data at {proc_rootdir}/{participant_id}/{session_id}/func
%     output from cmig_proc_bids_func
%   fsurf_rootdir: BIDS FreeSurfer derivatives root directory
%     FreeSurfer recons at {fsurf_rootdir}/{participant_id}_{session_id}
%     may be empty to skip FreeSurfer-dependent checks
%   output_rootdir: BIDS output root directory for QC results
%
% Optional Input ('key',value):
%   'skipTRs': number of initial TRs to omit from stdev calculation
%     if empty, set automatically from scan JSON via cmig_set_abcd_skipTRs
%     {default = []}
%   'skipTRs_GE_DV25': skipTRs for GE DV25/Orchestra
%     {default = 4}
%   'skipTRs_GE': skipTRs for GE other software
%     {default = 15}
%   'skipTRs_Philips': skipTRs for Philips
%     {default = 8}
%   'skipTRs_Siemens': skipTRs for Siemens
%     {default = 8}
%   'fMRI_topval_prctile': prctile for scaled fMRI image (manual QC)
%     {default = 99.9}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'verbose': [0|1] display status messages and warnings
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     form: {output_rootdir}/{participant_id}/{session_id}/func
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Prev Mod: 06/12/2026 by Don Hagler
% Prev Mod: 06/15/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 07/22/2026 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/24/2026 by Don Hagler
% Prev Mod: 09/29/2026 by Don Hagler
% Prev Mod: 10/01/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/20/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/12/2026 by Don Hagler
% Prev Mod: 11/13/2026 by Don Hagler
% Last Mod: 11/15/2026 by Don Hagler
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
  [series_qcstruct,session_row,errcode,errmsg] = qc_func(parms);
  if ~errcode
    if ~isempty(series_qcstruct)
      fname_series = sprintf('%s/%s_desc-bold_postqc.tsv',parms.outdir,parms.outstem);
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
  if isempty(errmsg), errmsg = cleanup_errmsg; end
end
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
    'skipTRs',[],[0,Inf],...
    'skipTRs_GE_DV25',4,[0,Inf],...
    'skipTRs_GE',15,[0,Inf],...
    'skipTRs_Philips',8,[0,Inf],...
    'skipTRs_Siemens',8,[0,Inf],...
    'fMRI_topval_prctile',99.9,[90,100],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });

  parms.outdir  = sprintf('%s/%s/%s/func',...
    parms.output_rootdir,parms.participant_id,parms.session_id);
  parms.tmpdir  = sprintf('%s/autoQC',parms.outdir);
  parms.outstem = sprintf('%s_%s',parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];
  pid = parms.participant_id;
  sid = parms.session_id;

  fprintf('%s: checking input for %s_%s (modality=func)\n',mfilename,pid,sid);

  parms.fnames_func = []; parms.fnames_info_func = [];
  [parms.fnames_func,parms.fnames_info_func] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'func','bold');

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

  for f=1:length(parms.fnames_func)
    fprintf('%s: func processed: %s\n',mfilename,parms.fnames_func{f});
  end

  if isempty(parms.fnames_func)
    errcode = 1;
    errmsg = 'no processed func data found';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;
  output_files = {};

  output_files{end+1} = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);
  if ~isempty(parms.fnames_func)
    output_files{end+1} = sprintf('%s/%s_desc-bold_postqc.tsv',parms.outdir,parms.outstem);
  end
  % the meanNorm sidecars, for the images that exist
  for fstem = {'desc-meanNorm_boldref','space-T1w_desc-meanNorm_boldref'}
    fname_img = sprintf('%s/%s_%s.nii.gz',parms.outdir,parms.outstem,fstem{1});
    if exist(fname_img,'file')
      output_files{end+1} = regexprep(fname_img,'\.nii\.gz$','.json');
    end
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

function [series_qcstruct,merged_session,errcode,errmsg] = qc_func(parms)
  series_qcstruct = {}; merged_session = []; errcode = 0; errmsg = [];

  n = length(parms.fnames_func);
  errcodes = zeros(n,1);
  errmsgs  = cell(n,1);
  scan_indices = [];

  stdev_row = []; scan_rms_vals = nan(n,1); fname_brain = [];
  scan_corr_vals = nan(n,1); scan_com_dists = nan(n,1);
  if n > 0
    [stdev_row,scan_rms_vals,scan_corr_vals,scan_com_dists,fname_brain,errcode,errmsg] = ...
      calc_fMRI_stdev(parms);
    if errcode, return; end
  end

  for f=1:n
    fname = parms.fnames_func{f};
    fname_info = parms.fnames_info_func{f};
    [qcinfo,errcodes(f),errmsgs{f}] = qc_func_scan(fname,fname_info,parms);
    if ~isempty(qcinfo)
      row = convert_func_row(qcinfo,parms.participant_id,parms.session_id);
      series_qcstruct{end+1} = add_proc_info(row,fname_info);
      scan_indices(end+1) = f;
    end
  end

  cutoff_row = [];
  if n > 0
    [cutoff_row,errcode,errmsg] = calc_fMRI_session_qc(parms);
    if errcode, return; end
  end

  % augment series rows with per-scan RMS of diff from inter-scan mean
  for k=1:length(series_qcstruct)
    f = scan_indices(k);
    series_qcstruct{k}.mr_y_qc__post__series__inter__diff_rms = scan_rms_vals(f);
    series_qcstruct{k}.mr_y_qc__post__series__inter__corr_score = scan_corr_vals(f);
    series_qcstruct{k}.mr_y_qc__post__series__inter__com_dist = scan_com_dists(f);
  end

  % spatial smoothness of each processed run, two arms: the processed image
  % itself (unmarked, no detrending, comparable to the raw QC value) and its
  % detrended residual (__dtr__, comparable to the raw QC detrended value)
  [fwhm_proc,fwhm_dt] = calc_func_fwhm(parms,fname_brain);
  axes_names = {'x','y','z'};
  for k=1:length(series_qcstruct)
    f = scan_indices(k);
    for a=1:3
      series_qcstruct{k}.(sprintf('mr_y_qc__post__series__%sdir_fwhm',...
        axes_names{a})) = fwhm_proc(f,a);
      series_qcstruct{k}.(sprintf('mr_y_qc__post__series__dtr__%sdir_fwhm',...
        axes_names{a})) = fwhm_dt(f,a);
    end
  end

  % NOTE: b0_rms is a property of the field map pair, not of the func data,
  %   so it is computed by cmig_postqc_bids_fmap
  b0_unwarp_indicator = cmig_get_b0_unwarp_indicator(parms.fnames_func,...
    parms.fnames_info_func,parms.proc_rootdir,parms.participant_id,...
    parms.session_id,'func');
  fmap_fallback_indicator = cmig_get_fmap_fallback_indicator(parms.fnames_func,...
    parms.fnames_info_func);

  % always build the session row, even if the stdev or cutoff step yielded
  %   nothing, so the session TSV is written and the module does not
  %   reprocess indefinitely
  merged_session = build_func_session_row(stdev_row,cutoff_row,fwhm_proc,fwhm_dt,...
    parms.participant_id,parms.session_id,b0_unwarp_indicator,...
    fmap_fallback_indicator);

  nonempty_msgs = errmsgs(~cellfun(@isempty,errmsgs));
  if ~isempty(nonempty_msgs)
    errmsg = strtrim(sprintf('%s ',nonempty_msgs{:}));
  end
  if all(errcodes), errcode = 1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auto_qcinfo,errcode,errmsg] = qc_func_scan(fname,fname_info,parms)
  auto_qcinfo = []; errcode = 0; errmsg = [];

  if parms.verbose
    fprintf('%s: running QC for func scan: %s\n',mfilename,fname);
  end

  auto_qcinfo = cmig_init_epi_qcinfo(fname,'func',parms.participant_id,parms.session_id);
  [auto_qcinfo,errcode,errmsg] = cmig_read_proc_json_epi(auto_qcinfo,fname,fname_info);
  if errcode, return; end
  % no registration check: at steady state the gradient-echo runs have too
  %   little gray/white contrast; cmig_postqc_bids_fmap checks the field maps
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stdev_row,scan_rms_vals,scan_corr_vals,scan_com_dists,fname_brain,errcode,errmsg] = calc_fMRI_stdev(parms)
  stdev_row = []; scan_rms_vals = nan(length(parms.fnames_func),1);
  scan_corr_vals = nan(length(parms.fnames_func),1);
  scan_com_dists = nan(length(parms.fnames_func),1);
  fname_brain = []; errcode = 0; errmsg = [];
  checkname = 'fMRI stdev and tSNR';

  if parms.verbose
    fprintf('%s: calculating fMRI stdev and tSNR...\n',mfilename);
    tic;
  end

  stdev_tmpdir = sprintf('%s/func_stdev',parms.tmpdir);
  mmil_mkdir(stdev_tmpdir);
  fnamelist_mgz = {};
  for f=1:length(parms.fnames_func)
    fname_mgz = cmig_import_to_mgz(parms.fnames_func{f},stdev_tmpdir,[],parms.forceflag);
    if isempty(fname_mgz)
      errmsg = sprintf('failed to import func file: %s',parms.fnames_func{f});
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
    fnamelist_mgz{end+1} = fname_mgz;
  end

  % the stdev/tSNR brain mask is a per-scan SynthStrip mask that
  % mmil_intrascan_stdev builds (cmig_create_brainmask, from the inter-scan
  % mean); fname_brain is set to that mask after the group loop below

  % split scans by phase-encode direction; forward/single-direction ->
  % base columns, reverse (PhaseEncodingDirection contains '-') -> _rev.
  % For a single-direction study (e.g. ABCD) the reverse group is empty
  % and the _rev columns are all NaN.
  revflags = get_func_revflags(parms);
  tag_list = {'inter_mean','inter_stdev','intra_stdev',...
              'inter_tSNR','intra_tSNR','inter_cov','intra_cov',...
              'inter_intra_stdev_ratio','inter_intra_tSNR_ratio',...
              'inter_tSNR_median','intra_tSNR_median','inter_intra_tSNR_median_ratio'};

  p = [];
  p.outstem   = 'func';
  p.skipTRs   = get_skipTRs(parms);
  p.plotflag  = 1;
  p.verbose   = parms.verbose;
  p.forceflag = parms.forceflag;

  % group specs: base (forward/single-direction) always; _rev when present
  groups = struct('suffix',{''},...
    'outdir',{sprintf('%s/func_stdev_results',parms.tmpdir)},'idx',{find(~revflags)});
  if any(revflags)
    groups(end+1) = struct('suffix','_rev',...
      'outdir',sprintf('%s/func_stdev_results_rev',parms.tmpdir),'idx',find(revflags));
  end
  stdev_outdir = groups(1).outdir;  % base (forward) group; used for meanNorm below

  % initialize qcinfo with base + _rev tags (NaN until a group fills them)
  stdev_qcinfo = [];
  stdev_qcinfo.VisitID = sprintf('%s_%s',parms.participant_id,parms.session_id);
  stdev_qcinfo.stype   = 'func_stdev';
  stdev_qcinfo.nscans  = length(fnamelist_mgz);
  % per-polarity scan counts, as postqc dwi reports for its scans
  stdev_qcinfo.nscans_for = nnz(~revflags);
  stdev_qcinfo.nscans_rev = nnz(revflags);
  stdev_qcinfo.fstem   = sprintf('%s_func_stdev',parms.outstem);
  for t=1:length(tag_list)
    stdev_qcinfo.(tag_list{t}) = NaN;
    stdev_qcinfo.(sprintf('%s_rev',tag_list{t})) = NaN;
  end
  align_tags = {'inter_corr_min','inter_com_dist_max'};
  for t=1:length(align_tags)
    stdev_qcinfo.(align_tags{t}) = NaN;
    stdev_qcinfo.(sprintf('%s_rev',align_tags{t})) = NaN;
  end

  any_results = false;
  for g=1:length(groups)
    gidx = groups(g).idx;
    if isempty(gidx), continue; end
    p.outdir = groups(g).outdir;
    a = mmil_parms2args(p);
    results = mmil_intrascan_stdev(fnamelist_mgz(gidx),a{:});
    % SynthStrip mask mmil_intrascan_stdev built from this group's inter-mean
    group_mask = sprintf('%s/func_inter_mean_mask.mgz',groups(g).outdir);
    % session-level QC jpgs only for the base group; per-scan jpgs always
    copy_stdev_jpgs(fnamelist_mgz(gidx),groups(g).outdir,parms,isempty(groups(g).suffix));
    scan_rms_vals(gidx) = calc_scan_diff_rms(fnamelist_mgz(gidx),groups(g).outdir,group_mask,parms);
    [scan_corr_vals(gidx),scan_com_dists(gidx)] = ...
      calc_scan_alignment(fnamelist_mgz(gidx),groups(g).outdir,group_mask);
    % session summaries: the least correlated and the farthest run of the group
    vals = scan_corr_vals(gidx);
    if any(~isnan(vals))
      stdev_qcinfo.(sprintf('inter_corr_min%s',groups(g).suffix)) = min(vals);
    end
    vals = scan_com_dists(gidx);
    if any(~isnan(vals))
      stdev_qcinfo.(sprintf('inter_com_dist_max%s',groups(g).suffix)) = max(vals);
    end
    if isempty(results), continue; end
    any_results = true;
    for t=1:length(tag_list)
      stdev_qcinfo.(sprintf('%s%s',tag_list{t},groups(g).suffix)) = ...
        get_result_val(results,tag_list{t});
    end
    % the inter-scan tSNR and its ratio to the intra-scan tSNR need 3 scans:
    %   with one the inter-scan stdev is 0 (Inf); with two it is a single
    %   difference per voxel, and the long tail of the voxelwise ratio sets
    %   the mean; n/a instead, for the medians too
    if length(gidx) < 3
      stdev_qcinfo.(sprintf('inter_tSNR%s',groups(g).suffix)) = NaN;
      stdev_qcinfo.(sprintf('inter_intra_tSNR_ratio%s',groups(g).suffix)) = NaN;
      stdev_qcinfo.(sprintf('inter_tSNR_median%s',groups(g).suffix)) = NaN;
      stdev_qcinfo.(sprintf('inter_intra_tSNR_median_ratio%s',groups(g).suffix)) = NaN;
    end
  end
  if parms.verbose, toc; end

  if ~any_results, return; end

  % return the base-group SynthStrip brain mask (used by b0dx RMS + meanNorm)
  tmp_mask = sprintf('%s/func_inter_mean_mask.mgz',stdev_outdir);
  if exist(tmp_mask,'file'), fname_brain = tmp_mask; end

  stdev_row = convert_func_row(stdev_qcinfo,parms.participant_id,parms.session_id);

  fname_mean = sprintf('%s/func_inter_mean.mgz',stdev_outdir);
  fname_mean_norm = sprintf('%s/%s_desc-meanNorm_boldref.nii.gz',parms.outdir,parms.outstem);
  if exist(fname_mean,'file') && (~exist(fname_mean_norm,'file') || parms.forceflag)
    if parms.verbose
      fprintf('%s: normalizing interscan mean for manual QC...\n',mfilename);
    end
    [vol_mean,M,mr_parms] = mmil_load_mgh(fname_mean);
    if ~isempty(fname_brain) && exist(fname_brain,'file')
      fname_mask = fname_brain;
    else
      fname_mask = sprintf('%s/func_inter_mean_mask.mgz',stdev_outdir);
    end
    [vol_mask] = fs_load_mgh(fname_mask);
    vals = vol_mean(vol_mask>0 & vol_mean>0);
    topval = prctile(vals,parms.fMRI_topval_prctile);
    vol_mean_norm = 100 * vol_mean / topval;
    vol_mean_norm(vol_mean_norm<0) = 0;
    vol_mean_norm(vol_mean_norm>100) = 100;
    fs_save_nifti(vol_mean_norm,fname_mean_norm,M,mr_parms);
  end
  if exist(fname_mean_norm,'file') && ~isempty(groups(1).idx)
    fnames_src = parms.fnames_func(groups(1).idx);
    info = [];
    info.Description = ['mean across runs of the processed BOLD time series means, '...
      'scaled to 100 at the top percentile within the brain mask'];
    info.Sources = cellfun(@(x) proc_uri(x,parms),fnames_src,'UniformOutput',false);
    info.SpatialReference = get_ref_space(fnames_src{1});
    write_json(info,fname_mean_norm,parms);
  end

  % resample normalized interscan mean to T1w space for manual QC
  fname_mean_norm_res = sprintf('%s/%s_space-T1w_desc-meanNorm_boldref.nii.gz',...
    parms.outdir,parms.outstem);
  if exist(fname_mean_norm,'file') && ~isempty(parms.fnames_t1w) &&...
    (~exist(fname_mean_norm_res,'file') || parms.forceflag)
    M_T1_to_EPI = get_M_T1_to_EPI(parms);
    if ~isempty(M_T1_to_EPI)
      if parms.verbose
        fprintf('%s: resampling normalized interscan mean to T1w space...\n',mfilename);
      end
      try
        [vol_mean_norm,M] = cmig_read_nifti(fname_mean_norm);
        [~,M_T1,volsz_T1] = cmig_read_nifti(parms.fnames_t1w{1},1);
        [vol_res,M_res] = mmil_resample_vol(vol_mean_norm,M,...
          'M_ref',M_T1,'nvox_ref',volsz_T1(1:3),...
          'M_reg',M_T1_to_EPI,'interpm',2,'bclamp',1);
        fs_save_nifti(vol_res,fname_mean_norm_res,M_res);
      catch me
        fprintf('%s: WARNING: failed to resample normalized interscan mean to T1w space: %s\n',...
          mfilename,me.message);
      end
    else
      fprintf('%s: WARNING: M_T1_to_EPI not found in proc JSON; normalized interscan mean not resampled to T1w space\n',...
        mfilename);
    end
  end
  if exist(fname_mean_norm_res,'file') && ~isempty(parms.fnames_t1w)
    info = [];
    info.Description = ['normalized mean across runs of the processed BOLD '...
      'time series means, resampled to the processed T1w image'];
    info.Sources = {proc_uri(fname_mean_norm,parms,1),proc_uri(parms.fnames_t1w{1},parms)};
    info.SpatialReference = proc_uri(parms.fnames_t1w{1},parms);
    write_json(info,fname_mean_norm_res,parms);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [session_row,errcode,errmsg] = calc_fMRI_session_qc(parms)
  session_row = []; errcode = 0; errmsg = [];
  n = length(parms.fnames_func);
  if n == 0, return; end

  if parms.verbose
    fprintf('%s: computing session-level func QC...\n',mfilename);
  end

  qcinfo = [];
  qcinfo.VisitID                = sprintf('%s_%s',parms.participant_id,parms.session_id);
  qcinfo.stype                  = 'func_session';
  qcinfo.nscans                 = n;
  qcinfo.fstem                  = sprintf('%s_func_session',parms.outstem);
  qcinfo.brain_cutoff           = NaN;
  qcinfo.brain_cutoff_dorsal    = NaN;
  qcinfo.brain_cutoff_ventral   = NaN;

  [qcinfo,errcode,errmsg] = cmig_detect_brain_cutoff(qcinfo,...
    parms.fnames_func,parms.fnames_info_func,'func',...
    'tmpdir',parms.tmpdir,'outdir',parms.outdir,'outstem',parms.outstem,...
    'fname_plot',sprintf('%s/%s_braincutoff.jpg',parms.outdir,parms.outstem),...
    'fspath',parms.fspath,'verbose',parms.verbose,'forceflag',parms.forceflag);
  if errcode, return; end

  session_row = convert_func_row(qcinfo,parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fwhm_proc,fwhm_dt] = calc_func_fwhm(parms,fname_brain)
  % FWHM per processed run, twice: undetrended (the processed image itself) and
  % detrended (its residual). Both use the SynthStrip mask mmil_intrascan_stdev
  % built from the interscan mean, so the two arms and all runs share one mask.
  % outstem keeps the two arms from colliding, since cmig_3dFWHM would
  % otherwise name both from the same input file stem
  n = length(parms.fnames_func);
  fwhm_proc = nan(n,3); fwhm_dt = nan(n,3);
  if isempty(fname_brain) || ~exist(fname_brain,'file')
    fprintf('%s: WARNING: no brain mask, skipping func FWHM\n',mfilename);
    return;
  end
  tmpdir = sprintf('%s/func_fwhm',parms.tmpdir);
  mmil_mkdir(tmpdir);
  arms = {'proc',0; 'dt',1};
  for f=1:n
    [~,fstem] = cmig_get_fileparts(parms.fnames_func{f});
    for a=1:size(arms,1)
      try
        fwhm = cmig_3dFWHM(parms.fnames_func{f},...
          'fname_mask',fname_brain,'outdir',tmpdir,...
          'outstem',sprintf('%s_%s',fstem,arms{a,1}),...
          'detrend_flag',arms{a,2},'verbose',0,'forceflag',parms.forceflag);
      catch me
        fprintf('%s: WARNING: FWHM (%s) failed for %s: %s\n',...
          mfilename,arms{a,1},fstem,me.message);
        continue;
      end
      if strcmp(arms{a,1},'proc')
        fwhm_proc(f,:) = fwhm;
      else
        fwhm_dt(f,:) = fwhm;
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function skipTRs = get_skipTRs(parms)
  if ~isempty(parms.skipTRs)
    skipTRs = parms.skipTRs;
    return;
  end
  skipTRs = 0;
  for f=1:length(parms.fnames_func)
    fname_info = parms.fnames_info_func{f};
    if isempty(fname_info) || ~exist(fname_info,'file')
      fname_info = regexprep(parms.fnames_func{f},'\.nii.*','.json');
    end
    if ~exist(fname_info,'file'), continue; end
    try
      scan_info = jsondecode(fileread(fname_info));
      p = mmil_parms2args(parms,{'skipTRs_GE_DV25','skipTRs_GE','skipTRs_Philips','skipTRs_Siemens'});
      skipTRs = cmig_set_abcd_skipTRs(scan_info,p{:});
    catch
    end
    break;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M_T1_to_EPI = get_M_T1_to_EPI(parms)
  M_T1_to_EPI = [];
  for f=1:length(parms.fnames_func)
    fname_info = parms.fnames_info_func{f};
    if isempty(fname_info) || ~exist(fname_info,'file')
      fname_info = regexprep(parms.fnames_func{f},'\.nii.*','.json');
    end
    if ~exist(fname_info,'file'), continue; end
    try
      json = jsondecode(fileread(fname_info));
      M_T1_to_EPI = cmig_get_reg_T1w(json);
    catch
    end
    if ~isempty(M_T1_to_EPI), break; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post = pe_toward_posterior(info)
  % true when the phase encoding runs toward posterior (AP): the
  %   PhaseEncodingDirection axis and sign through the image's vox2ras
  %   (Affine), so right in any output orientation (std_orient of
  %   cmig_proc_bids_func); without an Affine, or for an axis with no A/P
  %   component, the sign as in dcm2niix's LAS, where j points anterior
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

function revflags = get_func_revflags(parms)
  % reverse-encoded when the phase encoding runs toward posterior, j- in
  %   dcm2niix's LAS (pe_toward_posterior); missing/unreadable direction is
  %   treated as forward
  n = length(parms.fnames_func);
  revflags = false(1,n);
  for f=1:n
    fname_info = parms.fnames_info_func{f};
    if isempty(fname_info) || ~exist(fname_info,'file')
      fname_info = regexprep(parms.fnames_func{f},'\.nii.*','.json');
    end
    if ~exist(fname_info,'file'), continue; end
    try
      scan_info = jsondecode(fileread(fname_info));
      post = pe_toward_posterior(scan_info);
      if ~isempty(post)
        revflags(f) = post;
      end
    catch
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = get_result_val(results,tag)
  if ~isempty(results) && isfield(results,tag)
    v = results.(tag);
  else
    v = NaN;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = convert_func_row(qcinfo, participant_id, session_id)
  row.participant_id = participant_id;
  row.session_id     = session_id;
  switch qcinfo.stype
    case 'func'
      % series-identifying entity columns
      fstem = mmil_getfield(qcinfo,'fstem','');
      row.modality = 'bold';
      row.acq_time = qcinfo.StudyDate;
      fstem = rm_desc_preproc(fstem);
      row.filename = sprintf('func/%s.nii.gz', fstem);
      tok = regexp(fstem,'(?:^|_)task-([^_]+)','tokens','once');
      if ~isempty(tok), row.task = tok{1}; else row.task = 'n/a'; end
      tok = regexp(fstem,'(?:^|_)run-([^_]+)','tokens','once');
      if ~isempty(tok), row.run = tok{1}; else row.run = 'n/a'; end
      tok = regexp(fstem,'(?:^|_)acq-([^_]+)','tokens','once');
      if ~isempty(tok), row.acq = tok{1}; else row.acq = 'n/a'; end
      tok = regexp(fstem,'(?:^|_)rec-([^_]+)','tokens','once');
      if ~isempty(tok), row.rec = tok{1}; else row.rec = 'n/a'; end
    case 'func_stdev'
      row.mr_y_qc__post__aut__fmri__nscans_count              = qcinfo.nscans;
      row.mr_y_qc__post__aut__fmri__nscans__for_count          = qcinfo.nscans_for;
      row.mr_y_qc__post__aut__fmri__nscans__rev_count          = qcinfo.nscans_rev;
      row.mr_y_qc__post__aut__fmri__inter__signal_mean        = qcinfo.inter_mean;
      row.mr_y_qc__post__aut__fmri__inter__signal_std        = qcinfo.inter_stdev;
      row.mr_y_qc__post__aut__fmri__intra__signal_std        = qcinfo.intra_stdev;
      row.mr_y_qc__post__aut__fmri__inter__temp__snr_mean    = qcinfo.inter_tSNR;
      row.mr_y_qc__post__aut__fmri__intra__temp__snr_mean    = qcinfo.intra_tSNR;
      row.mr_y_qc__post__aut__fmri__inter__signal_cov        = qcinfo.inter_cov;
      row.mr_y_qc__post__aut__fmri__intra__signal_cov        = qcinfo.intra_cov;
      row.mr_y_qc__post__aut__fmri__inter__intra__std_ratio  = qcinfo.inter_intra_stdev_ratio;
      row.mr_y_qc__post__aut__fmri__inter__intra__temp__snr_ratio = qcinfo.inter_intra_tSNR_ratio;
      row.mr_y_qc__post__aut__fmri__inter__temp__snr_med     = qcinfo.inter_tSNR_median;
      row.mr_y_qc__post__aut__fmri__intra__temp__snr_med     = qcinfo.intra_tSNR_median;
      row.mr_y_qc__post__aut__fmri__inter__intra__temp__snr__med_ratio = qcinfo.inter_intra_tSNR_median_ratio;
      row.mr_y_qc__post__aut__fmri__inter__corr_min          = qcinfo.inter_corr_min;
      row.mr_y_qc__post__aut__fmri__inter__com__dist_max     = qcinfo.inter_com_dist_max;
      % reverse phase-encode direction, written only when the session has
      %   reverse func scans, so single-direction studies get no all-NaN columns
      if qcinfo.nscans_rev > 0
        row.mr_y_qc__post__aut__fmri__rev__inter__signal_mean      = qcinfo.inter_mean_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__signal_std       = qcinfo.inter_stdev_rev;
        row.mr_y_qc__post__aut__fmri__rev__intra__signal_std       = qcinfo.intra_stdev_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__temp__snr_mean   = qcinfo.inter_tSNR_rev;
        row.mr_y_qc__post__aut__fmri__rev__intra__temp__snr_mean   = qcinfo.intra_tSNR_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__signal_cov       = qcinfo.inter_cov_rev;
        row.mr_y_qc__post__aut__fmri__rev__intra__signal_cov       = qcinfo.intra_cov_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__intra__std_ratio  = qcinfo.inter_intra_stdev_ratio_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__intra__temp__snr_ratio = qcinfo.inter_intra_tSNR_ratio_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__temp__snr_med    = qcinfo.inter_tSNR_median_rev;
        row.mr_y_qc__post__aut__fmri__rev__intra__temp__snr_med    = qcinfo.intra_tSNR_median_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__intra__temp__snr__med_ratio = qcinfo.inter_intra_tSNR_median_ratio_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__corr_min     = qcinfo.inter_corr_min_rev;
        row.mr_y_qc__post__aut__fmri__rev__inter__com__dist_max = qcinfo.inter_com_dist_max_rev;
      end
    case 'func_session'
      row.mr_y_qc__post__aut__fmri__nscans_count               = qcinfo.nscans;
      row.mr_y_qc__post__aut__fmri__brain__cutoff_score        = qcinfo.brain_cutoff;
      row.mr_y_qc__post__aut__fmri__brain__cutoff__dorsal_score = qcinfo.brain_cutoff_dorsal;
      row.mr_y_qc__post__aut__fmri__brain__cutoff__ventral_score = qcinfo.brain_cutoff_ventral;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = add_proc_info(row,fname_info)
  % field map and between-scan start fields of the proc sidecar; n/a if absent
  info = [];
  try
    info = jsondecode(fileread(fname_info));
  catch
  end
  tags = {...
    'B0UnwarpFieldMapFallback','mr_y_qc__post__series__b0__fallback_indicator';...
    'B0UnwarpFieldMapShimMatch','mr_y_qc__post__series__b0__shim__match_indicator';...
    'B0UnwarpFieldMapDisplacement','mr_y_qc__post__series__b0__fmap_dist';...
    'B0UnwarpFieldMapCoverage','mr_y_qc__post__series__b0__fmap__cover_prop';...
    'BetweenScanStartDistance','mr_y_qc__post__series__inter__start_dist'};
  for t=1:size(tags,1)
    val = mmil_getfield(info,tags{t,1},NaN);
    if isempty(val) || ~(isnumeric(val) || islogical(val)), val = NaN; end
    row.(tags{t,2}) = double(val(1));
  end
  % coverage as text with 6 decimals in every row, as in the proc sidecar;
  %   the tsv writer prints an integer 1 as "1" and other values with 8 decimals
  fld = 'mr_y_qc__post__series__b0__fmap__cover_prop';
  if ~isnan(row.(fld)), row.(fld) = sprintf('%0.6f',row.(fld)); end
  val = mmil_getfield(info,'BetweenScanStart','n/a');
  if ~ischar(val) || isempty(val), val = 'n/a'; end
  row.mr_y_qc__post__series__inter__start_type = val;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_func_session_row(stdev_row,cutoff_row,fwhm_proc,fwhm_dt,participant_id,session_id,b0_unwarp_indicator,fmap_fallback_indicator)
  session_row = struct();
  session_row.participant_id = participant_id;
  session_row.session_id     = session_id;
  % NOTE: either part may be empty (e.g. all runs too short for the stdev
  %   step); its columns are then simply absent, as production's qcinfo did
  skip = {'participant_id','session_id'};
  rows = {stdev_row,cutoff_row};
  for r=1:length(rows)
    if isempty(rows{r}), continue; end
    flds = fieldnames(rows{r});
    for k=1:length(flds)
      if ~ismember(flds{k},skip)
        session_row.(flds{k}) = rows{r}.(flds{k});
      end
    end
  end
  session_row.mr_y_qc__post__aut__fmri__b0__unwarp_indicator = b0_unwarp_indicator;
  session_row.mr_y_qc__post__aut__fmri__b0__fallback_indicator = fmap_fallback_indicator;
  % session-level FWHM is the mean across runs of the per-run values, for both
  % arms; no per-run copies are written here
  axes_names = {'x','y','z'};
  for a=1:3
    session_row.(sprintf('mr_y_qc__post__aut__fmri__%sdir_fwhm',axes_names{a})) = ...
      mean(fwhm_proc(:,a),'omitnan');
    session_row.(sprintf('mr_y_qc__post__aut__fmri__dtr__%sdir_fwhm',axes_names{a})) = ...
      mean(fwhm_dt(:,a),'omitnan');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_stdev_jpgs(fnamelist_mgz,stdev_outdir,parms,sessflag)
  % sessflag: copy the session-level QC jpgs (default true); set false for
  % the reverse-direction group so its session jpgs do not overwrite the
  % base group's. Per-scan jpgs (keyed by unique fstem) are always copied.
  if ~exist('sessflag','var') || isempty(sessflag), sessflag = true; end
  session_map = {...
    'func_inter_mean',  'inter_mean';...
    'func_mask',        'inter_mask';...
    'func_inter_stdev', 'inter_std';...
    'func_intra_stdev', 'intra_std';...
    'func_inter_tSNR',  'inter_tsnr';...
    'func_intra_tSNR',  'intra_tsnr';...
  };
  if sessflag
    for j=1:size(session_map,1)
      src = sprintf('%s/%s.jpg',stdev_outdir,session_map{j,1});
      dst = sprintf('%s/%s_desc-%s.jpg',parms.outdir,parms.outstem,session_map{j,2});
      if exist(src,'file') && (~exist(dst,'file') || parms.forceflag)
        copyfile(src,dst);
      end
    end
  end
  scan_map = {'_mean','mean'; '_stdev','std'; '_tSNR','tsnr'};
  for f=1:length(fnamelist_mgz)
    [~,fstem] = fileparts(fnamelist_mgz{f});
    entities = regexprep(rm_desc_preproc(fstem),'_[^_]+$','');
    for j=1:size(scan_map,1)
      src = sprintf('%s/%s%s.jpg',stdev_outdir,fstem,scan_map{j,1});
      dst = sprintf('%s/%s_%s.jpg',parms.outdir,entities,scan_map{j,2});
      if exist(src,'file') && (~exist(dst,'file') || parms.forceflag)
        copyfile(src,dst);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_rms_vals = calc_scan_diff_rms(fnamelist_mgz,stdev_outdir,fname_brain,parms)
  % returns one value per element of fnamelist_mgz (group-local); the
  % caller places these at the corresponding global scan indices
  scan_rms_vals = nan(length(fnamelist_mgz),1);
  fname_inter_mean = sprintf('%s/func_inter_mean.mgz',stdev_outdir);
  if ~exist(fname_inter_mean,'file'), return; end
  vol_inter_mean = fs_load_mgh(fname_inter_mean);
  if ~isempty(fname_brain) && exist(fname_brain,'file')
    vol_mask = fs_load_mgh(fname_brain);
    mask = vol_mask > 0;
  else
    mask = true(size(vol_inter_mean));
  end
  for f=1:length(fnamelist_mgz)
    [~,fstem] = fileparts(fnamelist_mgz{f});
    fname_scan_mean = sprintf('%s/%s_mean.mgz',stdev_outdir,fstem);
    if ~exist(fname_scan_mean,'file'), continue; end
    vol_scan_mean = fs_load_mgh(fname_scan_mean);
    diff_vol = vol_scan_mean - vol_inter_mean;
    scan_rms_vals(f) = sqrt(mean(diff_vol(mask).^2));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scan_corr_vals,scan_com_dists] = calc_scan_alignment(fnamelist_mgz,stdev_outdir,fname_brain)
  % between-run alignment of each scan's temporal mean image, within the
  %   group brain mask, against the other scans of the group (leave one out):
  %   corr is the Pearson r with the voxelwise median of the others' means;
  %   com_dist is the distance (mm) from the intensity-weighted centre of mass
  %   (scanner coordinates, via the vox2ras) to the median of the others'
  %   centres of mass; NaN with fewer than 2 scans; with 2, both get the same
  % returns one value per element of fnamelist_mgz (group-local)
  nscans = length(fnamelist_mgz);
  scan_corr_vals = nan(nscans,1); scan_com_dists = nan(nscans,1);
  if nscans < 2, return; end
  fnames_mean = cell(nscans,1);
  for f=1:nscans
    [~,fstem] = fileparts(fnamelist_mgz{f});
    fnames_mean{f} = sprintf('%s/%s_mean.mgz',stdev_outdir,fstem);
  end
  mask = [];
  if ~isempty(fname_brain) && exist(fname_brain,'file')
    mask = fs_load_mgh(fname_brain) > 0;
  end
  [scan_corr_vals,scan_com_dists] = calc_image_alignment(fnames_mean,mask);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [corr_vals,com_dists] = calc_image_alignment(fnames,mask)
  % leave-one-out correlation and centre of mass distance of a set of
  %   images; a missing image is NaN and is left out of the others' references
  n = length(fnames);
  corr_vals = nan(n,1); com_dists = nan(n,1);
  vals = []; coms = nan(n,3); found = false(n,1); xyz = [];
  for f=1:n
    if ~exist(fnames{f},'file'), continue; end
    [vol,M] = fs_load_mgh(fnames{f});
    if isempty(mask), mask = true(size(vol)); end
    if isempty(xyz)
      ind = find(mask);
      [i,j,k] = ind2sub(size(mask),ind);
      xyz = (M(1:3,:) * [i j k ones(length(ind),1)]')';
      vals = nan(length(ind),n);
    end
    vals(:,f) = vol(ind);
    % negative intensities would pull the weighted centre outside the image
    w = max(vals(:,f),0);
    if sum(w) > 0
      coms(f,:) = sum(bsxfun(@times,xyz,w),1) / sum(w);
    end
    found(f) = true;
  end
  for f=1:n
    others = find(found & (1:n)' ~= f);
    if ~found(f) || isempty(others), continue; end
    ref = median(vals(:,others),2);
    x = vals(:,f) - mean(vals(:,f));
    y = ref - mean(ref);
    corr_vals(f) = sum(x.*y) / sqrt(sum(x.^2) * sum(y.^2));
    com_dists(f) = norm(coms(f,:) - median(coms(others,:),1));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = rm_desc_preproc(fstem)
  % the rawdata name of a proc file, which adds desc-preproc: the filename
  %   column joins with raw QC on it, and the jpg names carry the run entities
  fstem = regexprep(fstem,'_desc-preproc(?=_[^_]+$)','');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = proc_uri(fname,parms,postqc_flag)
  % BIDS URI of a file in the proc root (dataset cmig_proc),
  %   or with postqc_flag=1 in this postqc root (bids::)
  if ~exist('postqc_flag','var') || isempty(postqc_flag), postqc_flag = 0; end
  [~,fstem,fext] = fileparts(fname);
  if postqc_flag
    moddir = 'func'; dsname = '';
  else
    [~,moddir] = fileparts(fileparts(fname)); dsname = 'cmig_proc';
  end
  uri = sprintf('bids:%s:%s/%s/%s/%s%s',dsname,parms.participant_id,...
    parms.session_id,moddir,fstem,fext);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function space_ref = get_ref_space(fname_ref)
  % SpatialReference of a proc image, from its sidecar json; a bids:: URI
  %   there is local to the proc root, so restated as bids:cmig_proc:
  space_ref = 'orig';
  fname_json = regexprep(fname_ref,'\.nii(\.gz)?$','.json');
  if ~exist(fname_json,'file'), return; end
  try
    ref_info = jsondecode(fileread(fname_json));
  catch
    return;
  end
  if isfield(ref_info,'SpatialReference') && ~isempty(ref_info.SpatialReference)
    space_ref = regexprep(ref_info.SpatialReference,'^bids::','bids:cmig_proc:');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_json(info,fname_image,parms)
  % sidecar json for an output image
  fname_json = regexprep(fname_image,'\.nii(\.gz)?$','.json');
  if exist(fname_json,'file') && ~parms.forceflag, return; end
  % NOTE: PrettyPrint was added in R2021a (MATLAB 9.10)
  if verLessThan('matlab','9.10')
    json = jsonencode(info);
  else
    json = jsonencode(info,'PrettyPrint',true);
  end
  fid = fopen(fname_json,'wt');
  if fid<0
    fprintf('%s: WARNING: failed to open %s for writing\n',mfilename,fname_json);
    return;
  end
  fprintf(fid,'%s',json);
  fclose(fid);
return;
