function [outdir,errcode,errmsg] = cmig_analyze_bids_func(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_analyze_bids_func(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Analyze task and resting-state fMRI data from BIDS processed derivatives
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%   fsurf_rootdir: BIDS freesurfer derivatives root directory
%   output_rootdir: BIDS output root directory
%
% Optional Input ('key',value):
%   'task': task name (ABCD tasks only)
%     allowed: 'mid', 'nback', 'sst', 'rest'
%     {default = 'nback'}
%   'deconv_flag': [0|1] whether to run deconvolution
%     set to 0 automatically if no events files are found
%     {default = 1}
%   'stim_times_flag': [0|1|2|3] stimulus time file format for 3dDeconvolve
%     0 = .1D (old style)
%     1 = .txt (model each TR of trial as an event)
%     timing files from events_to_timing() in cmig_abcd_events_to_timing.py
%     2 = _block.txt and _dur.txt (model onset and duration)
%     3 = _block.txt and dur = 0 (model onset with impulse)
%     if empty, set based on task: nback=2, mid/sst=3
%     {default = []}
%   'block_gap': nback: maximum seconds between one event's offset and the
%     next event's onset for the two to belong to the same block
%     {default = 3}
%   'minfrac': fraction of a TR an event must fill to mark that TR
%     (stim_times_flag 0 and 1; SST uses its own minimum durations)
%     {default = 0.5}
%   'max_dur_std': stim_times_flag 2: maximum standard deviation (s) across
%     runs of a condition's block duration; the fit across runs takes their
%     mean, and a larger spread is an error
%     {default = 1}
%   'stim_times_model': name of stimulus HRF model for 3dDeconvolve
%     allowed: 'SPMG','SPMG1','SPMG2','GAM','TENT','TENTzero','CSPLIN','CSPLINzero','BLOCK'
%     overridden by task-specific stim_times_model_{task} parm if task is known
%     {default = 'SPMG1'}
%   'stim_times_model_mid': stim_times_model override for MID task
%     {default = 'SPMG1'}
%   'stim_times_model_sst': stim_times_model override for SST task
%     {default = 'SPMG2'}
%   'stim_times_model_nback': stim_times_model override for nback task
%     {default = 'SPMG2'}
%   'censor_thresh_deconv': framewise displacement threshold for deconvolution motion censoring (mm)
%     {default = 0.9}
%   'censor_thresh_corr': framewise displacement threshold for correlation analysis (mm)
%     {default = 0.2}
%   'contig_thresh_corr': minimum number of contiguous uncensored frames for
%     the correlation censor, censvec. was implicit, taken from the concat
%     function's own default
%     {default = 5}
%   'contig_thresh': minimum number of contiguous uncensored frames required
%     for the deconvolution; 0 disables the contiguity rule
%     {default = 0}
%   'min_numTRs': minimum number of TRs (before censoring) required to include
%     a run; a run with fewer frames is left out of every analysis
%     {default = 100}
%   'min_ntpoints': minimum number of TRs after censoring required to include
%     a run in correlation analysis
%     {default = 64}
%   'prefilt_flag': [0|1] apply bandpass filter to time series and motion data
%     before nuisance regression in correlation analysis
%     {default = 0}
%   'mc_resp_flag': [0|1] notch filter motion estimates for aliased respiration
%     signals; applied once, before all downstream use (deconvolution regressors
%     and censoring, concatenated time series, correlation censoring)
%     {default = 1}
%   'resp_low': respiration rate low frequency cut-off (respirations/minute)
%     {default = 18.6}
%   'resp_high': respiration rate high frequency cut-off (respirations/minute)
%     {default = 25.7}
%   'skipTRs': number of TRs to skip at the start of each run
%     if empty, set based on manufacturer and software version
%     {default = []}
%   'ntpoints_{task}': most TRs kept per run after skipTRs, for rest, mid,
%     nback and sst; a longer run loses TRs at its end, so runs have the same
%     length whatever the manufacturer's skipTRs. cmig_concat_bids_sessions
%     allocates the same lengths by default. Inf keeps every TR
%     {default = 375, 403, 362, 437}
%   'TR': scan repetition time (seconds)
%     if empty, read from BIDS JSON sidecar; required for deconvolution
%     {default = []}
%   'pydir': directory of the python scripts called for the run selection
%     and the stimulus timing files
%     if empty, $CMIG_PIPELINE_DIR/python; set it to test working copies or another
%     cmig_pipeline install
%     {default = []}
%   'ignore_delay_flag': [0|1] skip the run selection's eprime-scan delay
%     check, for a session-task with the REDCap eprime_mismatch_ok flag
%     {default = 0}
%   'exclude_noncompliant_flag': [0|1] the run selection treats runs with
%     compliant 0 in the rawdata _scans.tsv as invalid
%     {default = 0}
%   'overshoot_thresh': seconds a run's events may end past the end of its
%     scan before the run is excluded
%     {default = 5}
%   'max_delay': minutes allowed between a run's E-Prime start and its scan
%     {default = 12.5}
%   'short_flag': [0|1] use shortened condition/contrast names (BIDS-compliant)
%     {default = 1}
%   'norm_flag': [0|1] pre-normalize voxelwise data to percent signal change before extraction
%     normalizes each run NIfTI to mean=100 using 3dTstat + 3dcalc; applied before
%     extracting ico/ROI time series; affects both task and rest analyses
%     {default = 1}
%   'norm_thresh': threshold for masking low-signal voxels during pre-normalization
%     voxels with mean below threshold are set to 0 in normalized output
%     {default = 10}
%   'despike_flag': [0|1] run 3dDespike on each run before extraction
%     independent of norm_flag: if norm_flag=1, despiked images are normalized
%     (mean computed from despiked data); if norm_flag=0, despiked images are
%     used directly for extraction, without normalization
%     {default = 1}
%   'baseline_flag': [0|1] extract model-based polort-0 baseline for post-normalization
%     adds -bout to 3dDeconvolve (or -Rbeta to 3dREMLfit); extracts polort-0
%     coefficient for each run, averages across runs → baseline map saved in
%     the temporary directory
%     {default = 1}
%   'njobs': number of 3dDeconvolve sub-processes (-jobs), splitting the
%     voxels; need not match the CPUs a job has, a larger number sharing them
%     {default = 1}
%   'outliers_flag': [0|1] censor outlier conditions/contrasts based on sem rms
%     {default = 1}
%   'outliers_thresh': sem rms threshold for censoring outliers
%     {default = 5}
%   'prep_tseries_flag': [0|1] generate preprocessed tseries files
%     {default = 1}
%   'extract_ico_order': ico orders for the surface time series, for every
%     task; the ico5 series also feeds the vertexwise deconvolution
%     {default = [3,5]}
%   'corr_ico_order': ico orders for the vertexwise correlation matrices, and
%     for the desc-prep, taskpred, residual and trimmed ico time series; an
%     ico5 corrmat is 3.4 G and is not concatenated across sessions
%     {default = 3}
%   'trim_minutes': list of minute values for trimming (all tasks); e.g. [5]
%     {default = []}
%   'rest_trim_minutes': additional minute values applied only when task=rest; e.g. [10]
%     {default = []}
%   'trim_rand_flag': [0|1] randomly select trim frames; if 0, use first N
%     random frames span all of the visit's scans whatever the amount of
%     censoring, so visits with the same number of scans but different
%     magnitude or frequency of head motion are sampled alike; the first N
%     reach further into the visit the more frames are censored
%     {default = 1}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     will have form of {output_rootdir}/{participant_id}/{session_id}/func
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string
%
% Created:  04/21/2026 by Don Hagler
% Prev Mod: 05/14/2026 by Don Hagler
% Prev Mod: 05/15/2026 by Don Hagler
% Prev Mod: 05/26/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 06/01/2026 by Don Hagler
% Prev Mod: 06/02/2026 by Don Hagler
% Prev Mod: 06/04/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 06/07/2026 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/30/2026 by Don Hagler
% Prev Mod: 06/30/2026 by Don Hagler
% Prev Mod: 07/01/2026 by Don Hagler
% Prev Mod: 07/02/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/10/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Prev Mod: 08/28/2026 by Don Hagler
% Prev Mod: 09/03/2026 by Don Hagler
% Prev Mod: 09/07/2026 by Don Hagler
% Prev Mod: 09/08/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/14/2026 by Don Hagler
% Prev Mod: 09/15/2026 by Don Hagler
% Prev Mod: 09/16/2026 by Don Hagler
% Prev Mod: 09/17/2026 by Don Hagler
% Prev Mod: 09/18/2026 by Don Hagler
% Prev Mod: 09/23/2026 by Don Hagler
% Prev Mod: 09/25/2026 by Don Hagler
% Prev Mod: 09/28/2026 by Don Hagler
% Prev Mod: 10/03/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/08/2026 by Don Hagler
% Prev Mod: 10/13/2026 by Don Hagler
% Prev Mod: 10/16/2026 by Don Hagler
% Prev Mod: 10/19/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Prev Mod: 11/02/2026 by Don Hagler
% Prev Mod: 11/03/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/12/2026 by Don Hagler
% Prev Mod: 11/15/2026 by Don Hagler
% Prev Mod: 11/20/2026 by Don Hagler
% Prev Mod: 11/21/2026 by Don Hagler
% Last Mod: 11/22/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on cmig_analyze_bids_func_task
%   Created:  04/21/2026 by Don Hagler
%   Last Mod: 06/10/2026 by Don Hagler
% Based on cmig_analyze_bids_func_rest
%   Created:  09/04/2025 by Don Hagler
%   Last Mod: 06/10/2026 by Don Hagler
% Based on cmig_proc_bids_anat
%   Created:  08/04/2025 by Don Hagler
%   Last Mod: 08/14/2025 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,5), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

% check parameters
parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin);

% locate input data
[parms,errcode,errmsg] = check_input(parms);
if errcode
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate);
  return;
end

% check whether processing required
[parms,procflag] = check_output(parms);

% process data
if procflag
  % notch filter motion estimates for aliased respiration signals
  %   applied here only; downstream mc_resp_flag defaults stay false
  if parms.mc_resp_flag
    [errcode,tmp_errmsg,parms] = filter_motion_tsvs(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end
  end

  % remove the dummy TRs, and pre-normalize voxelwise data to percent signal
  %   change and/or despike; called unconditionally so the dummies come off
  %   even when norm_flag and despike_flag are both 0
  [errcode,tmp_errmsg,parms] = normalize_runs(parms);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end

  % calculate average time series for each ROI for each run
  [errcode,tmp_errmsg,parms] = extract_tseries(parms);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end

  if parms.deconv_flag
    % the deconvolution sees only the runs with an events file
    dparms = subset_deconv_runs(parms);

    % prepare TR, motion, stim, and GLT files for deconvolution
    [errcode,tmp_errmsg,dparms] = prep_deconv(dparms);
    if errcode && strncmp(tmp_errmsg,'not processed:',14)
      % a refusal goes first, where the rounds look for its prefix
      errmsg = cmig_append_errmsg(tmp_errmsg,errmsg);
    else
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    end
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end

    % run per-run design deconvolution to produce per-run design matrices
    [errcode,tmp_errmsg,dparms] = prep_design_per_run(dparms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end

    % the concat below indexes the designs by session run, and keeps every run,
    %   so an excluded run gets an empty entry. cmig_concat_fmri_roi_tseries
    %   loads a design only when its entry names a file, leaving that run's
    %   design zero, so it contributes nothing to taskpred or residual while
    %   the time series keeps its full length
    parms.fnames_design_per_run = cell(length(parms.fnames_data),1);
    parms.fnames_design_per_run(dparms.run_index) = dparms.fnames_design_per_run;
  end

  % concatenate ROI time series data across runs
  if parms.roi_flag
    % concatenate combined-atlas ROI time series (gpaseg, dskaseg, dstaseg, gpdstdskaseg)
    try
      concat_roi_tseries(parms);
    catch me
      [errcode,errmsg] = step_failed('concat_roi_tseries',me,errmsg,parms,startdate); return;
    end
    if parms.prep_tseries_flag
      try
        prep_roi_tseries(parms);
      catch me
        [errcode,errmsg] = step_failed('prep_roi_tseries',me,errmsg,parms,startdate); return;
      end
      if ~isempty(parms.trim_ntpoints_list)
        try
          trim_roi_tseries(parms);
        catch me
          [errcode,errmsg] = step_failed('trim_roi_tseries',me,errmsg,parms,startdate); return;
        end
      end
    end
  end

  % save vol_info for roi and ico
  if parms.roi_flag
    try
      cmig_save_vol_info(parms.participant_id,parms.session_id,...
        'outdir',parms.outdir,...
        'outstem',parms.outstem,...
        'forceflag',parms.forceflag);
    catch me
      [errcode,errmsg] = step_failed('cmig_save_vol_info',me,errmsg,parms,startdate); return;
    end
  end
  if parms.ico_flag
    % NOTE: a parms field added here, between the two deconv_flag blocks, is
    %   not seen by the second one unless refresh_deconv_parms carries it: the
    %   deconvolution works from a run-subset copy taken before this point.
    %   Omitting it fails at the read with Unrecognized field name
    for ico_order = parms.extract_ico_order
      try
        fname_vol_info = cmig_save_vol_info(parms.participant_id,parms.session_id,...
          'outdir',parms.outdir,...
          'outstem',parms.outstem,...
          'fname_tag',sprintf('atlas-ico%d',ico_order),'ico_order',ico_order,...
          'forceflag',parms.forceflag);
      catch me
        [errcode,errmsg] = step_failed('cmig_save_vol_info',me,errmsg,parms,startdate); return;
      end
      parms.(sprintf('fname_vol_info_ico%d',ico_order)) = fname_vol_info;
    end
  end

  % compute ROI-ROI correlation matrices and write tabulated data
  if parms.roi_flag
    try
      run_roi_corr(parms);
    catch me
      [errcode,errmsg] = step_failed('run_roi_corr',me,errmsg,parms,startdate); return;
    end
  end

  % concatenate ico time series data across runs and compute ico-ico correlations
  if parms.ico_flag
    for ico_order = parms.extract_ico_order
      try
        parms = concat_ico_tseries(parms,ico_order);
      catch me
        [errcode,errmsg] = step_failed('concat_ico_tseries',me,errmsg,parms,startdate); return;
      end
    end
    if parms.prep_tseries_flag
      try
        prep_ico_tseries(parms);
      catch me
        [errcode,errmsg] = step_failed('prep_ico_tseries',me,errmsg,parms,startdate); return;
      end
      if ~isempty(parms.trim_ntpoints_list)
        try
          trim_ico_tseries(parms);
        catch me
          [errcode,errmsg] = step_failed('trim_ico_tseries',me,errmsg,parms,startdate); return;
        end
      end
    end
    try
      calc_ico_corr(parms);
    catch me
      [errcode,errmsg] = step_failed('calc_ico_corr',me,errmsg,parms,startdate); return;
    end
  end

  if parms.deconv_flag
    % dparms still holds what the preparation above set up; it only needs the
    %   handful of fields the stages in between added to parms
    dparms = refresh_deconv_parms(dparms,parms);

    % run per-run deconvolution (always, for each run individually)
    [errcode,tmp_errmsg] = run_deconv_per_run(dparms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end

    % run concatenated deconvolution (only if multiple runs)
    [errcode,tmp_errmsg] = run_deconv_concat_runs(dparms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end

    % with one run there is no concatenated pass, so its outputs are the
    %   session's; copy them to the session stem
    [errcode,tmp_errmsg] = copy_single_run_outputs(dparms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate); return; end
  else
    dparms = [];
  end

  % write task motion tabdat, with the dof of the 3dDeconvolve fits
  if parms.roi_flag && ~strcmp(parms.task,'rest')
    errmsg = cmig_append_errmsg(errmsg,write_task_motion_tabdat(parms,dparms));
  end

end

% cleanup tmpdirs
if parms.cleanup_flag
  errmsg = cmig_append_errmsg(errmsg,cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir));
end

% write status file
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate,procflag);

outdir = parms.outdir;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = step_failed(step,me,errmsg,parms,startdate)
  % a step with no errcode output threw: record it and write the status
  errcode = 1;
  tmp_errmsg = sprintf('%s: %s',step,me.message);
  fprintf('%s: ERROR: %s\n',mfilename,tmp_errmsg);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,sprintf('task-%s_desc-analyze',parms.task),errcode,errmsg,parms.status_extra,startdate);
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
    'task','rest',[],...
    'short_flag',true,[false true],...
    'deconv_flag',true,[false true],...
    'stim_times_flag',[],[0:3],... % [] = set from task
    'block_gap',3,[0,Inf],... % s; nback: a longer gap starts a new block
    'minfrac',0.5,[0,1],... % fraction of a TR an event must fill to mark it (flags 0, 1)
    'max_dur_std',1,[0,Inf],... % s; flag 2: max std of a condition's block duration across runs
    'stim_times_model','SPMG1',{'SPMG','SPMG1','SPMG2','GAM','TENT','TENTzero','CSPLIN','CSPLINzero','BLOCK'},...
    'stim_times_model_mid','SPMG1',{'SPMG','SPMG1','SPMG2','GAM','TENT','TENTzero','CSPLIN','CSPLINzero','BLOCK'},...
    'stim_times_model_sst','SPMG2',{'SPMG','SPMG1','SPMG2','GAM','TENT','TENTzero','CSPLIN','CSPLINzero','BLOCK'},...
    'stim_times_model_nback','SPMG2',{'SPMG','SPMG1','SPMG2','GAM','TENT','TENTzero','CSPLIN','CSPLINzero','BLOCK'},...
    'censor_thresh_deconv',0.9,[0.001,Inf],...
    'censor_thresh_corr',0.2,[0.001,Inf],...
    'contig_thresh',0,[0,20],...
    'contig_thresh_corr',5,[0,20],...
    'min_numTRs',100,[1,Inf],... % min TRs before censoring to include a run
    'min_ntpoints',64,[1,Inf],... % min TRs after censoring to include run in corr
    'prefilt_flag',false,[false true],... % bandpass before nuisance regression in corr
    'mc_resp_flag',true,[false true],... % notch filter motion for aliased respiration
    'resp_low',18.6,[],... % respiration rate low cut-off (respirations/minute)
    'resp_high',25.7,[],... % respiration rate high cut-off (respirations/minute)
    'skipTRs',[],[0,1000],... % [] = set from manufacturer/software version
    'TR',[],[],... % [] = read from BIDS JSON sidecar
    'norm_flag',true,[false true],...
    'norm_thresh',10,[0,Inf],...
    'despike_flag',true,[false true],...
    'baseline_flag',true,[false true],...
    'outliers_flag',true,[false true],... % censor outlier conditions/contrasts based on sem rms
    'outliers_thresh',5,[0,Inf],... % sem rms threshold for censoring outliers
  ...
    'roi_flag',true,[false true],... % run ROI concat and deconvolution
    'ico_flag',true,[false true],... % run ico surface concat and deconvolution
    'vox_flag',true,[false true],... % run voxelwise deconvolution
    'prep_tseries_flag',true,[false true],... % generate preprocessed tseries files
    'trim_minutes',[],[],... % minute values for trimming (all tasks); e.g. [5]
    'rest_trim_minutes',[],[],... % additional minute values applied only when task=rest; e.g. [10]
    'trim_rand_flag',true,[false true],... % [0|1] randomly select trim frames; if 0, use first N
    'cleanup_flag',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'ntpoints_rest',375,[1,Inf],... % most TRs kept per run after skipTRs
    'ntpoints_mid',403,[1,Inf],...
    'ntpoints_nback',362,[1,Inf],...
    'ntpoints_sst',437,[1,Inf],...
    'skipTRs_GE_DV25',4,[],...  % skipTRs for GE DV25/Orchestra
    'skipTRs_GE',15,[],...      % skipTRs for GE other software
    'skipTRs_Philips',8,[],...  % skipTRs for Philips
    'skipTRs_Siemens',8,[],...  % skipTRs for Siemens
    'pydir',[],[],... % python scripts directory; [] = $CMIG_PIPELINE_DIR/python
    'ignore_delay_flag',false,[false true],... % REDCap eprime_mismatch_ok: skip the eprime-scan check
    'exclude_noncompliant_flag',false,[false true],... % compliant 0 runs invalid in the run selection
    'overshoot_thresh',5,[0,Inf],... % s the events may end past the scan end
    'max_delay',12.5,[0,Inf],... % min allowed between E-Prime run start and scan
    'col_run_tag','',[],... % per-run column tag, the release slot (r01, r02)
    'fnames_design_per_run',[],[],...
    'stim_labels',[],[],...
    'hemilist',{'lh','rh'},[],...
    'iresp_flag',false,[false true],... % output impulse response functions
    'ico_order',5,[0 7],...
    'extract_ico_order',[3,5],[],...
    'corr_ico_order',3,[0 7],...
    'extract_ico_smooth',[64,16],[],...
    'run_tag','',[],...  % prefix for per-run tmpdir output files
  ...
    'censor_flag',true,[false true],... % apply motion censoring
    'deriv_flag',true,[false true],... % include motion derivatives
    'errts_flag',false,[false true],... % output residual time series
    'reml_flag',false,[false true],... % use 3dREMLfit for pre-whitening
    'goforit_flag',true,[false true],... % use AFNI -GOFORIT
    'goforit_val',10,[],... % -GOFORIT threshold
    'njobs',1,[1,32],... % 3dDeconvolve -jobs
    'detrend',3,[0:3],... % polort order for 3dDeconvolve
    'minlag',0,[0,10],... % min lag in TRs (stim_times_flag=0 only)
    'maxlag',4,[0,30],... % max lag in TRs (stim_times_flag=0 only)
    'stim_times_nbasis',5,[1,100],... % basis functions for TENT/CSPLIN models
    'motion_labels',{'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'},[],...
    'nodata_nt',[],[],...
    'nodata_tr',[],[],...
    'deconv_tags',{...
      'stim_fnames','stim_labels','stim_times_flag','stim_times_model','stim_times_nbasis',...
      'fname_motion','fname_censor','fname_deriv','motion_labels',...
      'glt_fnames','glt_labels',...
      'censor_flag','deriv_flag','iresp_flag','errts_flag','reml_flag',...
      'goforit_flag','goforit_val','njobs','detrend','minlag','maxlag',...
      'nodata_nt','nodata_tr',...
      'baseline_flag',...
      'forceflag',...
    },[],...
  ...
    'modality','func',[],...
    'series_type','bold',[],...
    'status_extra',[],[],...
  });

  if ~isempty(parms.task) && ~ismember(parms.task,{'mid','sst','nback','rest'})
    error('unsupported task name: %s (must be mid, sst, nback, or rest)',parms.task);
  end
  parms.ntpoints = Inf;
  if ~isempty(parms.task)
    parms.ntpoints = parms.(sprintf('ntpoints_%s',parms.task));
  end

  % set stim_times_flag based on task if not specified by user
  if isempty(parms.stim_times_flag) && ~isempty(parms.task)
    switch lower(parms.task)
      case 'nback'
        parms.stim_times_flag = 2;
      case 'mid'
        parms.stim_times_flag = 3;
      case 'sst'
        parms.stim_times_flag = 3;
    end
  end


  % override stim_times_model with task-specific value if task is known
  if ~isempty(parms.task)
    task_model_field = sprintf('stim_times_model_%s',lower(parms.task));
    if isfield(parms,task_model_field)
      parms.stim_times_model = parms.(task_model_field);
    end
  end

  % set nreg_per_stim based on stim_times_model
  switch parms.stim_times_model
    case {'SPMG','SPMG2'}
      parms.nreg_per_stim = 2;
    case {'SPMG1','GAM','BLOCK'}
      parms.nreg_per_stim = 1;
    case {'TENT','TENTzero','CSPLIN','CSPLINzero'}
      parms.nreg_per_stim = parms.stim_times_nbasis;
  end

  % the condition order has to be settled before anything uses it: the time
  %   series extraction runs ahead of the deconvolution and would otherwise
  %   take the order the events file happens to list, which varies by session
  if ~isempty(parms.task)
    [parms.stim_labels,parms.conds_contrast,parms.condnames] = ...
      prep_contrasts(parms.task,parms.short_flag);
  end

  % the ico3 and ico5 time series are outputs for every task, rest included
  if isempty(parms.extract_ico_order)
    parms.extract_ico_order = [3,5];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];

  % set output directory
  parms.outdir = sprintf('%s/%s/%s/%s',parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);

  [parms,errcode,errmsg] = cmig_check_bids_func(parms,1);
  if errcode
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % warn if no events files found; disable deconvolution
  if all(cellfun(@isempty,parms.fnames_events))
    fprintf('%s: WARNING: no events files found; deconvolution will be skipped\n',mfilename);
    parms.deconv_flag = 0;
  end

  % the BIDS run number of each run, for the file names
  parms.run_numbers = bids_run_numbers(parms.fnames_data);

  % the runs to analyze, chosen by select_bids_func_runs, the selection the
  %   behavioral step uses too; every output uses them. A task with no events
  %   has no deconvolution and keeps all its runs
  if parms.deconv_flag
    [parms,errcode,errmsg] = select_runs(parms);
    if errcode
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
  end

  % a run too short to analyze (an aborted scan) is left out before the
  %   preprocessing, which some short runs fail (3dDespike needs 15 frames)
  [parms,errcode,errmsg] = drop_short_runs(parms);
  if errcode
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % a run with no events file has no design, so every stimulus column would be
  %   all zeros and 3dDeconvolve refuses the fit. Such a run is excluded from
  %   the deconvolution and kept for the data time series and correlations
  parms.ind_deconv_runs = find(~cellfun(@isempty,parms.fnames_events));
  if parms.deconv_flag && length(parms.ind_deconv_runs) < length(parms.fnames_data)
    fprintf('%s: WARNING: excluding %d of %d runs from deconvolution, no events file\n',...
      mfilename,length(parms.fnames_data)-length(parms.ind_deconv_runs),...
      length(parms.fnames_data));
  end

  [parms.fspath,parms.fs_version,parms.fs_build,parms.fs_complete_flag,errcode,errmsg] = ...
    cmig_check_bids_fsurf(parms.participant_id,parms.session_id,parms.fsurf_rootdir);
  parms.status_extra = struct('fs_version',parms.fs_version,...
    'fs_build',parms.fs_build,'fs_complete_flag',parms.fs_complete_flag);
  if errcode
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  % volume-only recon is fatal: surface sampling requires a complete recon
  if ~parms.fs_complete_flag
    errcode = 1;
    errmsg = sprintf('not processed: freesurfer recon %s is volume-only; surface sampling requires a complete recon',parms.fspath);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % set skipTRs and TR from first scan's info JSON if not provided
  if isempty(parms.skipTRs) || isempty(parms.TR)
    try
      scan_info = jsondecode(fileread(parms.fnames_info{1}));
    catch me
      errcode = 1;
      errmsg = sprintf('failed to load json file %s:\n%s',parms.fnames_info{1},me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    if isempty(parms.skipTRs)
      % fname_events lets it fall back to the dummy event when the sidecar
      %   carries no usable Manufacturer, as the release's mmps_mproc does not
      fname_events = [];
      if ~isempty(parms.fnames_events), fname_events = parms.fnames_events{1}; end
      parms.skipTRs = cmig_set_abcd_skipTRs(scan_info,...
        'skipTRs_GE_DV25',parms.skipTRs_GE_DV25,...
        'skipTRs_GE',parms.skipTRs_GE,...
        'skipTRs_Philips',parms.skipTRs_Philips,...
        'skipTRs_Siemens',parms.skipTRs_Siemens,...
        'fname_events',fname_events);
    end
    if isempty(parms.TR)
      parms.TR = scan_info.RepetitionTime;
    end
  end

  % compute combined trim minutes and ntpoints lists
  all_trim_minutes = parms.trim_minutes(:)';
  if strcmp(parms.task,'rest') && ~isempty(parms.rest_trim_minutes)
    all_trim_minutes = [all_trim_minutes, parms.rest_trim_minutes(:)'];
  end
  parms.all_trim_minutes = all_trim_minutes;
  if ~isempty(all_trim_minutes)
    parms.trim_ntpoints_list = round(all_trim_minutes * 60 / parms.TR);
  else
    parms.trim_ntpoints_list = [];
  end

  % seed for the random frame selection, fixed per visit but differing
  %   between visits, so visits with the same number of valid frames do not
  %   draw the same frames; the same for every atlas and trim length
  parms.trim_seed = visit_seed(parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seed = visit_seed(participant_id,session_id)
  % rolling hash of the visit ID string, below 2^31 (exact in double)
  str = sprintf('%s_%s',participant_id,session_id);
  seed = 0;
  for c = double(str)
    seed = mod(seed*31 + c,2^31-1);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dparms = refresh_deconv_parms(dparms,parms)
  % the deconvolution's parms is taken before the time series stages, which add
  %   these; none is run-length, so they carry over as they are. Copying
  %   this way rather than merging the deconvolution's own fields back into
  %   parms keeps the run-length fields of each on their own run list
  fields = {};
  for ico_order = parms.extract_ico_order
    fields{end+1} = sprintf('fname_vol_info_ico%d',ico_order);
    fields{end+1} = sprintf('fname_ico%d_concat',ico_order);
  end
  for f = 1:length(fields)
    fld = fields{f};
    if isfield(parms,fld), dparms.(fld) = parms.(fld); end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_numbers = deconv_run_numbers(parms,nruns)
  % the session run number of each run handed to the deconvolution; the two
  %   differ when a run was excluded for having no events file
  if isfield(parms,'run_numbers') && length(parms.run_numbers) == nruns
    run_numbers = parms.run_numbers;
  else
    run_numbers = 1:nruns;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_numbers = bids_run_numbers(fnames)
  % the run- entity of each file name; a name without one gets its position
  run_numbers = zeros(1,length(fnames));
  for r = 1:length(fnames)
    tok = regexp(fnames{r},'_run-(\d+)_','tokens','once');
    if ~isempty(tok)
      run_numbers(r) = str2double(tok{1});
    else
      run_numbers(r) = r;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_py,errcode,errmsg] = find_python_script(parms,fname)
  % a python script in pydir ($CMIG_PIPELINE_DIR/python if empty), and python3 to run it
  errcode = 0; errmsg = [];
  pydir = parms.pydir;
  if isempty(pydir), pydir = sprintf('%s/python',cmig_getenv('CMIG_PIPELINE_DIR')); end
  fname_py = sprintf('%s/%s',pydir,fname);
  if ~exist(fname_py,'file')
    errcode = 1;
    errmsg = sprintf('python script not found: %s',fname_py);
    return;
  end
  [status,result] = mmil_unix('which python3');
  if status || ~isempty(regexp(result,'not found','once'))
    errcode = 1;
    errmsg = sprintf('no python3 found to run %s',fname);
    return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = select_runs(parms)
  % keeps the selected runs, in slot order, in every per-run field. The
  %   selection is the behavioral step's (cmig_analyze_bids_behav, from
  %   rawdata), read from its {outstem}_desc-select_runs.tsv in outdir, when
  %   that step has succeeded there for this task; otherwise
  %   cmig_select_bids_func_runs.py runs on the proc tree, the events copies
  %   beside the bold, and writes that file. It is the record compiled into
  %   mr_y_qc__post__series. With events but no run selected, all runs are
  %   kept for the data and the deconvolution is off
  errcode = 0; errmsg = [];
  fname_sel = sprintf('%s/%s_%s_task-%s_desc-select_runs.tsv',...
    parms.outdir,parms.participant_id,parms.session_id,parms.task);
  fname_behav = sprintf('%s/%s_%s_task-%s_desc-behav_status.tsv',...
    parms.outdir,parms.participant_id,parms.session_id,parms.task);
  if exist(fname_sel,'file') && status_errcode(fname_behav) == 0
    fprintf('%s: run selection of the behavioral step: %s\n',mfilename,fname_sel);
  else
    [fname_py,errcode,errmsg] = find_python_script(parms,'cmig_select_bids_func_runs.py');
    if errcode, return; end
    cmd = sprintf(['python3 %s --participant_id %s --session_id %s --task %s',...
      ' --raw_rootdir %s --outdir %s --overshoot_thresh %g --max_delay %g',...
      ' --ignore_delay_flag %d --exclude_noncompliant_flag %d',...
      ' --forceflag 1 --verbose 1'],...
      fname_py,parms.participant_id,parms.session_id,parms.task,...
      parms.proc_rootdir,parms.outdir,parms.overshoot_thresh,parms.max_delay,...
      parms.ignore_delay_flag,parms.exclude_noncompliant_flag);
    [status,result] = mmil_unix(cmd);
    fprintf('%s',result);
    if status
      errcode = 1;
      errmsg = sprintf('run selection failed:\n%s',result);
      return;
    end
  end
  if ~exist(fname_sel,'file')
    errcode = 1;
    errmsg = sprintf('run selection wrote no %s',fname_sel);
    return;
  end
  % one row per bold run; the selected ones carry a slot
  lines = regexp(fileread(fname_sel),'\r?\n','split');
  lines = lines(~cellfun(@isempty,lines));
  % CollapseDelimiters false: an empty field (an unselected run's slot)
  %   must keep its column
  header = strsplit(lines{1},'\t','CollapseDelimiters',false);
  i_run  = find(strcmp(header,'run'));
  i_sel  = find(strcmp(header,'mr_y_qc__post__series__tfmri__select_indicator'));
  i_slot = find(strcmp(header,'mr_y_qc__post__series__tfmri__slot_num'));
  if isempty(i_run) || isempty(i_sel) || isempty(i_slot)
    errcode = 1;
    errmsg = sprintf('unexpected columns in %s',fname_sel);
    return;
  end
  sel_runs = []; sel_slots = [];
  for k = 2:length(lines)
    vals = strsplit(lines{k},'\t','CollapseDelimiters',false);
    if length(vals) < max([i_run,i_sel,i_slot]), continue; end
    if str2double(vals{i_sel}) == 1
      sel_runs(end+1) = str2double(vals{i_run});
      sel_slots(end+1) = str2double(vals{i_slot});
    end
  end
  [~,order] = sort(sel_slots);
  sel_runs = sel_runs(order);
  nruns = length(parms.fnames_data);
  % a run selected for the behavioral metrics only (not usable) is not among
  %   the session's runs
  ind = [];
  for k = 1:length(sel_runs)
    i = find(parms.run_numbers == sel_runs(k),1);
    if isempty(i)
      fprintf('%s: WARNING: selected run %d not among the session''s runs; left out\n',...
        mfilename,sel_runs(k));
      continue;
    end
    ind(end+1) = i;
  end
  if isempty(ind)
    fprintf('%s: WARNING: no run selected for task %s; data from all %d runs, no deconvolution\n',...
      mfilename,parms.task,nruns);
    parms.deconv_flag = 0;
    return;
  end
  if length(ind) < nruns
    fprintf('%s: using %d of %d runs (BIDS runs %s), as selected\n',...
      mfilename,length(ind),nruns,sprintf('%d ',parms.run_numbers(ind)));
  end
  run_fields = {'fnames_data','fnames_events','fnames_info','fnames_motion',...
                'scan_info','run_numbers'};
  for f = 1:length(run_fields)
    fld = run_fields{f};
    if isfield(parms,fld) && length(parms.(fld)) == nruns
      parms.(fld) = parms.(fld)(ind);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = drop_short_runs(parms)
  % leave out the runs with fewer than min_numTRs frames
  errcode = 0; errmsg = [];
  nruns = length(parms.fnames_data);
  nframes = zeros(1,nruns);
  for i = 1:nruns
    try
      [~,~,volsz] = cmig_read_nifti(parms.fnames_data{i},1);
    catch me
      errcode = 1;
      errmsg = sprintf('failed to read header of %s:\n%s',parms.fnames_data{i},me.message);
      return;
    end
    nframes(i) = volsz(4);
  end
  ind = find(nframes >= parms.min_numTRs);
  if length(ind) == nruns, return; end
  for i = setdiff(1:nruns,ind)
    fprintf('%s: WARNING: run %d has %d frames, fewer than min_numTRs (%d); left out\n',...
      mfilename,parms.run_numbers(i),nframes(i),parms.min_numTRs);
  end
  if isempty(ind)
    errcode = 1;
    errmsg = sprintf('not processed: every run has fewer than min_numTRs (%d) frames',parms.min_numTRs);
    return;
  end
  run_fields = {'fnames_data','fnames_events','fnames_info','fnames_motion',...
                'scan_info','run_numbers'};
  for f = 1:length(run_fields)
    fld = run_fields{f};
    if isfield(parms,fld) && length(parms.(fld)) == nruns
      parms.(fld) = parms.(fld)(ind);
    end
  end
  if parms.deconv_flag && all(cellfun(@isempty,parms.fnames_events))
    fprintf('%s: WARNING: no events files left; deconvolution will be skipped\n',mfilename);
    parms.deconv_flag = 0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errcode = status_errcode(fname)
  % the errcode of a status file (cmig_write_status), or NaN when the file is
  %   absent or has none
  errcode = NaN;
  if ~exist(fname,'file'), return; end
  lines = regexp(fileread(fname),'\r?\n','split');
  lines = lines(~cellfun(@isempty,lines));
  if length(lines) < 2, return; end
  header = strsplit(lines{1},'\t','CollapseDelimiters',false);
  vals = strsplit(lines{2},'\t','CollapseDelimiters',false);
  i = find(strcmp(header,'errcode'),1);
  if isempty(i) || i > length(vals), return; end
  errcode = str2double(vals{i});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dparms = subset_deconv_runs(parms)
  % a copy of parms holding only the runs that have an events file, so the
  %   deconvolution and everything derived from the design see those runs and
  %   the ~10 open-coded nruns = length(parms.fnames_data) sites follow
  dparms = parms;
  ind = parms.ind_deconv_runs;
  % the surviving runs keep their BIDS run numbers, so the per-run output stems
  %   still name the run the data came from; run_index maps per-run results
  %   back to the positions in parms
  if isfield(parms,'run_numbers') && length(parms.run_numbers) == length(parms.fnames_data)
    dparms.run_numbers = parms.run_numbers(ind);
  else
    dparms.run_numbers = ind(:)';
  end
  dparms.run_index = ind(:)';
  dparms.ind_deconv_runs = 1:length(ind);
  if length(ind) == length(parms.fnames_data), return; end
  run_fields = {'fnames_data','fnames_events','fnames_info','fnames_motion',...
                'fnames_norm','fnames_extract','scan_info'};
  for f = 1:length(run_fields)
    fld = run_fields{f};
    if isfield(dparms,fld) && length(dparms.(fld)) == length(parms.fnames_data)
      dparms.(fld) = dparms.(fld)(ind);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;

  % set output file stem and temporary directory
  if ~isempty(parms.task) && ~iscell(parms.task)
    parms.outstem = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  else
    parms.outstem = sprintf('%s_%s_task',parms.participant_id,parms.session_id);
  end
  if ~isempty(parms.task) && ~iscell(parms.task)
    parms.tmpdir = sprintf('%s/analyze_func_%s',parms.outdir,parms.task);
  else
    parms.tmpdir = sprintf('%s/analyze_func',parms.outdir);
  end

  % check whether output directory exists
  if ~exist(parms.outdir,'dir') || parms.forceflag
    procflag = 1;
  else
    % create list of expected output files
    output_files = {};

    % concat mat files
    % the ico entries follow the parameters that decide what is produced, or
    %   anything switched off is missing forever and every session reprocesses
    suffix_list = {'atlas-gpdstdskaseg_tseries',...
      'atlas-gpaseg_tseries','atlas-dskaseg_tseries','atlas-dstaseg_tseries'};
    for ico_order = parms.extract_ico_order
      suffix_list{end+1} = sprintf('atlas-ico%d_tseries',ico_order);
      suffix_list{end+1} = sprintf('atlas-ico%d_volinfo',ico_order);
    end
    for ico_order = parms.corr_ico_order
      suffix_list{end+1} = sprintf('atlas-ico%d_corrmat',ico_order);
    end
    if parms.prep_tseries_flag
      prep_atlases = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
      for ico_order = prep_ico_orders(parms)
        prep_atlases{end+1} = sprintf('ico%d',ico_order);
      end
      for pi = 1:length(prep_atlases)
        suffix_list{end+1} = sprintf('atlas-%s_desc-prep_tseries',prep_atlases{pi});
      end
      if parms.deconv_flag
        for pi = 1:length(prep_atlases)
          suffix_list{end+1} = sprintf('atlas-%s_desc-taskpred_tseries',prep_atlases{pi});
          suffix_list{end+1} = sprintf('atlas-%s_desc-residual_tseries',prep_atlases{pi});
        end
      end
      % trimmed tseries and corr
      if ~isempty(parms.trim_ntpoints_list)
        trim_roi_atlases = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
        trim_flags = [0];
        if parms.deconv_flag, trim_flags = [0,1,2]; end
        for ti = 1:length(parms.all_trim_minutes)
          trim_min = parms.all_trim_minutes(ti);
          for flag = trim_flags
            ttag = trim_desc_tag(flag,trim_min);
            for ra = 1:length(trim_roi_atlases)
              suffix_list{end+1} = sprintf('atlas-%s_%stseries',trim_roi_atlases{ra},ttag);
            end
            for ico_order = prep_ico_orders(parms)
              suffix_list{end+1} = sprintf('atlas-ico%d_%stseries',ico_order,ttag);
            end
          end
        end
      end
    end
    % trimmed ico corr (flag=0 only; consistent with non-trim ico corr in suffix_list above)
    if ~isempty(parms.trim_ntpoints_list)
      for ti = 1:length(parms.all_trim_minutes)
        trim_min = parms.all_trim_minutes(ti);
        % bids_tagged_instem inserts tag before last word, so
        %   '{outstem}_atlas-ico3_corrmat' + '_desc-trim5' → '{outstem}_atlas-ico3_desc-trim5_corrmat'
        ttag_s = trim_corr_tag(0,trim_min);
        ttag_s = ttag_s(2:end); % strip leading underscore
        for ico_order = parms.corr_ico_order
          suffix_list{end+1} = sprintf('atlas-ico%d_%s_corrmat',ico_order,ttag_s);
        end
      end
    end
    for s=1:length(suffix_list)
      fname_out = sprintf('%s/%s_%s.mat',parms.outdir,parms.outstem,suffix_list{s});
      output_files{end+1} = fname_out;
    end

    % deconv stats mat files (beta/sem for each condition, run, and hemisphere)
    condnames = {}; nruns = 0;
    if parms.deconv_flag && ~isempty(parms.task)
      condnames = parms.condnames;
      % only the runs with an events file are deconvolved, so only those have
      %   per-run stats; counting all of them would leave procflag always set
      deconv_runs = parms.ind_deconv_runs;
      nruns = length(deconv_runs);
      % the per-run files are named by BIDS run
      deconv_rnums = deconv_runs;
      if isfield(parms,'run_numbers') && length(parms.run_numbers) == length(parms.fnames_data)
        deconv_rnums = parms.run_numbers(deconv_runs);
      end
      % per-run stats (always produced), and the concat stats (only produced
      %   when nruns > 1)
      run_entities = {};
      for ri = 1:nruns
        run_entities{end+1} = sprintf('_run-%02d',deconv_rnums(ri));
      end
      if nruns > 1, run_entities{end+1} = ''; end
      for ri = 1:length(run_entities)
        run_entity = run_entities{ri};
        for h = 1:length(parms.hemilist)
          bids_hemi = cmig_get_bids_hemi(parms.hemilist{h});
          fstem_ico = sprintf('%s/%s%s_hemi-%s_atlas-ico5',parms.outdir,parms.outstem,run_entity,bids_hemi);
          for c = 1:length(condnames)
            output_files{end+1} = sprintf('%s_desc-%s_beta.mat',fstem_ico,condnames{c});
            output_files{end+1} = sprintf('%s_desc-%s_sem.mat',fstem_ico,condnames{c});
            if parms.ico_flag
              output_files{end+1} = sprintf('%s_desc-%s_beta.json',fstem_ico,condnames{c});
              output_files{end+1} = sprintf('%s_desc-%s_sem.json',fstem_ico,condnames{c});
            end
          end
          if parms.ico_flag
            output_files{end+1} = sprintf('%s_mask.mat',fstem_ico);
            output_files{end+1} = sprintf('%s_mask.func.gii',fstem_ico);
            output_files{end+1} = sprintf('%s_mask.json',fstem_ico);
          end
        end
        if parms.vox_flag
          for c = 1:length(condnames)
            output_files{end+1} = sprintf('%s/%s%s_desc-%s_beta.json',parms.outdir,parms.outstem,run_entity,condnames{c});
            output_files{end+1} = sprintf('%s/%s%s_desc-%s_sem.json',parms.outdir,parms.outstem,run_entity,condnames{c});
          end
        end
      end
    end

    % tabulated data TSV files
    if parms.roi_flag && ~isempty(parms.task)
      bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
      atlaslist = {'aseg','dsk','dst','gp'};
      % the fstat and deconv tabdat are written without gp
      deconv_atlaslist = {'aseg','dsk','dst'};
      % motion tabdat (one file per task)
      output_files{end+1} = sprintf('%s/%s_motionqc.tsv',parms.outdir,bids_base);
      % roi corrmat per combined atlas (flag=0 always; flags 1,2 when deconv)
      combined_atlases = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
      for ca_i = 1:length(combined_atlases)
        output_files{end+1} = sprintf('%s/%s_atlas-%s_corrmat.tsv',parms.outdir,bids_base,combined_atlases{ca_i});
      end
      for a = 1:length(atlaslist)
        output_files{end+1} = sprintf('%s/%s_atlas-%s_var.tsv', parms.outdir,bids_base,atlaslist{a});
      end
      if parms.deconv_flag
        for ca_i = 1:length(combined_atlases)
          output_files{end+1} = sprintf('%s/%s_atlas-%s_desc-taskpred_corrmat.tsv',parms.outdir,bids_base,combined_atlases{ca_i});
          output_files{end+1} = sprintf('%s/%s_atlas-%s_desc-residual_corrmat.tsv',parms.outdir,bids_base,combined_atlases{ca_i});
        end
        for a = 1:length(atlaslist)
          output_files{end+1} = sprintf('%s/%s_atlas-%s_desc-taskpred_var.tsv',  parms.outdir,bids_base,atlaslist{a});
          output_files{end+1} = sprintf('%s/%s_atlas-%s_desc-residual_var.tsv',  parms.outdir,bids_base,atlaslist{a});
        end
        % fstat and deconv tabdat per-run
        for ri = 1:nruns
          run_str = sprintf('run-%02d',deconv_rnums(ri));
          for a = 1:length(deconv_atlaslist)
            output_files{end+1} = sprintf('%s/%s_%s_atlas-%s_fstat.tsv',parms.outdir,bids_base,run_str,deconv_atlaslist{a});
            for c = 1:length(condnames)
              output_files{end+1} = sprintf('%s/%s_%s_atlas-%s_desc-%s_deconv.tsv',...
                parms.outdir,bids_base,run_str,deconv_atlaslist{a},condnames{c});
            end
          end
        end
        % session fstat and deconv tabdat (copy_single_run_outputs writes
        %   them for one run)
        if nruns > 0
          for a = 1:length(deconv_atlaslist)
            output_files{end+1} = sprintf('%s/%s_atlas-%s_fstat.tsv',parms.outdir,bids_base,deconv_atlaslist{a});
            for c = 1:length(condnames)
              output_files{end+1} = sprintf('%s/%s_atlas-%s_desc-%s_deconv.tsv',...
                parms.outdir,bids_base,deconv_atlaslist{a},condnames{c});
            end
          end
        end
      end
      % trimmed corrmat TSVs
      % NOTE: trimmed data is not tabulated, rest and task alike: it stays in
      %   the concatenated data, as in release 7.0
      if ~isempty(parms.trim_ntpoints_list)
        trim_flags = [0];
        if parms.deconv_flag, trim_flags = [0,1,2]; end
        for ti = 1:length(parms.all_trim_minutes)
          trim_min = parms.all_trim_minutes(ti);
          for flag = trim_flags
            tdesc = trim_tab_desc(flag,trim_min);
            for ca_i = 1:length(combined_atlases)
              output_files{end+1} = sprintf('%s/%s_atlas-%s_desc-%s_corrmat.tsv',parms.outdir,bids_base,combined_atlases{ca_i},tdesc);
            end
          end
        end
      end
    end

    % check whether output files exist
    for i=1:length(output_files)
      % cmig_write_fmri_corrmat writes no _corrmat.tsv for an all-NaN
      %   correlation
      if ~isempty(regexp(output_files{i},'_corrmat\.tsv$','once')) && ...
          ~corrmat_tsv_expected(output_files{i})
        continue;
      end
      if ~exist(output_files{i},'file')
        [~,name,ext] = fileparts(output_files{i});
        missing_file = [name,ext];
        fprintf('%s: missing output file %s\n',mfilename,missing_file);
        procflag = 1;
        break;
      end
    end
  end

  % create temporary directory (and output directory)
  if procflag
    clear_tmpdir_on_new_selection(parms);
    mmil_mkdir(parms.tmpdir);
    write_tmpdir_selection(parms);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flag = corrmat_tsv_expected(fname_tsv)
  % false when the matching _corrmat.mat exists with an all-NaN corrmat
  flag = true;
  fname_mat = regexprep(fname_tsv,'\.tsv$','.mat');
  if ~exist(fname_mat,'file'), return; end
  try
    tmp = load(fname_mat,'corrmat');
  catch
    return;
  end
  if ~isfield(tmp,'corrmat'), return; end
  if size(tmp.corrmat,1) < 1 || all(isnan(tmp.corrmat(1,:)))
    flag = false;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = selection_record(parms)
  % the runs whose files fill tmpdir, by BIDS run, and which have events
  str = sprintf('runs %s; events %s; deconv %d; ntpoints %g',...
    strtrim(sprintf('%d ',parms.run_numbers)),...
    strtrim(sprintf('%d ',parms.run_numbers(parms.ind_deconv_runs))),...
    parms.deconv_flag,parms.ntpoints);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clear_tmpdir_on_new_selection(parms)
  % tmpdir holds files named by position in the run list (stim, motion, GLT
  %   and deconvolution files), which a rerun skips when they exist; after a
  %   change of selection they belong to other runs, so tmpdir is removed.
  %   With the same selection it is kept, for a restart on partial files
  fname = sprintf('%s/selected_runs.txt',parms.tmpdir);
  if ~exist(fname,'file'), return; end
  prev = strtrim(fileread(fname));
  if strcmp(prev,selection_record(parms)), return; end
  fprintf('%s: run selection changed (%s, now %s); removing %s\n',...
    mfilename,prev,selection_record(parms),parms.tmpdir);
  rmdir(parms.tmpdir,'s');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tmpdir_selection(parms)
  fid = fopen(sprintf('%s/selected_runs.txt',parms.tmpdir),'wt');
  if fid < 0, return; end
  fprintf(fid,'%s\n',selection_record(parms));
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,parms] = extract_tseries(parms)
  errcode = 0; errmsg = [];
  try
    cmig_extract_bids_func_tseries(parms.fnames_extract,parms.fnames_info,parms.fspath,parms.tmpdir,...
      'fs_version',parms.fs_version,...
      'extract_ico_order',parms.extract_ico_order,...
      'extract_ico_smooth',parms.extract_ico_smooth,...
      'forceflag',parms.forceflag);
  catch me
    errcode = 1;
    errmsg = me.message;
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,parms] = normalize_runs(parms)
  errcode = 0; errmsg = [];
  nruns = length(parms.fnames_data);
  parms.fnames_norm = cell(nruns,1);
  for r = 1:nruns
    fname_in = parms.fnames_data{r};
    [~,fstem] = cmig_get_fileparts(fname_in);
    fname_skip      = sprintf('%s/%s_skip.nii.gz',parms.tmpdir,fstem);
    fname_despike   = sprintf('%s/%s_despike.nii.gz',parms.tmpdir,fstem);
    fname_mean      = sprintf('%s/%s_norm_mean.nii.gz',parms.tmpdir,fstem);
    fname_norm      = sprintf('%s/%s_norm.nii.gz',parms.tmpdir,fstem);
    % the dummy TRs, and any TRs past ntpoints, come off as an AFNI
    %   sub-brick selector on the first program in the chain, so no path
    %   writes a file just to drop them; sel is cleared once a program has
    %   applied it
    [sel,errcode,errmsg] = run_selector(fname_in,parms.skipTRs,parms.ntpoints);
    if errcode
      errmsg = sprintf('run %d: %s',r,errmsg);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    if parms.despike_flag
      if ~exist(fname_despike,'file') || parms.forceflag
        fprintf('%s: running 3dDespike for run %d...\n',mfilename,r);
        cmd = sprintf('3dDespike -prefix %s ''%s%s''',...
          fname_despike,fname_in,sel);
        [status,result] = mmil_unix(cmd);
        if status
          errcode = 1;
          errmsg = sprintf('3dDespike failed for run %d: %s',r,result);
          return;
        end
      end
      fname_for_norm = fname_despike;
      sel = '';
    else
      fname_for_norm = fname_in;
    end
    if ~parms.norm_flag
      if isempty(sel)
        parms.fnames_norm{r} = fname_for_norm;
      else
        % nothing upstream has written a file, so the cut needs one of its own
        if ~exist(fname_skip,'file') || parms.forceflag
          fprintf('%s: selecting TRs %s from run %d...\n',...
            mfilename,sel,r);
          cmd = sprintf('3dTcat -prefix %s ''%s%s''',...
            fname_skip,fname_for_norm,sel);
          [status,result] = mmil_unix(cmd);
          if status
            errcode = 1;
            errmsg = sprintf('3dTcat failed for run %d: %s',r,result);
            return;
          end
        end
        parms.fnames_norm{r} = fname_skip;
      end
      continue;
    end
    if ~exist(fname_mean,'file') || parms.forceflag
      fprintf('%s: computing mean for run %d normalization...\n',mfilename,r);
      cmd = sprintf('3dTstat -mean -prefix %s ''%s%s''',...
        fname_mean,fname_for_norm,sel);
      [status,result] = mmil_unix(cmd);
      if status
        errcode = 1;
        errmsg = sprintf('3dTstat failed for run %d: %s',r,result);
        return;
      end
    end
    if ~exist(fname_norm,'file') || parms.forceflag
      fprintf('%s: normalizing run %d to percent signal change...\n',mfilename,r);
      cmd = sprintf('3dcalc -a ''%s%s'' -b %s -expr ''100*a/b*ispositive(b-%g)'' -prefix %s',...
        fname_for_norm,sel,fname_mean,parms.norm_thresh,fname_norm);
      [status,result] = mmil_unix(cmd);
      if status
        errcode = 1;
        errmsg = sprintf('3dcalc normalization failed for run %d: %s',r,result);
        return;
      end
    end
    parms.fnames_norm{r} = fname_norm;
  end
  parms.fnames_extract = parms.fnames_norm;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sel,errcode,errmsg] = run_selector(fname,skipTRs,ntpoints)
  % AFNI sub-brick selector keeping TRs skipTRs+1 to skipTRs+ntpoints (as
  %   many as the run has), or '' for all of them
  errcode = 0; errmsg = []; sel = '';
  % the frame count from the NIfTI-1 header, dim[0..7] as int16 at byte 40;
  %   gzip -f passes an uncompressed file through; its stderr is dropped
  %   (in sh, as mmil_unix runs csh), since od closing the pipe makes it
  %   report a broken pipe
  [status,result] = mmil_unix(sprintf(...
    'sh -c "gzip -dcf %s 2>/dev/null | od -An -t d2 -j 40 -N 16"',fname));
  dims = sscanf(result,'%d');
  if status || numel(dims)~=8 || dims(1)<1 || dims(1)>7
    errcode = 1;
    errmsg = sprintf('failed to read the NIfTI dims of %s: %s',fname,result);
    return;
  end
  nv = 1;
  if dims(1)>=4, nv = dims(5); end
  last = min(nv,skipTRs+ntpoints) - 1;
  if skipTRs > 0 || last < nv-1
    sel = sprintf('[%d..%d]',skipTRs,last);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = fname_only(fpath)
  [~,name,ext] = fileparts(fpath);
  fname = [name ext];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function concat_roi_tseries(parms)
  nmotion = length(parms.motion_labels) * (1 + parms.deriv_flag);
  combined_list = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
  roitypes_list = {{'gp','aseg'}, {'dsk','aseg'}, {'dst','aseg'}, []};
  for i = 1:length(combined_list)
    atlas = combined_list{i};
    fname_out = sprintf('%s/%s_atlas-%s_tseries.mat',parms.outdir,parms.outstem,atlas);
    cmig_concat_fmri_roi_tseries(parms.fnames_extract,parms.fnames_motion,...
      parms.fnames_events,parms.tmpdir,...
      'fname_out',fname_out,'skipTRs',parms.skipTRs,...
      'TR',parms.TR,'fnames_design',parms.fnames_design_per_run,'nmotion',nmotion,...
      'stim_labels',parms.stim_labels,...
      'nreg_per_stim',parms.nreg_per_stim,'detrend',parms.detrend,...
      'censor_thresh',parms.censor_thresh_corr,...
      'contig_thresh',parms.contig_thresh_corr,...
      'censor_thresh_deconv',parms.censor_thresh_deconv,...
      'contig_thresh_deconv',parms.contig_thresh,...
      'roitypes',roitypes_list{i},'forceflag',parms.forceflag);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_roi_corr(parms,task_residual_flag,trim_min)
  if nargin < 3, trim_min = 0; end
  combined_list = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
  if trim_min > 0
    tag = trim_corr_tag(task_residual_flag,trim_min);
  else
    tag = corr_tag(task_residual_flag);
  end
  for i = 1:length(combined_list)
    atlas = combined_list{i};
    if parms.prep_tseries_flag
      if trim_min > 0
        ptag = trim_desc_tag(task_residual_flag,trim_min);
      else
        ptag = prep_desc_tag(task_residual_flag);
      end
      instem = sprintf('%s_atlas-%s_%stseries',parms.outstem,atlas,ptag);
      actual_task_residual_flag = 0;
      actual_preproc_flag = 0;
    else
      instem = sprintf('%s_atlas-%s_tseries',parms.outstem,atlas);
      actual_task_residual_flag = task_residual_flag;
      actual_preproc_flag = 1;
    end
    p = [];
    p.indir = parms.outdir;
    p.outdir = parms.outdir;
    p.instem = instem;
    p.infostem = sprintf('%s_volinfo',parms.outstem);
    p.outstem = bids_tagged_instem(sprintf('%s_atlas-%s_corrmat',parms.outstem,atlas),tag);
    p.task_residual_flag = actual_task_residual_flag;
    p.preproc_flag = actual_preproc_flag;
    if ~parms.prep_tseries_flag && trim_min > 0
      p.trim_ntpoints = round(trim_min * 60 / parms.TR);
      p.trim_rand_flag = parms.trim_rand_flag;
      p.trim_seed = parms.trim_seed;
    end
    p.censor_thresh = parms.censor_thresh_corr;
    p.min_numTRs = parms.min_numTRs;
    p.min_ntpoints = parms.min_ntpoints;
    p.prefilt_flag = parms.prefilt_flag;
    p.forceflag = parms.forceflag;
    a = mmil_parms2args(p);
    cmig_calc_fmri_roi_corr(a{:});
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = write_deconv_tabdat(parms)
  errmsg = [];
  % load condnames saved during prep_deconv
  fname_condnames = sprintf('%s/deconv_condnames.mat',parms.tmpdir);
  if ~exist(fname_condnames,'file')
    errmsg = 'condnames file not found, skipping deconv tabdat';
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end
  tmp = load(fname_condnames,'condnames');
  condnames = tmp.condnames;
  if isempty(condnames), return; end

  % deconv roi mat and run_tag for column names
  [~, col_run_tag, fname_deconv_roi] = deconv_roi_fname(parms);
  if ~exist(fname_deconv_roi,'file')
    errmsg = sprintf('deconv roi mat not found: %s',fname_deconv_roi);
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end

  [~,roi_instem] = cmig_get_fileparts(fname_deconv_roi);
  atlaslist = {'aseg','dsk','dst'};
  for c = 1:length(condnames)
    for a = 1:length(atlaslist)
      cmig_write_tabdat_tfmri('indir',parms.tmpdir,...
        'outdir',parms.outdir,...
        'participant_id',parms.participant_id,...
        'session_id',parms.session_id,...
        'colstem','mr_y_tfmri',...
        'atlas',atlaslist{a},...
        'fs_version',parms.fs_version,...
        'instem',roi_instem,...
        'condidx',c,...
        'condname',condnames{c},...
        'task',parms.task,...
        'run_tag',col_run_tag,...
        'fname_run_tag',parms.run_tag,...
        'forceflag',parms.forceflag);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = concat_ico_tseries(parms,ico_order)
  % NOTE: this runs between the two deconv_flag blocks, so a field set here
  %   reaches the second one only if refresh_deconv_parms carries it
  fname_field = sprintf('fname_ico%d_concat',ico_order);
  parms.(fname_field) = sprintf('%s/%s_atlas-ico%d_tseries.mat',parms.outdir,parms.outstem,ico_order);
  nmotion = length(parms.motion_labels) * (1 + parms.deriv_flag);
  cmig_concat_fmri_ico_tseries(parms.fnames_extract,parms.fnames_motion,...
    parms.fnames_events,parms.tmpdir,...
    'fname_out',parms.(fname_field),'skipTRs',parms.skipTRs,...
    'ico_order',ico_order,...
    'TR',parms.TR,'fnames_design',parms.fnames_design_per_run,'nmotion',nmotion,...
    'stim_labels',parms.stim_labels,...
    'nreg_per_stim',parms.nreg_per_stim,'detrend',parms.detrend,...
    'censor_thresh',parms.censor_thresh_corr,...
    'contig_thresh',parms.contig_thresh_corr,...
    'censor_thresh_deconv',parms.censor_thresh_deconv,...
    'contig_thresh_deconv',parms.contig_thresh,...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_ico_corr(parms)
  if parms.prep_tseries_flag && parms.deconv_flag
    flag_list = [0,1,2];
  else
    flag_list = [0];
  end
  for ico_order = parms.corr_ico_order
    fname_field = sprintf('fname_ico%d_concat',ico_order);
    if ~isfield(parms,fname_field) || isempty(parms.(fname_field)), continue; end
    for task_residual_flag = flag_list
      if parms.prep_tseries_flag
        ptag = prep_desc_tag(task_residual_flag);
        indir = parms.outdir;
        instem = sprintf('%s_atlas-ico%d_%stseries',parms.outstem,ico_order,ptag);
        actual_preproc_flag = 0;
      else
        [indir,instem] = fileparts(parms.(fname_field));
        actual_preproc_flag = 1;
      end
      tag = corr_tag(task_residual_flag);
      p = [];
      p.indir = indir;
      p.instem = instem;
      p.infostem = sprintf('%s_volinfo',parms.outstem);
      p.outdir = parms.outdir;
      p.outstem = bids_tagged_instem(sprintf('%s_atlas-ico%d_corrmat',parms.outstem,ico_order),tag);
      p.preproc_flag = actual_preproc_flag;
      if actual_preproc_flag
        p.roi_instem = sprintf('%s_atlas-gpdstdskaseg_tseries',parms.outstem);
      end
      p.censor_thresh = parms.censor_thresh_corr;
      p.min_numTRs = parms.min_numTRs;
      p.min_ntpoints = parms.min_ntpoints;
      p.prefilt_flag = parms.prefilt_flag;
      p.forceflag = parms.forceflag;
      a = mmil_parms2args(p);
      cmig_calc_fmri_ico_corr(a{:});
    end
    % trimmed ico corr
    if ~isempty(parms.trim_ntpoints_list)
      trim_flags = [0];
      if parms.deconv_flag, trim_flags = [0,1,2]; end
      for task_residual_flag = trim_flags
        for ti = 1:length(parms.all_trim_minutes)
          trim_min = parms.all_trim_minutes(ti);
          ttag = trim_corr_tag(task_residual_flag,trim_min);
          if parms.prep_tseries_flag
            ptag = trim_desc_tag(task_residual_flag,trim_min);
            p2_indir = parms.outdir;
            p2_instem = sprintf('%s_atlas-ico%d_%stseries',parms.outstem,ico_order,ptag);
            p2_preproc_flag = 0;
          else
            [p2_indir,~] = fileparts(parms.(fname_field));
            p2_instem = sprintf('%s_atlas-ico%d_tseries',parms.outstem,ico_order);
            p2_preproc_flag = 1;
          end
          p2 = [];
          p2.indir = p2_indir;
          p2.instem = p2_instem;
          p2.infostem = sprintf('%s_volinfo',parms.outstem);
          p2.outdir = parms.outdir;
          p2.outstem = bids_tagged_instem(sprintf('%s_atlas-ico%d_corrmat',parms.outstem,ico_order),ttag);
          p2.preproc_flag = p2_preproc_flag;
          if p2_preproc_flag
            p2.roi_instem = sprintf('%s_atlas-gpdstdskaseg_tseries',parms.outstem);
            p2.trim_ntpoints = round(trim_min * 60 / parms.TR);
            p2.trim_rand_flag = parms.trim_rand_flag;
            p2.trim_seed = parms.trim_seed;
          end
          p2.censor_thresh = parms.censor_thresh_corr;
          p2.min_numTRs = parms.min_numTRs;
          p2.min_ntpoints = parms.min_ntpoints;
          p2.prefilt_flag = parms.prefilt_flag;
          p2.forceflag = parms.forceflag;
          a2 = mmil_parms2args(p2);
          cmig_calc_fmri_ico_corr(a2{:});
        end
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prep_roi_tseries(parms)
  combined_list = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
  flag_list = [0];
  if parms.deconv_flag, flag_list = [0,1,2]; end
  for i = 1:length(combined_list)
    atlas = combined_list{i};
    instem = sprintf('%s_atlas-%s_tseries',parms.outstem,atlas);
    for task_residual_flag = flag_list
      ptag = prep_desc_tag(task_residual_flag);
      outstem = sprintf('%s_atlas-%s_%stseries',parms.outstem,atlas,ptag);
      cmig_prep_fmri_roi_tseries(...
        'indir',parms.outdir,'instem',instem,...
        'outdir',parms.outdir,'outstem',outstem,...
        'task_residual_flag',task_residual_flag,...
        'prefilt_flag',parms.prefilt_flag,...
        'TR',parms.TR,...
        'forceflag',parms.forceflag);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prep_ico_tseries(parms)
  flag_list = [0];
  if parms.deconv_flag, flag_list = [0,1,2]; end
  for ico_order = prep_ico_orders(parms)
    instem = sprintf('%s_atlas-ico%d_tseries',parms.outstem,ico_order);
    roi_instem = sprintf('%s_atlas-gpdstdskaseg_tseries',parms.outstem);
    for task_residual_flag = flag_list
      ptag = prep_desc_tag(task_residual_flag);
      outstem = sprintf('%s_atlas-ico%d_%stseries',parms.outstem,ico_order,ptag);
      cmig_prep_fmri_ico_tseries(...
        'indir',parms.outdir,'instem',instem,...
        'roi_instem',roi_instem,...
        'outdir',parms.outdir,'outstem',outstem,...
        'task_residual_flag',task_residual_flag,...
        'prefilt_flag',parms.prefilt_flag,...
        'TR',parms.TR,...
        'forceflag',parms.forceflag);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ico_orders = prep_ico_orders(parms)
  % the desc-prep, taskpred, residual and trimmed ico time series are made
  %   only for the orders with correlations, their one reader
  ico_orders = intersect(parms.extract_ico_order,parms.corr_ico_order);
  ico_orders = ico_orders(:)';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tag = prep_desc_tag(task_residual_flag)
  switch task_residual_flag
    case 1,  tag = 'desc-taskpred_';
    case 2,  tag = 'desc-residual_';
    otherwise, tag = 'desc-prep_';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lbl = trim_label(minutes)
  lbl = sprintf('trim%d',minutes);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tag = trim_corr_tag(task_residual_flag,minutes)
  lbl = trim_label(minutes);
  switch task_residual_flag
    case 1,  tag = sprintf('_desc-taskpred%s',lbl);
    case 2,  tag = sprintf('_desc-residual%s',lbl);
    otherwise, tag = sprintf('_desc-%s',lbl);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tag = trim_desc_tag(task_residual_flag,minutes)
  lbl = trim_label(minutes);
  switch task_residual_flag
    case 1,  tag = sprintf('desc-taskpred%s_',lbl);
    case 2,  tag = sprintf('desc-residual%s_',lbl);
    otherwise, tag = sprintf('desc-%s_',lbl);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = trim_tab_desc(task_residual_flag,minutes)
  lbl = trim_label(minutes);
  switch task_residual_flag
    case 1,  desc = sprintf('taskpred%s',lbl);
    case 2,  desc = sprintf('residual%s',lbl);
    otherwise, desc = lbl;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trim_roi_tseries(parms)
  combined_list = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
  flag_list = [0];
  if parms.deconv_flag, flag_list = [0,1,2]; end
  for ti = 1:length(parms.all_trim_minutes)
    trim_min = parms.all_trim_minutes(ti);
    trim_n = parms.trim_ntpoints_list(ti);
    for i = 1:length(combined_list)
      atlas = combined_list{i};
      for task_residual_flag = flag_list
        ptag = prep_desc_tag(task_residual_flag);
        instem = sprintf('%s_atlas-%s_%stseries',parms.outstem,atlas,ptag);
        ttag = trim_desc_tag(task_residual_flag,trim_min);
        outstem = sprintf('%s_atlas-%s_%stseries',parms.outstem,atlas,ttag);
        cmig_trim_fmri_roi_tseries(...
          'indir',parms.outdir,'instem',instem,...
          'outdir',parms.outdir,'outstem',outstem,...
          'trim_ntpoints',trim_n,...
          'trim_rand_flag',parms.trim_rand_flag,...
          'trim_seed',parms.trim_seed,...
          'forceflag',parms.forceflag);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trim_ico_tseries(parms)
  flag_list = [0];
  if parms.deconv_flag, flag_list = [0,1,2]; end
  for ti = 1:length(parms.all_trim_minutes)
    trim_min = parms.all_trim_minutes(ti);
    trim_n = parms.trim_ntpoints_list(ti);
    for ico_order = prep_ico_orders(parms)
      for task_residual_flag = flag_list
        ptag = prep_desc_tag(task_residual_flag);
        instem = sprintf('%s_atlas-ico%d_%stseries',parms.outstem,ico_order,ptag);
        ttag = trim_desc_tag(task_residual_flag,trim_min);
        outstem = sprintf('%s_atlas-ico%d_%stseries',parms.outstem,ico_order,ttag);
        cmig_trim_fmri_ico_tseries(...
          'indir',parms.outdir,'instem',instem,...
          'outdir',parms.outdir,'outstem',outstem,...
          'trim_ntpoints',trim_n,...
          'trim_rand_flag',parms.trim_rand_flag,...
          'trim_seed',parms.trim_seed,...
          'forceflag',parms.forceflag);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_roi_corr(parms)
  indir_tab     = parms.outdir;
  corr_instem   = sprintf('%s_atlas-gpdstdskaseg_corrmat',parms.outstem);
  gpnet_instem  = sprintf('%s_atlas-gpnet_corrmat',parms.outstem);
  gpaseg_instem = sprintf('%s_atlas-gpnetaseg_corrmat',parms.outstem);
  pid = parms.participant_id; sid = parms.session_id;

  % motion tabdat colstem depends on task
  if strcmp(parms.task,'rest')
    mot_colstem = 'mr_y_qc__mot__rsfmri';
  else
    mot_colstem = sprintf('mr_y_qc__mot__tfmri__%s',parms.task);
  end

  % compute correlations for actual signal
  calc_roi_corr(parms,0);
  calc_between_network_corr(parms,0);
  calc_network_subcort_corr(parms,0);

  % write motion tabdat (once, using base gpdstdskaseg_corrmat); a task's is
  %   written after the deconvolution, which its dof comes from
  if strcmp(parms.task,'rest')
    cmig_write_tabdat_fmri_motion('indir',indir_tab,'instem',corr_instem,...
      'outdir',parms.outdir,...
      'task',parms.task,'colstem',mot_colstem,...
      'participant_id',pid,'session_id',sid,'forceflag',parms.forceflag);
  end

  % write tabdat for each signal type; add task/resid if deconvolution was run
  flaglist = 0;
  if parms.deconv_flag
    for flag = [1,2]
      calc_roi_corr(parms,flag);
      calc_between_network_corr(parms,flag);
      calc_network_subcort_corr(parms,flag);
    end
    flaglist = [0,1,2];
  end
  combined_atlases = {'gpaseg','dskaseg','dstaseg','gpdstdskaseg'};
  for flag = flaglist
    tag = corr_tag(flag);
    desc = tag_to_desc(tag);
    col_tag = '';
    switch desc
      case 'taskpred', col_tag = '__task';
      case 'residual', col_tag = '__resid';
    end
    % colstems depend on task
    if strcmp(parms.task,'rest')
      var_colstem  = 'mr_y_rsfmri__var';
      net_colstem  = 'mr_y_rsfmri__corr__gpnet';
      aseg_colstem = 'mr_y_rsfmri__corr__gpnet__aseg';
    else
      var_colstem  = sprintf('mr_y_tfmri__%s__var%s',parms.task,col_tag);
      net_colstem  = sprintf('mr_y_tfmri__%s__corr__gpnet%s',parms.task,col_tag);
      aseg_colstem = sprintf('mr_y_tfmri__%s__corr__gpnet__aseg%s',parms.task,col_tag);
    end
    for ca_i = 1:length(combined_atlases)
      atlas = combined_atlases{ca_i};
      cmig_write_fmri_corrmat('indir',parms.outdir,...
        'instem',bids_tagged_instem(sprintf('%s_atlas-%s_corrmat',parms.outstem,atlas),tag),...
        'outdir',parms.outdir,...
        'task',parms.task,'desc',desc,'atlas',atlas,...
        'participant_id',pid,'session_id',sid,'forceflag',parms.forceflag);
    end

    cmig_write_tabdat_fmri_var('indir',indir_tab,...
      'instem',bids_tagged_instem(corr_instem,tag),...
      'outdir',parms.outdir,...
      'task',parms.task,'desc',desc,...
      'colstem',var_colstem,...
      'tag','','roitypes',{'aseg','dsk','dst','gp'},...
      'participant_id',pid,'session_id',sid,'forceflag',parms.forceflag);

    cmig_write_tabdat_fmri_network_corr('indir',parms.outdir,...
      'instem',bids_tagged_instem(gpnet_instem,tag),'outdir',parms.outdir,...
      'task',parms.task,'desc',desc,...
      'colstem',net_colstem,...
      'participant_id',pid,'session_id',sid,'forceflag',parms.forceflag);

    cmig_write_tabdat_fmri_network_aseg_corr('indir',parms.outdir,...
      'instem',bids_tagged_instem(gpaseg_instem,tag),'outdir',parms.outdir,...
      'task',parms.task,'desc',desc,...
      'colstem',aseg_colstem,...
      'participant_id',pid,'session_id',sid,'forceflag',parms.forceflag);
  end

  % trimmed correlations and corrmat TSVs
  % NOTE: no trimmed tabdat, rest and task alike: the trimmed data stays in
  %   the concatenated data, as in release 7.0
  if ~isempty(parms.trim_ntpoints_list)
    trim_flaglist = [0];
    if parms.deconv_flag, trim_flaglist = [0,1,2]; end
    for ti = 1:length(parms.all_trim_minutes)
      trim_min = parms.all_trim_minutes(ti);
      for flag = trim_flaglist
        calc_roi_corr(parms,flag,trim_min);
        calc_between_network_corr(parms,flag,trim_min);
        calc_network_subcort_corr(parms,flag,trim_min);
      end
      for flag = trim_flaglist
        tdesc = trim_tab_desc(flag,trim_min);
        ttag = trim_corr_tag(flag,trim_min);
        for ca_i = 1:length(combined_atlases)
          atlas = combined_atlases{ca_i};
          cmig_write_fmri_corrmat('indir',parms.outdir,...
            'instem',bids_tagged_instem(sprintf('%s_atlas-%s_corrmat',parms.outstem,atlas),ttag),...
            'outdir',parms.outdir,...
            'task',parms.task,'desc',tdesc,'atlas',atlas,...
            'participant_id',pid,'session_id',sid,'forceflag',parms.forceflag);
        end
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = write_task_motion_tabdat(parms,dparms)
  errmsg = [];
  % the dof is that of the ROI fits, per run and over the concatenated runs,
  %   as 3dDeconvolve gives it for the t statistics: the uncensored time
  %   points less the regressors. With one run, the session's is the run's.
  %   Without a deconvolution there is none, and no dof columns
  nruns = length(parms.fnames_data);
  dofvec = NaN;
  dofvec_run = nan(1,nruns);
  deconv_meanfd = [];
  if ~isempty(dparms)
    [deconv_meanfd,tmp_errmsg] = deconv_motion_stats(dparms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    run_numbers = deconv_run_numbers(dparms,length(dparms.fnames_data));
    for r = 1:length(run_numbers)
      fname_info = sprintf('%s/run%02d_roi_deconv_3dDeconv_stats_info.mat',...
        parms.tmpdir,run_numbers(r));
      [dofvec_run(dparms.run_index(r)),tmp_errmsg] = read_deconv_dof(fname_info);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    end
    if length(run_numbers) > 1
      [dofvec,tmp_errmsg] = read_deconv_dof(sprintf('%s/roi_deconv_3dDeconv_stats_info.mat',...
        parms.tmpdir));
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    else
      dofvec = dofvec_run(dparms.run_index(1));
    end
  end
  cmig_write_tabdat_fmri_motion('indir',parms.outdir,...
    'instem',sprintf('%s_atlas-gpdstdskaseg_corrmat',parms.outstem),...
    'outdir',parms.outdir,...
    'task',parms.task,'colstem',sprintf('mr_y_qc__mot__tfmri__%s',parms.task),...
    'participant_id',parms.participant_id,'session_id',parms.session_id,...
    'dofvec',dofvec,'dofvec_run',dofvec_run,...
    'deconv_meanfd',deconv_meanfd,...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [meanfd,errmsg] = deconv_motion_stats(parms)
  % mean framewise displacement over the frames the session's fit
  %   kept: its runs (an events file, >= min_numTRs frames), without the
  %   dummy TRs, uncensored in fname_censor (FD < censor_thresh_deconv). FD
  %   is computed per run, as for the per-run columns, so each run's first
  %   frame (FD 0 by definition) is left out
  meanfd = NaN; errmsg = [];
  if isempty(parms.fname_censor) || ~exist(parms.fname_censor,'file')
    errmsg = sprintf('deconv censor file not found: %s',parms.fname_censor);
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end
  kept = load(parms.fname_censor);
  fd = []; mask = [];
  for r = 1:length(parms.fnames_motion_1d)
    tmp = load(parms.fnames_motion_1d{r});
    if size(tmp,2) > 6, tmp = tmp(:,2:7); end
    [~,tfd] = mmil_censor_motion(tmp,parms.censor_thresh_deconv,parms.contig_thresh);
    fd = cat(1,fd,tfd(:));
    mask = cat(1,mask,[false;true(length(tfd)-1,1)]);
  end
  if length(kept) ~= length(fd)
    errmsg = sprintf('deconv censor file %s has %d frames, the motion %d',...
      parms.fname_censor,length(kept),length(fd));
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end
  ind = find(mask & kept(:) ~= 0);
  if isempty(ind), return; end
  meanfd = mean(fd(ind));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dof,errmsg] = read_deconv_dof(fname_info)
  % the largest t statistic dof, as the old pipeline took it
  dof = NaN; errmsg = [];
  if ~exist(fname_info,'file')
    errmsg = sprintf('deconv stats info file not found: %s',fname_info);
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    return;
  end
  tmp = load(fname_info,'stats_info');
  si = tmp.stats_info;
  dofs = [];
  for k = 1:length(si)
    if strcmp(si(k).type,'tstat') && numel(si(k).dofs) == 1
      dofs(end+1) = si(k).dofs;
    end
  end
  if ~isempty(dofs), dof = max(dofs); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_between_network_corr(parms,task_residual_flag,trim_min)
  if nargin < 3, trim_min = 0; end
  if trim_min > 0
    tag = trim_corr_tag(task_residual_flag,trim_min);
  else
    tag = corr_tag(task_residual_flag);
  end
  % reads from gpaseg_corrmat
  p = [];
  p.indir = parms.outdir;
  p.outdir = parms.outdir;
  p.instem = bids_tagged_instem(sprintf('%s_atlas-gpaseg_corrmat',parms.outstem),tag);
  p.infostem = sprintf('%s_volinfo',parms.outstem);
  p.outstem = bids_tagged_instem(sprintf('%s_atlas-gpnet_corrmat',parms.outstem),tag);
  p.network_flag = 3;
  p.forceflag = parms.forceflag;
  a = mmil_parms2args(p);
  cmig_calc_fmri_network_corr(a{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_network_subcort_corr(parms,task_residual_flag,trim_min)
  if nargin < 3, trim_min = 0; end
  if trim_min > 0
    tag = trim_corr_tag(task_residual_flag,trim_min);
  else
    tag = corr_tag(task_residual_flag);
  end
  % reads from gpaseg_corrmat
  p = [];
  p.indir = parms.outdir;
  p.outdir = parms.outdir;
  p.instem = bids_tagged_instem(sprintf('%s_atlas-gpaseg_corrmat',parms.outstem),tag);
  p.infostem = sprintf('%s_volinfo',parms.outstem);
  p.outstem = bids_tagged_instem(sprintf('%s_atlas-gpnetaseg_corrmat',parms.outstem),tag);
  p.network_flag = 1;
  if parms.fs_version > 700
    p.roinames = {...
      'Left-Cerebellum-Cortex'...
      'Left-Thalamus'...
      'Left-Caudate'...
      'Left-Putamen'...
      'Left-Pallidum'...
      'Brain-Stem'...
      'Left-Hippocampus'...
      'Left-Amygdala'...
      'Left-Accumbens-area'...
      'Left-VentralDC'...
      'Right-Cerebellum-Cortex'...
      'Right-Thalamus'...
      'Right-Caudate'...
      'Right-Putamen'...
      'Right-Pallidum'...
      'Right-Hippocampus'...
      'Right-Amygdala'...
      'Right-Accumbens-area'...
      'Right-VentralDC'...
    };
  else
    p.roinames = {...
      'Left-Cerebellum-Cortex'...
      'Left-Thalamus-Proper'...
      'Left-Caudate'...
      'Left-Putamen'...
      'Left-Pallidum'...
      'Brain-Stem'...
      'Left-Hippocampus'...
      'Left-Amygdala'...
      'Left-Accumbens-area'...
      'Left-VentralDC'...
      'Right-Cerebellum-Cortex'...
      'Right-Thalamus-Proper'...
      'Right-Caudate'...
      'Right-Putamen'...
      'Right-Pallidum'...
      'Right-Hippocampus'...
      'Right-Amygdala'...
      'Right-Accumbens-area'...
      'Right-VentralDC'...
    };
  end
  p.forceflag = parms.forceflag;
  a = mmil_parms2args(p);
  cmig_calc_fmri_network_corr(a{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tag = corr_tag(task_residual_flag)
  switch task_residual_flag
    case 1,  tag = '_desc-taskpred';
    case 2,  tag = '_desc-residual';
    otherwise, tag = '';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = concat_tfmri_stats(parms)
  % thin wrapper over cmig_concat_tfmri_stats, gathering the arguments the
  %   stand-alone needs from parms
  p = [];
  p.fname_info    = parms.fname_vol_info_ico5;
  p.nruns         = length(parms.fnames_data);
  p.hemilist      = parms.hemilist;
  p.ico_order     = parms.ico_order;
  p.run_tag       = parms.run_tag;
  p.iresp_flag    = parms.iresp_flag;
  p.baseline_flag = parms.baseline_flag;
  p.mask_flag     = true;
  p.outliers_flag = parms.outliers_flag;
  p.stim_fnames   = parms.stim_fnames;
  p.nreg_per_stim = parms.nreg_per_stim;
  p.detrend       = parms.detrend;
  p.fname_motion  = parms.fname_motion;
  p.motion_labels = parms.motion_labels;
  p.deriv_flag    = parms.deriv_flag;
  p.forceflag     = parms.forceflag;
  a = mmil_parms2args(p);
  [~,errcode,errmsg] = cmig_concat_tfmri_stats(parms.participant_id,...
    parms.session_id,parms.output_rootdir,parms.task,a{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = convert_tfmri_stats(parms)
  errcode = 0; errmsg = [];

  % get condnames (conditions + contrasts) from deconvolution output
  fname_condnames = sprintf('%s/deconv_condnames.mat',parms.tmpdir);
  if ~exist(fname_condnames,'file')
    error('condnames file not found: %s',fname_condnames);
  end
  tmp = load(fname_condnames,'condnames');
  condnames = tmp.condnames;
  if isempty(condnames)
    error('no condition names found in %s',fname_condnames);
  end

  % set stat labels
  statlist = {'beta','sem'};
  if parms.iresp_flag
    statlist = [statlist, {'amp','lat'}];
  end

  bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  if ~isempty(parms.run_tag)
    run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
  else
    run_entity = '';
  end

  fprintf('%s: saving tfmri stats in GIfTI format...\n',mfilename);
  % the ico vertices are fsaverage's first 10*4^n+2
  [spatial_ref,ico_desc] = cmig_fsaverage_ref(parms.ico_order);
  sources = proc_bold_uris(parms);
  for hemii = 1:length(parms.hemilist)
    hemi = parms.hemilist{hemii};
    % coverage mask paired with the stats (1 covered, 0 not); float values
    %   with no label table, so .func.gii
    bids_hemi = cmig_get_bids_hemi(hemi);
    fname_mat = sprintf('%s/%s%s_hemi-%s_atlas-ico5_mask.mat',parms.outdir,bids_base,run_entity,bids_hemi);
    fname_gii = sprintf('%s/%s%s_hemi-%s_atlas-ico5_mask.func.gii',parms.outdir,bids_base,run_entity,bids_hemi);
    if ~exist(fname_mat,'file')
      error('coverage mask file not found: %s',fname_mat);
    end
    if ~exist(fname_gii,'file') || parms.forceflag
      tmp_mat = load(fname_mat,'measmat');
      vals = mmil_colvec(tmp_mat.measmat);  % [nverts, 1]
      fname_tmp = sprintf('%s/%s_mask_%s_tmp.mgz',...
        parms.tmpdir,parms.participant_id,hemi);
      fs_save_mgh(vals,fname_tmp);
      fs_mri_convert(fname_tmp,fname_gii);
      delete(fname_tmp);
      strip_gii_metadata(fname_gii);
      fprintf('%s: wrote %s\n',mfilename,fname_gii);
    end
    write_stats_json(regexprep(fname_gii,'\.func\.gii$','.json'),...
      sprintf('Vertices with data in every fitted run (1) or not (0), on %s',ico_desc),...
      sources,spatial_ref,parms.forceflag);
    for c = 1:length(condnames)
      condname = condnames{c};
      for s = 1:length(statlist)
        stat = statlist{s};
        bids_hemi = cmig_get_bids_hemi(hemi);
        fname_mat = sprintf('%s/%s%s_hemi-%s_atlas-ico5_desc-%s_%s.mat',parms.outdir,bids_base,run_entity,bids_hemi,condname,stat);
        fname_gii = sprintf('%s/%s%s_hemi-%s_atlas-ico5_desc-%s_%s.func.gii',parms.outdir,bids_base,run_entity,bids_hemi,condname,stat);
        if ~exist(fname_mat,'file')
          error('deconv stats file not found: %s',fname_mat);
        end
        write_stats_json(regexprep(fname_gii,'\.func\.gii$','.json'),...
          sprintf('%s, on %s',stats_description(stat,condname,parms.baseline_flag),ico_desc),...
          sources,spatial_ref,parms.forceflag);
        if exist(fname_gii,'file') && ~parms.forceflag, continue; end
        % load measmat [1, nverts] and convert to GIfTI via temp mgz
        %   as a column: a [1, nverts] mgz has width 1, and mri_convert
        %   writes a 1-vertex GIfTI
        tmp_mat = load(fname_mat,'measmat');
        vals = mmil_colvec(tmp_mat.measmat);  % [nverts, 1]
        fname_tmp = sprintf('%s/%s_%s_%s_%s_tmp.mgz',...
          parms.tmpdir,parms.participant_id,condname,stat,hemi);
        fs_save_mgh(vals,fname_tmp);
        fs_mri_convert(fname_tmp,fname_gii);
        delete(fname_tmp);
        strip_gii_metadata(fname_gii);
        fprintf('%s: wrote %s\n',mfilename,fname_gii);
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,parms] = prep_deconv(parms)
  errcode = 0; errmsg = [];
  parms.stim_fnames = {};
  parms.fname_motion = [];
  parms.fname_censor = [];
  parms.fname_deriv = [];
  parms.glt_fnames = {};
  parms.glt_labels = {};

  nruns = length(parms.fnames_data);

  % get TR from scan JSON sidecar if not provided by user
  if isempty(parms.TR)
    for r = 1:nruns
      if ~isempty(parms.fnames_info{r}) && exist(parms.fnames_info{r},'file')
        try
          tmp_json = jsondecode(fileread(parms.fnames_info{r}));
          if isfield(tmp_json,'RepetitionTime')
            parms.TR = tmp_json.RepetitionTime;
            break;
          end
        catch
        end
      end
    end
  end
  if isempty(parms.TR)
    errcode = 1; errmsg = 'could not get TR from scan JSON'; return;
  end

  % convert motion TSV files to per-run AFNI .1D format
  parms.fnames_motion_1d = cell(nruns,1);
  for r = 1:nruns
    parms.fnames_motion_1d{r} = prep_motion_1d(parms.fnames_motion{r},parms.tmpdir,r,...
      parms.skipTRs,parms.ntpoints,parms.forceflag);
  end
  [parms.fname_motion,parms.fname_censor,parms.fname_deriv] = prep_motion(...
    parms.fnames_motion_1d,parms.tmpdir,parms.censor_thresh_deconv,...
    parms.contig_thresh,parms.forceflag);

  % the condition and contrast definitions come from check_parms, which needs
  %   them before the time series extraction
  condnames = parms.condnames;

  % convert events TSV files to per-condition stim timing files (one line per run)
  [parms.stim_fnames,parms.stim_fnames_per_run,nevents_by_run_cond,errcode,errmsg] = prep_stims(parms,parms.stim_labels);
  if errcode, return; end

  % condnames is the conditions followed by the contrasts, and only a
  %   condition can be absent, so the contrast entries stay false
  nconds = size(nevents_by_run_cond,2);
  empty_conds = false(1,length(condnames));
  empty_conds(1:nconds) = all(nevents_by_run_cond == 0, 1);
  % per-run flags, with the session run number of each row so that a reader
  %   can match them to a run_tag
  empty_conds_run = false(size(nevents_by_run_cond,1),length(condnames));
  empty_conds_run(:,1:nconds) = (nevents_by_run_cond == 0);
  run_numbers = deconv_run_numbers(parms,size(nevents_by_run_cond,1));

  % a contrast is only as measured as the conditions it references: an absent
  %   condition's beta is zero, so subtracting it is a no-op and the contrast
  %   comes out equal to the remaining main effect rather than obviously wrong
  for i = 1:length(parms.conds_contrast)
    k = nconds + i;
    if k > length(condnames), break; end
    members = parms.conds_contrast(i).groupA;
    if ~isempty(parms.conds_contrast(i).groupB)
      members = [members, parms.conds_contrast(i).groupB];
    end
    if isempty(members), continue; end
    [~,ic] = ismember(members,parms.stim_labels);
    ic = ic(ic>0);
    if isempty(ic), continue; end
    empty_conds(k) = any(empty_conds(ic));
    empty_conds_run(:,k) = any(empty_conds_run(:,ic),2);
  end

  % generate GLT files
  stimdir = sprintf('%s/stims',parms.tmpdir);
  glt_outstem = sprintf('%s/glt',stimdir);
  [parms.glt_fnames,parms.glt_labels] = create_glt_files(glt_outstem,nruns,parms,parms.fname_motion);

  % save condnames for concat_tfmri_stats and convert_tfmri_stats
  fname_condnames = sprintf('%s/deconv_condnames.mat',parms.tmpdir);
  % a file written before empty_conds existed is rewritten, so the readers do
  %   not silently fall back on a stale one
  write_flag = ~exist(fname_condnames,'file') || parms.forceflag;
  if ~write_flag
    tmp_cn = whos('-file',fname_condnames);
    write_flag = ~all(ismember({'empty_conds','empty_conds_run','run_numbers'},...
      {tmp_cn.name}));
  end
  % a file with other condition names (e.g. nback c, now cue) is rewritten too
  if ~write_flag
    tmp_cn = load(fname_condnames,'condnames');
    write_flag = ~isequal(tmp_cn.condnames,condnames);
  end
  if write_flag
    save(fname_condnames,'condnames','empty_conds','empty_conds_run',...
      'run_numbers','-v7.3');
  end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,parms] = prep_design_per_run(parms)
  errcode = 0; errmsg = [];
  nruns = length(parms.fnames_data);
  stimdir = sprintf('%s/stims',parms.tmpdir);
  parms.fnames_design_per_run = cell(nruns,1);
  run_numbers = deconv_run_numbers(parms,nruns);
  for r = 1:nruns
    rnum = run_numbers(r);
    rparms = parms;
    rparms.run_tag = sprintf('run%02d',rnum);
    rparms.fnames_data = {parms.fnames_data{r}};
    % NOTE: fnames_extract must be subset with fnames_data; the deconv
    %   functions size their loop from fnames_data and index fnames_extract,
    %   so leaving it whole makes every run read run 1
    rparms.fnames_extract = {parms.fnames_extract{r}};
    rparms.stim_fnames = parms.stim_fnames_per_run{r};
    [rparms.fname_motion,rparms.fname_censor,rparms.fname_deriv] = prep_motion_run(r,parms);
    glt_outstem = sprintf('%s/run%02d_glt',stimdir,rnum);
    [rparms.glt_fnames,rparms.glt_labels] = create_glt_files(glt_outstem,1,parms,rparms.fname_motion);
    [errcode,errmsg,rparms] = run_deconv_design(rparms);
    if errcode, return; end
    parms.fnames_design_per_run{r} = rparms.fname_design;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [glt_fnames,glt_labels] = create_glt_files(glt_outstem,nruns,parms,fname_motion)
  % motion_flag: 0 = no motion regressors, 1 = 6 (motion only), 2 = 12 (motion + derivatives)
  %   (see mmil_create_glt_files docstring); must match actual design so GLT
  %   zero-padding lines up with the real regressor columns
  if isempty(fname_motion)
    motion_flag = 0;
  elseif parms.deriv_flag
    motion_flag = 2;
  else
    motion_flag = 1;
  end
  [glt_fnames,glt_labels] = mmil_create_glt_files(glt_outstem,...
    'nruns',nruns,...
    'stim_labels',parms.stim_labels,...
    'stim_times_flag',parms.stim_times_flag,...
    'stim_times_model',parms.stim_times_model,...
    'detrend',parms.detrend,...
    'motion_flag',motion_flag,...
    'conds_contrast',parms.conds_contrast,...
    'contrasts_flag',false,...
    'forceflag',parms.forceflag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,parms] = run_deconv_design(parms)
  errcode = 0; errmsg = [];

  % skip if stim files were not prepared
  if isempty(parms.stim_fnames)
    return;
  end

  nruns = length(parms.fnames_data);
  TR = parms.TR;

  % total time points from per-run ico mgz files, already post-skipTRs
  total_nt = 0;
  for r = 1:nruns
    [~,instem] = cmig_get_fileparts(parms.fnames_extract{r});
    fname_mgz = sprintf('%s/%s-ico%d-lh.mgz',parms.tmpdir,instem,parms.ico_order);
    if ~exist(fname_mgz,'file')
      errcode = 1; errmsg = sprintf('ico lh file not found: %s',fname_mgz); return;
    end
    [~,~,~,volsz] = mmil_load_mgh(fname_mgz);
    total_nt = total_nt + volsz(4);
  end
  nodata_nt = total_nt;

  % run design-only deconvolution to produce HRF-convolved design matrix
  if ~isempty(parms.run_tag)
    outstem_design = sprintf('%s/%s_design_3dDeconv',parms.tmpdir,parms.run_tag);
  else
    outstem_design = sprintf('%s/design_3dDeconv',parms.tmpdir);
  end
  parms.fname_design = sprintf('%s.xmat.1D',outstem_design);
  tparms = parms;
  tparms.censor_flag = false;
  tparms.fname_censor = [];
  tparms.nodata_nt = nodata_nt;
  tparms.nodata_tr = parms.TR;
  args = mmil_parms2args(tparms,parms.deconv_tags);
  [~,deconv_errcode] = mmil_run_3dDeconv({},...
    'outstem',outstem_design,...
    'design_flag',true,...
    args{:});
  if deconv_errcode
    errcode = 1; errmsg = 'mmil_run_3dDeconv failed for design'; return;
  end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = run_deconv_ico(parms)
  errcode = 0; errmsg = [];
  if isempty(parms.stim_fnames), return; end
  TR = parms.TR;
  nruns = length(parms.fnames_data);
  deconv_args = mmil_parms2args(parms,parms.deconv_tags);
  fnames_ico_nii_lh = cell(nruns,1);
  fnames_ico_nii_rh = cell(nruns,1);
  for r = 1:nruns
    [~,instem] = cmig_get_fileparts(parms.fnames_extract{r});
    fname_mgz_lh = sprintf('%s/%s-ico%d-lh.mgz',parms.tmpdir,instem,parms.ico_order);
    fname_mgz_rh = sprintf('%s/%s-ico%d-rh.mgz',parms.tmpdir,instem,parms.ico_order);
    if ~exist(fname_mgz_lh,'file')
      errcode = 1; errmsg = sprintf('ico lh file not found: %s',fname_mgz_lh); return;
    end
    if ~exist(fname_mgz_rh,'file')
      errcode = 1; errmsg = sprintf('ico rh file not found: %s',fname_mgz_rh); return;
    end
    if ~isempty(parms.run_tag)
      fnames_ico_nii_lh{r} = sprintf('%s/%s_ico%d_lh.nii.gz',parms.tmpdir,parms.run_tag,parms.ico_order);
      fnames_ico_nii_rh{r} = sprintf('%s/%s_ico%d_rh.nii.gz',parms.tmpdir,parms.run_tag,parms.ico_order);
    else
      fnames_ico_nii_lh{r} = sprintf('%s/run%02d_ico%d_lh.nii.gz',parms.tmpdir,r,parms.ico_order);
      fnames_ico_nii_rh{r} = sprintf('%s/run%02d_ico%d_rh.nii.gz',parms.tmpdir,r,parms.ico_order);
    end
    ico_mgz_to_nii(fname_mgz_lh,fnames_ico_nii_lh{r},TR,parms.forceflag);
    ico_mgz_to_nii(fname_mgz_rh,fnames_ico_nii_rh{r},TR,parms.forceflag);
  end
  for hemii = 1:length(parms.hemilist)
    hemi = parms.hemilist{hemii};
    [errcode,errmsg] = save_mask_ico(parms,hemi);
    if errcode, return; end
    switch hemi
      case 'lh', fnames_ico_nii = fnames_ico_nii_lh;
      case 'rh', fnames_ico_nii = fnames_ico_nii_rh;
    end
    if ~isempty(parms.run_tag)
      outstem_ico = sprintf('%s/%s_ico%d_deconv_%s_3dDeconv',parms.tmpdir,parms.run_tag,parms.ico_order,hemi);
    else
      outstem_ico = sprintf('%s/ico%d_deconv_%s_3dDeconv',parms.tmpdir,parms.ico_order,hemi);
    end
    [tmp_out,deconv_errcode] = mmil_run_3dDeconv(fnames_ico_nii,...
      deconv_args{:},'outstem',outstem_ico);
    if deconv_errcode
      errcode = 1; errmsg = sprintf('mmil_run_3dDeconv failed for ico hemi=%s',hemi); return;
    end
    deconv_brik_to_mgh(tmp_out.glm_results_brik,...
      sprintf('%s.mgh',outstem_ico),parms.forceflag);
    save_deconv_stats_info(tmp_out.glm_results_brik,...
      sprintf('%s_stats_info.mat',outstem_ico),parms.forceflag);
    if parms.baseline_flag
      fname_rbeta = '';
      if parms.reml_flag && isfield(tmp_out,'glm_Rbeta_brik') && ~isempty(tmp_out.glm_Rbeta_brik)
        fname_rbeta = tmp_out.glm_Rbeta_brik;
      end
      save_baseline_ico(parms,sprintf('%s.mgh',outstem_ico),fname_rbeta,hemi,nruns);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = run_deconv_roi(parms)
  errcode = 0; errmsg = [];
  if isempty(parms.stim_fnames), return; end
  TR = parms.TR;
  nruns = length(parms.fnames_data);
  deconv_args = mmil_parms2args(parms,parms.deconv_tags);
  fnames_roi_nii = cell(nruns,1);
  roinames = {};
  for r = 1:nruns
    [~,instem] = cmig_get_fileparts(parms.fnames_extract{r});
    fname_roi_mat = sprintf('%s/%s_roi_data.mat',parms.tmpdir,instem);
    if ~exist(fname_roi_mat,'file')
      errcode = 1; errmsg = sprintf('roi mat file not found: %s',fname_roi_mat); return;
    end
    if ~isempty(parms.run_tag)
      fnames_roi_nii{r} = sprintf('%s/%s_roi_data.nii.gz',parms.tmpdir,parms.run_tag);
    else
      fnames_roi_nii{r} = sprintf('%s/run%02d_roi_data.nii.gz',parms.tmpdir,r);
    end
    roi_mat_to_nii(fname_roi_mat,fnames_roi_nii{r},TR,parms.forceflag);
    if isempty(roinames)
      tmp_roi = load(fname_roi_mat,'roinames');
      roinames = tmp_roi.roinames;
    end
  end
  if ~isempty(parms.run_tag)
    outstem_roi = sprintf('%s/%s_roi_deconv_3dDeconv',parms.tmpdir,parms.run_tag);
  else
    outstem_roi = sprintf('%s/roi_deconv_3dDeconv',parms.tmpdir);
  end
  [tmp_out_roi,deconv_errcode] = mmil_run_3dDeconv(fnames_roi_nii,...
    deconv_args{:},'outstem',outstem_roi);
  if deconv_errcode
    errcode = 1; errmsg = 'mmil_run_3dDeconv failed for roi'; return;
  end
  deconv_brik_to_mat(tmp_out_roi.glm_results_brik,...
    sprintf('%s.mat',outstem_roi),roinames,parms.forceflag);
  save_deconv_stats_info(tmp_out_roi.glm_results_brik,...
    sprintf('%s_stats_info.mat',outstem_roi),parms.forceflag);
  if parms.baseline_flag
    fname_rbeta = '';
    if parms.reml_flag && isfield(tmp_out_roi,'glm_Rbeta_brik') && ~isempty(tmp_out_roi.glm_Rbeta_brik)
      fname_rbeta = tmp_out_roi.glm_Rbeta_brik;
    end
    save_baseline_roi(parms,sprintf('%s.mat',outstem_roi),fname_rbeta,roinames,nruns);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = run_deconv_vox(parms)
  errcode = 0; errmsg = [];
  if isempty(parms.stim_fnames), return; end
  TR = parms.TR;
  nruns = length(parms.fnames_data);
  deconv_args = mmil_parms2args(parms,parms.deconv_tags);
  fnames_vox_nii = cell(nruns,1);
  for r = 1:nruns
    [~,instem] = cmig_get_fileparts(parms.fnames_extract{r});
    if ~isempty(parms.run_tag)
      fnames_vox_nii{r} = sprintf('%s/%s_%s_TR.nii.gz',parms.tmpdir,parms.run_tag,instem);
    else
      fnames_vox_nii{r} = sprintf('%s/run%02d_%s_TR.nii.gz',parms.tmpdir,r,instem);
    end
    if ~exist(fnames_vox_nii{r},'file') || parms.forceflag
      [~,outstem] = fileparts(fnames_vox_nii{r});
      outstem = regexprep(outstem,'\.nii$','');
      cmd = sprintf('3dTcat -prefix %s.nii.gz -session %s %s', ...
        outstem, parms.tmpdir, parms.fnames_extract{r});
      mmil_unix(cmd);
      mmil_unix(sprintf('cd %s && 3drefit -TR %.6f %s', parms.tmpdir, TR, fnames_vox_nii{r}));
    end
  end
  if ~isempty(parms.run_tag)
    outstem_vox = sprintf('%s/%s_vox_deconv_3dDeconv',parms.tmpdir,parms.run_tag);
  else
    outstem_vox = sprintf('%s/vox_deconv_3dDeconv',parms.tmpdir);
  end
  [tmp_out_vox,deconv_errcode] = mmil_run_3dDeconv(fnames_vox_nii,...
    deconv_args{:},'outstem',outstem_vox);
  if deconv_errcode
    errcode = 1; errmsg = 'mmil_run_3dDeconv failed for voxelwise'; return;
  end
  deconv_brik_to_nii(tmp_out_vox.glm_results_brik,...
    sprintf('%s.nii.gz',outstem_vox),parms.forceflag);
  save_deconv_stats_info(tmp_out_vox.glm_results_brik,...
    sprintf('%s_stats_info.mat',outstem_vox),parms.forceflag);
  deconv_brik_to_mgh(tmp_out_vox.glm_results_brik,...
    sprintf('%s.mgh',outstem_vox),parms.forceflag);
  fname_rbeta = '';
  if parms.reml_flag && isfield(tmp_out_vox,'glm_Rbeta_brik') && ~isempty(tmp_out_vox.glm_Rbeta_brik)
    fname_rbeta = tmp_out_vox.glm_Rbeta_brik;
  end
  [errcode,errmsg] = write_tfmri_stats_vox(parms,sprintf('%s.mgh',outstem_vox),fname_rbeta,nruns);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = write_tfmri_stats_vox(parms,fname_mgh,fname_rbeta,nruns)
  % write the voxelwise beta and sem of each condition and contrast to outdir,
  %   one nii.gz per condition and stat on the proc func grid, named as the
  %   ico5 outputs of cmig_concat_tfmri_stats without the hemi and atlas
  %   entities, so they can be resampled to atlas and concatenated
  % NOTE: the same column layout, censoring and post-normalization as
  %   cmig_concat_tfmri_stats: condnames index the GLT block after
  %   cmig_deconv_col_offset columns, outlier and empty conditions are left
  %   NaN, and with baseline_flag the beta and sem are scaled to percent of the
  %   model's polort-0 baseline, which is also written, to tmpdir
  % NOTE: impulse response amp and lat (iresp_flag) are not written here
  errcode = 0; errmsg = [];
  fname_condnames = sprintf('%s/deconv_condnames.mat',parms.tmpdir);
  if ~exist(fname_condnames,'file')
    errcode = 1;
    errmsg = sprintf('condnames file not found: %s',fname_condnames);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  tmp = load(fname_condnames,'condnames');
  condnames = tmp.condnames;
  npars = length(condnames);
  if npars == 0, return; end
  empty_conds = cmig_get_empty_conds(fname_condnames,parms.run_tag);
  if length(empty_conds) ~= npars, empty_conds = false(1,npars); end
  excluded = false(1,npars);
  if parms.outliers_flag
    if ~isempty(parms.run_tag)
      fname_outliers = sprintf('%s/%s_deconv_outliers.mat',parms.tmpdir,parms.run_tag);
    else
      fname_outliers = sprintf('%s/deconv_outliers.mat',parms.tmpdir);
    end
    if exist(fname_outliers,'file')
      tmp = load(fname_outliers,'excluded');
      if length(tmp.excluded) == npars, excluded = tmp.excluded; end
    end
  end
  col_offset = cmig_deconv_col_offset(parms,nruns);

  bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  if ~isempty(parms.run_tag)
    run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
  else
    run_entity = '';
  end

  [vol,M,~,volsz] = mmil_load_mgh(fname_mgh);
  nvox = prod(volsz(1:3));
  vol = reshape(vol,[nvox,volsz(4)]);
  if size(vol,2) < 1+col_offset+npars*2
    errcode = 1;
    errmsg = sprintf('voxelwise deconv output %s has %d frames, fewer than the %d conditions need',...
      fname_mgh,size(vol,2),npars);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % model-based polort-0 baseline, averaged over runs, as save_baseline_ico
  baseline = [];
  if parms.baseline_flag
    baseline = zeros(nvox,1,'single');
    if ~isempty(fname_rbeta) && exist(fname_rbeta,'file')
      fname_tmp = strrep(fname_rbeta,'+orig.BRIK','_rbetatmp.mgz');
      mmil_BRIK2mgh(fname_rbeta,fname_tmp,parms.forceflag);
      [rvol,~,~,rvolsz] = mmil_load_mgh(fname_tmp);
      delete(fname_tmp);
      rvol = reshape(rvol,[nvox,rvolsz(4)]);
      for r = 1:nruns
        sb = (r-1)*(parms.detrend+1) + 1;
        if sb <= size(rvol,2), baseline = baseline + single(rvol(:,sb)); end
      end
    else
      for r = 1:nruns
        col = 2 + (r-1)*(parms.detrend+1)*2;
        if col <= size(vol,2), baseline = baseline + single(vol(:,col)); end
      end
    end
    baseline = baseline / nruns;
    fname_out = sprintf('%s/%s%s_baseline.nii.gz',parms.tmpdir,bids_base,run_entity);
    write_vol_nii(baseline,volsz,M,fname_out,parms);
  end

  fprintf('%s: saving voxelwise tfmri stats...\n',mfilename);
  statlist = {'beta','sem'};
  % the proc bold of the runs fitted; the outputs are on their grid
  sources = proc_bold_uris(parms);
  for c = 1:npars
    is_outlier = excluded(c) || empty_conds(c);
    for s = 1:2
      fname_out = sprintf('%s/%s%s_desc-%s_%s.nii.gz',...
        parms.outdir,bids_base,run_entity,condnames{c},statlist{s});
      write_stats_json(regexprep(fname_out,'\.nii\.gz$','.json'),...
        stats_description(statlist{s},condnames{c},parms.baseline_flag),...
        sources,sources{1},parms.forceflag);
      if exist(fname_out,'file') && ~parms.forceflag, continue; end
      vals = NaN(nvox,1,'single');
      if ~is_outlier
        vals = single(vol(:,1+col_offset+(c-1)*2+s));
        if ~isempty(baseline)
          vals = 100 * vals ./ (baseline + eps);
        end
      end
      write_vol_nii(vals,volsz,M,fname_out,parms);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_vol_nii(vals,volsz,M,fname_out,parms)
  if exist(fname_out,'file') && ~parms.forceflag, return; end
  [~,fstem] = cmig_get_fileparts(fname_out);
  fname_tmp = sprintf('%s/%s_tmp.mgz',parms.tmpdir,fstem);
  fs_save_mgh(reshape(vals,volsz(1:3)),fname_tmp,M);
  fs_mri_convert(fname_tmp,fname_out,'forceflag',1);
  delete(fname_tmp);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function strip_gii_metadata(fname)
  % remove the UserName and Date metadata that mri_convert writes into a GIfTI
  str = fileread(fname);
  str = regexprep(str,...
    '\s*<MD>\s*<Name><!\[CDATA\[(UserName|Date)\]\]></Name>\s*<Value>.*?</Value>\s*</MD>','');
  fid = fopen(fname,'w');
  if fid < 0, error('failed to open %s for writing',fname); end
  fwrite(fid,str,'char');
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uris = proc_bold_uris(parms)
  % BIDS URIs of the proc bold of the runs fitted, in the proc root
  %   (dataset cmig_proc)
  uris = cell(1,length(parms.fnames_data));
  for r = 1:length(parms.fnames_data)
    [~,fstem,fext] = cmig_get_fileparts(parms.fnames_data{r});
    uris{r} = sprintf('bids:cmig_proc:%s/%s/func/%s%s',...
      parms.participant_id,parms.session_id,fstem,fext);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = stats_description(stat,condname,baseline_flag)
  switch stat
    case 'beta', desc = sprintf('Task fMRI beta for %s',condname);
    case 'sem',  desc = sprintf('Standard error of the task fMRI beta for %s',condname);
    case 'amp',  desc = sprintf('Task fMRI impulse response amplitude for %s',condname);
    case 'lat',  desc = sprintf('Task fMRI impulse response latency for %s',condname);
    otherwise,   desc = sprintf('Task fMRI %s for %s',stat,condname);
  end
  if baseline_flag && ismember(stat,{'beta','sem'})
    desc = sprintf('%s, in percent of the baseline signal',desc);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stats_json(fname_json,description,sources,spatial_reference,forceflag)
  if exist(fname_json,'file') && ~forceflag, return; end
  info = [];
  info.Description = description;
  info.Sources = sources;
  info.SpatialReference = spatial_reference;
  % NOTE: PrettyPrint needs R2021a (9.10)
  if verLessThan('matlab','9.10')
    json = jsonencode(info);
  else
    json = jsonencode(info,'PrettyPrint',true);
  end
  fid = fopen(fname_json,'wt');
  if fid < 0, error('failed to open %s for writing',fname_json); end
  fprintf(fid,'%s',json);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = run_deconv_per_run(parms)
  errcode = 0; errmsg = [];
  nruns = length(parms.fnames_data);
  stimdir = sprintf('%s/stims',parms.tmpdir);
  % r indexes the runs handed to the deconvolution; rnum is the run it was in
  %   the session, which differs when a run was excluded for having no events
  run_numbers = deconv_run_numbers(parms,nruns);
  for r = 1:nruns
    rnum = run_numbers(r);
    fprintf('%s: running per-run deconvolution for run %d/%d...\n',mfilename,rnum,nruns);

    % build per-run rparms: single-run data, stim files, motion files, and output stem
    rparms = parms;
    rparms.run_tag = sprintf('run%02d',rnum);
    % the columns take the release slot, the runs being in slot order; the
    %   file names and temp files keep the BIDS run
    rparms.col_run_tag = sprintf('r%02d',r);
    rparms.outstem = sprintf('%s_run-%02d',parms.outstem,rnum);
    rparms.fnames_data = {parms.fnames_data{r}};
    % NOTE: fnames_extract must be subset with fnames_data; the deconv
    %   functions size their loop from fnames_data and index fnames_extract,
    %   so leaving it whole makes every run read run 1
    rparms.fnames_extract = {parms.fnames_extract{r}};
    rparms.stim_fnames = parms.stim_fnames_per_run{r};
    [rparms.fname_motion,rparms.fname_censor,rparms.fname_deriv] = prep_motion_run(r,parms);

    % regenerate GLT files with nruns=1 (sets correct number of baseline polynomial regressors)
    glt_outstem = sprintf('%s/run%02d_glt',stimdir,rnum);
    [rparms.glt_fnames,rparms.glt_labels] = create_glt_files(glt_outstem,1,parms,rparms.fname_motion);

    % run per-run design deconvolution (nodata_nt computed from single-run mgz)
    [errcode,tmp_errmsg,rparms] = run_deconv_design(rparms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end

    if parms.roi_flag
      [errcode,tmp_errmsg] = run_deconv_roi(rparms);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, return; end
      errmsg = cmig_append_errmsg(errmsg,postnorm_deconv_roi(rparms));
      errmsg = cmig_append_errmsg(errmsg,write_deconv_tabdat(rparms));
      cmig_write_tabdat_tfmri_fstat('indir',rparms.tmpdir,...
        'instem','roi_deconv_3dDeconv_norm',...
        'outdir',rparms.outdir,...
        'colstem','mr_y_tfmri','task',rparms.task,'run_tag',rparms.run_tag,...
        'col_run_tag',rparms.col_run_tag,...
        'participant_id',rparms.participant_id,'session_id',rparms.session_id,...
        'forceflag',rparms.forceflag);

    end
    if parms.ico_flag
      [errcode,tmp_errmsg] = run_deconv_ico(rparms);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, return; end
      [errcode,tmp_errmsg] = concat_tfmri_stats(rparms);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, return; end
      [errcode,tmp_errmsg] = convert_tfmri_stats(rparms);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, return; end
    end
    if parms.vox_flag
      [errcode,tmp_errmsg] = run_deconv_vox(rparms);
      errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
      if errcode, return; end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = copy_single_run_outputs(parms)
  % run_deconv_concat_runs does nothing for a single run, so nothing writes
  %   the bare-stem outputs that the concat-format files and the tabulated
  %   data are read from, and the session drops out of tfmri_stats. The one
  %   run's fit IS the session result, as the release treated it, so its
  %   outputs are copied rather than refitted
  % NOTE: copied rather than produced by calling the writers again with the
  %   run entity suppressed. concat_tfmri_stats reads the ico baseline file
  %   under the run entity and, finding none, leaves baseline_hemi empty and
  %   skips the post-normalization without complaint, so that route would
  %   quietly write un-normalized session betas
  errcode = 0; errmsg = [];
  if length(parms.fnames_data) ~= 1, return; end
  run_numbers = deconv_run_numbers(parms,1);
  tag = sprintf('_run-%02d',run_numbers(1));
  % restricted to this task: outdir is shared across tasks, and another task
  %   in the session may legitimately have per-run outputs of its own
  d = dir(sprintf('%s/*task-%s*%s*',parms.outdir,parms.task,tag));
  ncopied = 0;
  for k = 1:length(d)
    if d(k).isdir, continue; end
    % the per-run deconv and fstat TSVs carry the per-run columns (r01);
    %   their session TSVs are written below
    if ~isempty(regexp(d(k).name,'_(deconv|fstat)\.tsv$','once')), continue; end
    fname_src = sprintf('%s/%s',parms.outdir,d(k).name);
    fname_out = sprintf('%s/%s',parms.outdir,strrep(d(k).name,tag,''));
    if strcmp(fname_src,fname_out), continue; end
    if exist(fname_out,'file') && ~parms.forceflag, continue; end
    [success,msg] = copyfile(fname_src,fname_out);
    if ~success
      errcode = 1;
      errmsg = sprintf('failed to copy %s to %s: %s',fname_src,fname_out,msg);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    ncopied = ncopied + 1;
  end
  fprintf('%s: one deconvolved run, copied %d outputs to the session stem\n',...
    mfilename,ncopied);

  % the session deconv and fstat TSVs: the run's values under the session
  %   column names
  if parms.roi_flag
    rparms = parms;
    rparms.run_tag = sprintf('run%02d',run_numbers(1));
    [~,~,fname_src] = deconv_roi_fname(rparms);
    sparms = parms;
    sparms.run_tag = '';
    sparms.col_run_tag = '';
    [~,~,fname_out] = deconv_roi_fname(sparms);
    if ~exist(fname_src,'file')
      errcode = 1;
      errmsg = sprintf('per-run roi deconv not found: %s',fname_src);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    if ~exist(fname_out,'file') || parms.forceflag
      [success,msg] = copyfile(fname_src,fname_out);
      if ~success
        errcode = 1;
        errmsg = sprintf('failed to copy %s to %s: %s',fname_src,fname_out,msg);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
    end
    errmsg = cmig_append_errmsg(errmsg,write_deconv_tabdat(sparms));
    cmig_write_tabdat_tfmri_fstat('indir',sparms.tmpdir,...
      'instem','roi_deconv_3dDeconv_norm',...
      'outdir',sparms.outdir,...
      'colstem','mr_y_tfmri','task',sparms.task,...
      'participant_id',sparms.participant_id,'session_id',sparms.session_id,...
      'forceflag',sparms.forceflag);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = run_deconv_concat_runs(parms)
  errcode = 0; errmsg = [];
  if length(parms.fnames_data) <= 1, return; end

  if parms.roi_flag
    % run roi deconvolution on concatenated data
    [errcode,tmp_errmsg] = run_deconv_roi(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end
    errmsg = cmig_append_errmsg(errmsg,postnorm_deconv_roi(parms));
    errmsg = cmig_append_errmsg(errmsg,write_deconv_tabdat(parms));
    cmig_write_tabdat_tfmri_fstat('indir',parms.tmpdir,...
      'instem','roi_deconv_3dDeconv_norm',...
      'outdir',parms.outdir,...
      'colstem','mr_y_tfmri','task',parms.task,...
      'participant_id',parms.participant_id,'session_id',parms.session_id,...
      'forceflag',parms.forceflag);

  end

  if parms.ico_flag
    % run ico surface deconvolution on concatenated data
    [errcode,tmp_errmsg] = run_deconv_ico(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end

    % concatenate ico surface stats across runs
    [errcode,tmp_errmsg] = concat_tfmri_stats(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end

    % convert ico surface stats to GIfTI format
    [errcode,tmp_errmsg] = convert_tfmri_stats(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end
  end

  if parms.vox_flag
    % run voxelwise deconvolution on concatenated data
    [errcode,tmp_errmsg] = run_deconv_vox(parms);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, return; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,parms] = filter_motion_tsvs(parms)
  % notch filter motion estimates for aliased respiration signals
  %   writes filtered copies of the motion TSV files to tmpdir and replaces
  %   parms.fnames_motion so all downstream consumers inherit filtered motion
  %   the first skipTRs rows are left unfiltered, since every consumer drops
  %   them; filtering after the trim keeps the dummy TRs out of the filter
  errcode = 0; errmsg = [];
  nruns = length(parms.fnames_motion);
  fnames_filt = cell(nruns,1);
  for r=1:nruns
    fname_in = parms.fnames_motion{r};
    [~,fstem] = fileparts(fname_in);
    fname_out = sprintf('%s/%s_respfilt_trim.tsv',parms.tmpdir,fstem);
    if ~exist(fname_out,'file') || parms.forceflag
      fprintf('%s: filtering respiratory signals from motion estimates for run %d...\n',...
        mfilename,r);
      raw = mmil_readtext(fname_in,'\t');
      col_names = raw(1,:);
      vals = cell2mat(raw(2:end,:));
      % motion parameter columns follow mmil_load_motion_tsv conventions:
      %   first column is frame index; if more than 7 columns, params are 3:8
      if size(vals,2)==7
        ind_mp = 2:7;
      elseif size(vals,2)>7
        ind_mp = 3:8;
      else
        errcode = 1;
        errmsg = sprintf('unexpected number of columns (%d) in motion file %s',...
          size(vals,2),fname_in);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
      ind_t = parms.skipTRs+1:size(vals,1);
      vals(ind_t,ind_mp) = abcd_resp_filter_motion(vals(ind_t,ind_mp),parms.TR,...
        parms.resp_low,parms.resp_high);
      fid = fopen(fname_out,'wt');
      if fid<0
        errcode = 1;
        errmsg = sprintf('failed to open %s for writing',fname_out);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
      fprintf(fid,'%s\n',strjoin(col_names,'\t'));
      for t=1:size(vals,1)
        fprintf(fid,'%g',vals(t,1));
        fprintf(fid,'\t%g',vals(t,2:end));
        fprintf(fid,'\n');
      end
      fclose(fid);
    end
    fnames_filt{r} = fname_out;
  end
  parms.fnames_motion = fnames_filt;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_1d = prep_motion_1d(fname_tsv,outdir,run_idx,skipTRs,ntpoints,forceflag)
  [~,fstem] = fileparts(fname_tsv);
  fname_1d = sprintf('%s/run%02d_%s.1D',outdir,run_idx,fstem);
  if exist(fname_1d,'file') && ~forceflag, return; end
  motion = mmil_load_motion_tsv(fname_tsv);
  % the TRs normalize_runs keeps
  motion = motion(skipTRs+1:min(end,skipTRs+ntpoints),:);
  fid = fopen(fname_1d,'wt');
  if fid < 0, error('failed to open %s for writing',fname_1d); end
  for t = 1:size(motion,1)
    fprintf(fid,'%g ',motion(t,:));
    fprintf(fid,'\n');
  end
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_motion,fname_censor,fname_deriv] = prep_motion(fnames_1d,outdir,censor_thresh,contig_thresh,forceflag,prefix)
  if ~exist('prefix','var') || isempty(prefix), prefix = 'concat'; end
  fname_motion = sprintf('%s/%s_motion.1D',outdir,prefix);
  fname_censor = sprintf('%s/%s_censor.1D',outdir,prefix);
  fname_deriv  = sprintf('%s/%s_deriv.1D',outdir,prefix);
  if exist(fname_motion,'file') && exist(fname_censor,'file') && exist(fname_deriv,'file') && ~forceflag
    return;
  end
  % prep_motion_1d already wrote these in mmil_load_motion_tsv's output order,
  %   dx,dy,dz,rx,ry,rz with clockwise rotations, which is what
  %   mmil_censor_motion, mmil_motion_stats and mmil_write_motion_1D all take.
  %   Reading them back with mmil_load_motion_1D would apply that reorder and
  %   the rotation sign flip a second time, leaving rotations in columns 1 to 3
  %   so that mmil_motion_stats scaled the translations by the head radius.
  motion_data = [];
  run_ind = [];
  for r = 1:length(fnames_1d)
    tmp = load(fnames_1d{r});
    if size(tmp,2) > 6, tmp = tmp(:,2:7); end
    motion_data = cat(1,motion_data,tmp);
    run_ind = cat(1,run_ind,r*ones(size(tmp,1),1));
  end
  % censoring and derivatives from the motion as estimated
  censor_tseries = mmil_censor_motion(motion_data,censor_thresh,contig_thresh);
  [~,~,deriv_data] = mmil_motion_stats(motion_data,50,true,[],false);
  % regressors demeaned per run over the kept frames, so the polort-0
  %   coefficient is the kept-frame signal level
  motion_data = demean_kept(motion_data,run_ind,~censor_tseries);
  deriv_data = demean_kept(deriv_data,run_ind,~censor_tseries);
  mmil_write_motion_1D(motion_data,fname_motion,'forceflag',1);
  fid = fopen(fname_censor,'w');
  if fid < 0, error('failed to open %s for writing',fname_censor); end
  fprintf(fid,'%d\n',~censor_tseries);
  fclose(fid);
  mmil_write_motion_1D(deriv_data,fname_deriv,'forceflag',1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = demean_kept(data,run_ind,kept)
  % subtract each run's mean over its kept frames; a run with none kept
  %   uses all its frames
  kept = logical(kept(:));
  for r = unique(run_ind(:))'
    ind = (run_ind(:) == r);
    ind_mean = ind & kept;
    if ~any(ind_mean), ind_mean = ind; end
    data(ind,:) = bsxfun(@minus,data(ind,:),mean(data(ind_mean,:),1));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_motion,fname_censor,fname_deriv] = prep_motion_run(r,parms)
  % labelled by the run it was in the session, so a changed run selection does
  %   not collide with temp files from an earlier pass
  run_numbers = deconv_run_numbers(parms,length(parms.fnames_data));
  prefix = sprintf('run%02d',run_numbers(r));
  [fname_motion,fname_censor,fname_deriv] = prep_motion(...
    {parms.fnames_motion_1d{r}},parms.tmpdir,parms.censor_thresh_deconv,...
    parms.contig_thresh,parms.forceflag,prefix);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ico_mgz_to_nii(fname_mgz,fname_nii,TR,forceflag)
  % the dummy TRs are already gone, removed in normalize_runs
  if exist(fname_nii,'file') && ~forceflag, return; end
  [vol,M,mr_parms,volsz] = mmil_load_mgh(fname_mgz);
  vol = reshape(vol,[2,3,volsz(1)/6,volsz(4)]);
  mr_parms(1) = TR * 1000;
  fs_save_nifti(vol,fname_nii,M,mr_parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_mat_to_nii(fname_mat,fname_nii,TR,forceflag)
  % the dummy TRs are already gone, removed in normalize_runs
  if exist(fname_nii,'file') && ~forceflag, return; end
  tmp = load(fname_mat,'data');
  data = tmp.data';  % [nrois, ntpoints]
  nrois = size(data,1);
  ntpoints = size(data,2);
  % zero-pad nrois to multiple of 6 and reshape to pseudo-volume so AFNI
  % treats it as a BRIK rather than a 1D file (ny=1,nz=1 triggers 1D mode)
  pad = ceil(nrois/6)*6;
  data_pad = zeros(pad,ntpoints);
  data_pad(1:nrois,:) = data;
  vol = reshape(data_pad,[2,3,pad/6,ntpoints]);
  mr_parms = zeros(4,1);
  mr_parms(1) = TR*1000;
  fs_save_nifti(vol,fname_nii,[],mr_parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deconv_brik_to_mgh(fname_brik,fname_mgh,forceflag)
  if exist(fname_mgh,'file') && ~forceflag, return; end
  fname_tmp = strrep(fname_brik,'+orig.BRIK','_tmp.mgz');
  mmil_BRIK2mgh(fname_brik,fname_tmp,forceflag);
  [vol,M,mr_parms] = mmil_load_mgh(fname_tmp);
  vol = convert_tstat_to_sem(vol);
  fs_save_mgh(vol,fname_mgh,M,mr_parms);
  delete(fname_tmp);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deconv_brik_to_nii(fname_brik,fname_nii,forceflag)
  if exist(fname_nii,'file') && ~forceflag, return; end
  fname_tmp = strrep(fname_brik,'+orig.BRIK','_tmp.mgz');
  fname_tmp_mgh = strrep(fname_brik,'+orig.BRIK','_tmpsem.mgh');
  mmil_BRIK2mgh(fname_brik,fname_tmp,forceflag);
  [vol,M,mr_parms] = mmil_load_mgh(fname_tmp);
  vol = convert_tstat_to_sem(vol);
  fs_save_mgh(vol,fname_tmp_mgh,M,mr_parms);
  fs_mri_convert(fname_tmp_mgh,fname_nii,'forceflag',forceflag);
  delete(fname_tmp);
  delete(fname_tmp_mgh);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deconv_brik_to_mat(fname_brik,fname_mat,roinames,forceflag)
  if exist(fname_mat,'file') && ~forceflag, return; end
  fname_tmp = strrep(fname_brik,'+orig.BRIK','_tmp.mgz');
  mmil_BRIK2mgh(fname_brik,fname_tmp,forceflag);
  [vol,~,~,volsz] = mmil_load_mgh(fname_tmp);
  nroi = length(roinames);
  nframes = volsz(4);
  vol = reshape(vol,[prod(volsz(1:3)),nframes]);
  vol = vol(1:nroi,:);
  vol = convert_tstat_to_sem(vol);
  save(fname_mat,'vol','roinames','-v7.3');
  delete(fname_tmp);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_deconv_stats_info(fname_brik,fname_info,forceflag)
  if ~exist('forceflag','var') || isempty(forceflag), forceflag = false; end
  if exist(fname_info,'file') && ~forceflag, return; end
  if ~exist(fname_brik,'file')
    error('BRIK file not found: %s',fname_brik);
  end
  stats_info = mmil_read_HEAD_statsinfo(fname_brik);
  save(fname_info,'stats_info');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stim_fnames,stim_fnames_per_run,nevents_by_run_cond,errcode,errmsg] = prep_stims(parms,stim_labels)
  % the per-run timing files come from events_to_timing() in
  %   cmig_abcd_events_to_timing.py, the python version of the
  %   abcd_extract_eprime_* timing files, given this module's skipTRs, each
  %   run's TR count and minfrac; they are joined here into the per-run and
  %   concatenated files 3dDeconvolve reads
  errcode = 0; errmsg = [];
  stim_fnames = {};
  stim_fnames_per_run = {};
  nevents_by_run_cond = [];

  if isempty(stim_labels)
    errcode = 1;
    errmsg = 'stim_labels is empty; cannot create stim timing files';
    return;
  end
  nconds = length(stim_labels);
  nruns = length(parms.fnames_data);

  stimdir = sprintf('%s/stims',parms.tmpdir);
  etdir = sprintf('%s/events_to_timing',stimdir);
  mmil_mkdir(etdir);

  [fname_py,errcode,errmsg] = find_python_script(parms,'cmig_abcd_events_to_timing.py');
  if errcode, return; end

  % TRs per run in the data deconvolved, the skipTRs frames already dropped
  ntpoints_per_run = zeros(nruns,1);
  for r = 1:nruns
    motion_r = mmil_load_motion_1D(parms.fnames_motion_1d{r});
    ntpoints_per_run(r) = size(motion_r,1);
  end

  % write every format for each run; a run without an events file has none
  for r = 1:nruns
    fname_events = parms.fnames_events{r};
    if isempty(fname_events) || ~exist(fname_events,'file'), continue; end
    cmd = sprintf(['python3 %s --fname_events %s --TR %g --numTRs %d',...
      ' --minfrac %g --skipTRs %d --block_gap %g --task %s',...
      ' --outdir %s --outstem run%02d --forceflag 1 --verbose 0'],...
      fname_py,fname_events,parms.TR,ntpoints_per_run(r)+parms.skipTRs,...
      parms.minfrac,parms.skipTRs,parms.block_gap,parms.task,etdir,r);
    [status,result] = mmil_unix(cmd);
    if ~isempty(result), fprintf('%s',result); end
    if status
      errcode = 1;
      errmsg = sprintf('events_to_timing failed for %s:\n%s',fname_events,result);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
  end

  % the reading and writing below call error() on a wrong-length .1D or a
  %   file that cannot be opened; caught here, so the module still writes
  %   its status with the message
  try
    % read back, per run and condition, the times for the chosen format, the
    %   block onsets (counting the events) and the mean block duration; a
    %   condition with no events, or a run with no responses, which
    %   events_to_timing() writes nothing for, has none
    times_by_run_cond = cell(nruns,nconds);
    dur_by_run_cond = nan(nruns,nconds);
    nevents_by_run_cond = zeros(nruns,nconds);
    for r = 1:nruns
      for c = 1:nconds
        fstem = sprintf('%s/run%02d_%s',etdir,r,stim_labels{c});
        onsets = read_stim_times([fstem '_block.txt']);
        nevents_by_run_cond(r,c) = length(onsets);
        if isempty(onsets), continue; end
        dur_by_run_cond(r,c) = read_stim_times([fstem '_block_dur.txt']);
        switch parms.stim_times_flag
          case 0
            times_by_run_cond{r,c} = read_stim_vec([fstem '.1D'],ntpoints_per_run(r));
          case 1
            times_by_run_cond{r,c} = read_stim_times([fstem '.txt']);
          otherwise
            times_by_run_cond{r,c} = onsets;
        end
      end
    end

    % a condition with no event is not measured, so its beta must be reported as
    %   missing rather than as the zero the all-zero regressor gives; counted per
    %   run, since the per-run outputs carry the same zeros as the session ones
    for c = 1:nconds
      if all(nevents_by_run_cond(:,c) == 0)
        fprintf('%s: WARNING: no %s events in any run; its beta will be reported as missing\n',...
          mfilename,stim_labels{c});
      end
    end

    % the fit across runs takes one duration per condition, the mean of the
    %   runs'; runs whose durations disagree by too much stop it, as
    %   mmil_3dDeconv did
    if parms.stim_times_flag == 2
      for c = 1:nconds
        durs = dur_by_run_cond(:,c);
        durs = durs(~isnan(durs));
        if length(durs) > 1 && std(durs) >= parms.max_dur_std
          errcode = 1;
          errmsg = sprintf(['not processed: %s block durations differ across runs (%s s):',...
            ' std %0.2f s >= %g s'],stim_labels{c},...
            strtrim(sprintf('%0.2f ',durs)),std(durs),parms.max_dur_std);
          fprintf('%s: ERROR: %s\n',mfilename,errmsg);
          return;
        end
      end
    end

    switch parms.stim_times_flag
      case 0
        fext = '.1D';
      case 1
        fext = '.txt';
      otherwise
        fext = '_block.txt';
    end
    stim_fnames = cell(nconds,1);
    stim_fnames_per_run = cell(nruns,1);
    for r = 1:nruns
      stim_fnames_per_run{r} = cell(nconds,1);
    end
    for c = 1:nconds
      % all runs: .1D vectors stacked, AFNI local_times one line per run
      fname_stim = sprintf('%s/concat_stim_%s%s',stimdir,stim_labels{c},fext);
      stim_fnames{c} = fname_stim;
      if ~exist(fname_stim,'file') || parms.forceflag
        fid = fopen(fname_stim,'wt');
        if fid < 0, error('failed to open %s for writing',fname_stim); end
        for r = 1:nruns
          write_stim_run(fid,times_by_run_cond{r,c},parms.stim_times_flag,ntpoints_per_run(r));
        end
        fclose(fid);
      end
      if parms.stim_times_flag == 2
        durs = dur_by_run_cond(:,c);
        write_stim_dur(strrep(fname_stim,'.txt','_dur.txt'),durs(~isnan(durs)),parms.forceflag);
      end
      % each run on its own
      for r = 1:nruns
        fname_stim_run = sprintf('%s/run%02d_stim_%s%s',stimdir,r,stim_labels{c},fext);
        stim_fnames_per_run{r}{c} = fname_stim_run;
        if ~exist(fname_stim_run,'file') || parms.forceflag
          fid = fopen(fname_stim_run,'wt');
          if fid < 0, error('failed to open %s for writing',fname_stim_run); end
          write_stim_run(fid,times_by_run_cond{r,c},parms.stim_times_flag,ntpoints_per_run(r));
          fclose(fid);
        end
        if parms.stim_times_flag == 2
          dur = dur_by_run_cond(r,c);
          write_stim_dur(strrep(fname_stim_run,'.txt','_dur.txt'),...
            dur(~isnan(dur)),parms.forceflag);
        end
      end
    end
  catch me
    errcode = 1;
    errmsg = sprintf('stimulus timing files: %s',me.message);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = read_stim_times(fname)
  % the numbers in an events_to_timing() file; '*' or a missing file is none
  vals = [];
  if ~exist(fname,'file'), return; end
  toks = strsplit(strtrim(fileread(fname)));
  toks = toks(~cellfun(@isempty,toks) & ~strcmp(toks,'*'));
  if isempty(toks), return; end
  vals = str2double(toks);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vec = read_stim_vec(fname,ntpoints)
  % an events_to_timing() .1D mask, one value per TR after the skipped ones;
  %   its empty-condition placeholder keeps the skipped TRs, so an all-zero
  %   file of another length is taken as no events
  vec = zeros(ntpoints,1);
  vals = read_stim_times(fname);
  if isempty(vals) || ~any(vals), return; end
  if length(vals) ~= ntpoints
    error('%s has %d TRs, not %d',fname,length(vals),ntpoints);
  end
  vec = vals(:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stim_run(fid,times,stim_times_flag,ntpoints)
  % one run of a stim file: a .1D vector, or an AFNI local_times line
  if stim_times_flag == 0
    if isempty(times), times = zeros(ntpoints,1); end
    fprintf(fid,'%d\n',times);
  else
    write_stim_line(fid,times);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stim_dur(fname,durs,forceflag)
  % the one block duration of a stim file, the mean of its runs', 0 for none
  if exist(fname,'file') && ~forceflag, return; end
  if isempty(durs), dur_val = 0; else dur_val = mean(durs); end
  fid = fopen(fname,'wt');
  if fid < 0, error('failed to open %s for writing',fname); end
  fprintf(fid,'%0.4f\n',dur_val);
  fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stim_labels,conds_contrast,condnames] = prep_contrasts(task,short_flag)
  % rest has no conditions, so it returns the empty outputs rather than being
  %   an error the caller has to guard against
  % stim_labels is [] rather than {} to match the parms default: mmil_args2parms
  %   indexes into any cell value it is handed, so an empty cell errors there
  stim_labels = []; conds_contrast = []; condnames = {};
  if ~exist('short_flag','var') || isempty(short_flag), short_flag = true; end
  switch upper(task)
    case 'REST'
      return;
    case 'MID'
      [stim_labels,conds_contrast] = abcd_set_contrasts_mid(0);
      [stim_labels_out,conds_out] = abcd_set_contrasts_mid(short_flag);
    case 'SST'
      [stim_labels,conds_contrast] = abcd_set_contrasts_sst(0);
      [stim_labels_out,conds_out] = abcd_set_contrasts_sst(short_flag);
      % correct_go_vs_fixation is the correct_go condition itself, SEM
      %   included, and the condition is written too; dropped here rather
      %   than in abcd_set_contrasts_sst, which the ABCD code shares
      i_drop = find(strcmp({conds_contrast.name},'correct_go_vs_fixation'));
      conds_contrast(i_drop) = [];
      conds_out(i_drop) = [];
    case 'NBACK'
      [stim_labels,conds_contrast] = abcd_set_contrasts_nback(0);
      [stim_labels_out,conds_out] = abcd_set_contrasts_nback(short_flag);
      % the cue condition is named cue, not c; renamed here rather than in
      %   abcd_set_contrasts_nback, which the ABCD code shares
      stim_labels_out(strcmp(stim_labels_out,'c')) = {'cue'};
    otherwise
      error('unsupported task: %s; must be mid, sst, nback, or rest',task);
  end
  condnames = cat(2,stim_labels_out,{conds_out.name});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol = convert_tstat_to_sem(vol)
  sz = size(vol);
  nframes = sz(end);
  vol = reshape(vol,[],nframes);
  for i = 2:2:nframes
    vol(:,i+1) = vol(:,i) ./ (vol(:,i+1) + eps);
  end
  vol = reshape(vol,sz);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fname_deconv_roi, col_run_tag, fname_deconv_roi_norm] = deconv_roi_fname(parms)
  % the raw 3dDeconvolve output and, separately, what postnorm_deconv_roi
  %   makes of it; keeping them apart leaves the raw columns readable
  if ~isempty(parms.run_tag)
    fstem = sprintf('%s/%s_roi_deconv_3dDeconv',parms.tmpdir,parms.run_tag);
    if ~isempty(parms.col_run_tag)
      col_run_tag = parms.col_run_tag;
    else
      col_run_tag = strrep(parms.run_tag,'run','r');
    end
  else
    fstem = sprintf('%s/roi_deconv_3dDeconv',parms.tmpdir);
    col_run_tag = '';
  end
  fname_deconv_roi = sprintf('%s.mat',fstem);
  fname_deconv_roi_norm = sprintf('%s_norm.mat',fstem);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_stim_line(fid, onsets)
  if isempty(onsets)
    fprintf(fid,'*\n');
  else
    fprintf(fid,'%0.3f ',onsets);
    fprintf(fid,'\n');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function desc = tag_to_desc(tag)
  switch tag
    case '_desc-taskpred', desc = 'taskpred';
    case '_desc-residual', desc = 'residual';
    otherwise,             desc = '';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function instem_out = bids_tagged_instem(instem,tag)
  if isempty(tag)
    instem_out = instem;
  else
    instem_out = regexprep(instem,'_([^_]+)$',[tag '_$1']);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = save_mask_ico(parms,hemi)
  % coverage mask of the deconvolved runs: the vertices covered in every run,
  %   from the -ico{N}-mask-{hemi}.mgz written by cmig_extract_fmri_tseries
  errcode = 0; errmsg = [];
  bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  if ~isempty(parms.run_tag)
    run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
  else
    run_entity = '';
  end
  bids_hemi = cmig_get_bids_hemi(hemi);
  fname_out = sprintf('%s/%s%s_hemi-%s_atlas-ico5_mask.mat',...
    parms.outdir,bids_base,run_entity,bids_hemi);
  if exist(fname_out,'file') && ~parms.forceflag, return; end
  mask = [];
  for r = 1:length(parms.fnames_extract)
    [~,instem] = cmig_get_fileparts(parms.fnames_extract{r});
    fname_mask = sprintf('%s/%s-ico%d-mask-%s.mgz',parms.tmpdir,instem,parms.ico_order,hemi);
    if ~exist(fname_mask,'file')
      errcode = 1; errmsg = sprintf('ico mask file not found: %s',fname_mask);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    vals = mmil_colvec(fs_load_mgh(fname_mask)) > 0;
    if isempty(mask)
      mask = vals;
    else
      mask = mask & vals;
    end
  end
  measmat = single(mask');
  save(fname_out,'measmat','-v7');
  fprintf('%s: wrote %s\n',mfilename,fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_baseline_ico(parms, fname_mgh, fname_rbeta, hemi, nruns)
  bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  if ~isempty(parms.run_tag)
    run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
  else
    run_entity = '';
  end
  bids_hemi = cmig_get_bids_hemi(hemi);
  fname_out = sprintf('%s/%s%s_hemi-%s_atlas-ico5_baseline.mat',...
    parms.tmpdir,bids_base,run_entity,bids_hemi);
  if exist(fname_out,'file') && ~parms.forceflag, return; end

  if ~isempty(fname_rbeta) && exist(fname_rbeta,'file')
    % REML: polort-0 for run r is at sub-brick (r-1)*(detrend+1) (0-indexed)
    fname_tmp = strrep(fname_rbeta,'+orig.BRIK','_rbetatmp.mgz');
    mmil_BRIK2mgh(fname_rbeta,fname_tmp,parms.forceflag);
    [vol,~,~,volsz] = mmil_load_mgh(fname_tmp);
    delete(fname_tmp);
    nverts = prod(volsz(1:3));
    vol = reshape(vol,[nverts,volsz(4)]);
    baseline = zeros(nverts,1,'single');
    for r = 1:nruns
      sb = (r-1)*(parms.detrend+1) + 1;
      if sb <= size(vol,2)
        baseline = baseline + single(vol(:,sb));
      end
    end
  else
    % non-REML: polort-0 Coef for run r is at 1-indexed column 2+(r-1)*(detrend+1)*2
    [vol,~,~,volsz] = mmil_load_mgh(fname_mgh);
    nverts = prod(volsz(1:3));
    vol = reshape(vol,[nverts,volsz(4)]);
    baseline = zeros(nverts,1,'single');
    for r = 1:nruns
      col = 2 + (r-1)*(parms.detrend+1)*2;
      if col <= size(vol,2)
        baseline = baseline + single(vol(:,col));
      end
    end
  end
  baseline = baseline / nruns;
  measmat = baseline';
  save(fname_out,'measmat','-v7');
  fprintf('%s: wrote %s\n',mfilename,fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function errmsg = postnorm_deconv_roi(parms)
  errmsg = [];
  % the slice below is unconditional, since the readers of the output assume
  %   it; post-normalizing (baseline_flag) and censoring (outliers_flag) are
  %   what the flags gate
  [fname_deconv_roi,~,fname_out_roi] = deconv_roi_fname(parms);
  if ~exist(fname_deconv_roi,'file'), return; end
  if exist(fname_out_roi,'file') && ~parms.forceflag, return; end

  tmp = load(fname_deconv_roi);
  vol = tmp.vol;  % [nroi, nframes]
  roinames = tmp.roinames;

  % locate start of GLT (contrast) block; the offset comes from
  % cmig_deconv_col_offset, which concat_tfmri_stats also calls, so the two
  % cannot drift apart
  nruns = length(parms.fnames_data);
  col_offset = cmig_deconv_col_offset(parms,nruns);
  glt_start = 2 + col_offset;
  if glt_start > size(vol,2), return; end

  % slice down to placeholder col 1 (F-stat, read by cmig_write_tabdat_tfmri_fstat)
  % plus the GLT block only, so that downstream readers which assume
  % col_offset==0 (e.g. cmig_write_tabdat_tfmri's extract_deconv, which
  % hardcodes beta=v(:,2*c)/sem=v(:,2*c+1)) index the correct columns
  % by construction
  vol = vol(:,[1,glt_start:end]);

  % apply percent-of-baseline normalization to the (now-leading) GLT block
  if parms.baseline_flag
    bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
    if ~isempty(parms.run_tag)
      run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
    else
      run_entity = '';
    end
    fname_baseline = sprintf('%s/%s%s_roi_baseline.mat',parms.tmpdir,bids_base,run_entity);
    if exist(fname_baseline,'file')
      tmp_bl = load(fname_baseline,'baseline');
      baseline = tmp_bl.baseline;  % [nroi, 1]
      for col = 2:size(vol,2)
        vol(:,col) = 100 * vol(:,col) ./ (baseline + eps);
      end
      % a baseline of 0 means the ROI had no data, as on the ico path
      vol(baseline==0,2:end) = NaN;
      fprintf('%s: applied baseline post-normalization to roi deconv mat\n',mfilename);
    else
      errmsg = sprintf('baseline file not found, roi deconv not post-normalized: %s',...
        fname_baseline);
      fprintf('%s: WARNING: %s\n',mfilename,errmsg);
    end
  end

  % an ROI with no data in a run of the fit (all NaN, which 3dDeconvolve
  %   reads as 0) gets beta and sem NaN
  nodata = roi_nodata(parms,size(vol,1));
  if any(nodata)
    vol(nodata,2:end) = NaN;
    fprintf('%s: %d ROIs with no data in a run set to NaN\n',mfilename,nnz(nodata));
  end

  % censor outlier conditions/contrasts: NaN both beta and sem columns for
  % any condition/contrast whose sem rms across ROIs exceeds outliers_thresh
  % (cf. reference/ABCD_Analyze_taskBOLD_Exam.m's censor_outliers); save the
  % excluded set so concat_tfmri_stats can propagate the same censoring to
  % ico outputs
  if parms.outliers_flag
    fname_condnames = sprintf('%s/deconv_condnames.mat',parms.tmpdir);
    if exist(fname_condnames,'file')
      tmp_cond = load(fname_condnames,'condnames');
      condnames = tmp_cond.condnames;
      npars = length(condnames);
      excluded = false(1,npars);
      if 1+2*npars <= size(vol,2)
        for c = 1:npars
          % over the ROIs with a value, the NaN ones left out
          sem_vals = vol(:,2*c+1);
          sem_vals = sem_vals(~isnan(sem_vals));
          sem_rms = sqrt(mean(sem_vals.^2));
          if sem_rms > parms.outliers_thresh
            excluded(c) = true;
            vol(:,2*c)   = NaN;
            vol(:,2*c+1) = NaN;
            fprintf('%s: censored outlier condition %s (sem rms=%0.3f > thresh=%0.3f)\n',...
              mfilename,condnames{c},sem_rms,parms.outliers_thresh);
          end
        end
        outliers_thresh = parms.outliers_thresh;
        if ~isempty(parms.run_tag)
          fname_outliers = sprintf('%s/%s_deconv_outliers.mat',parms.tmpdir,parms.run_tag);
        else
          fname_outliers = sprintf('%s/deconv_outliers.mat',parms.tmpdir);
        end
        save(fname_outliers,'condnames','excluded','outliers_thresh','-v7.3');
      end
    end
  end

  save(fname_out_roi,'vol','roinames','-v7.3');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_baseline_roi(parms, fname_mat, fname_rbeta, roinames, nruns)
  bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  if ~isempty(parms.run_tag)
    run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
  else
    run_entity = '';
  end
  fname_out = sprintf('%s/%s%s_roi_baseline.mat',parms.tmpdir,bids_base,run_entity);
  if exist(fname_out,'file') && ~parms.forceflag, return; end

  nroi = length(roinames);
  if ~isempty(fname_rbeta) && exist(fname_rbeta,'file')
    % REML: polort-0 for run r is at sub-brick (r-1)*(detrend+1) (0-indexed)
    fname_tmp = strrep(fname_rbeta,'+orig.BRIK','_rbetatmp.mgz');
    mmil_BRIK2mgh(fname_rbeta,fname_tmp,parms.forceflag);
    [vol,~,~,volsz] = mmil_load_mgh(fname_tmp);
    delete(fname_tmp);
    vol = reshape(vol,[prod(volsz(1:3)),volsz(4)]);
    vol = vol(1:nroi,:);
    baseline = zeros(nroi,1,'single');
    for r = 1:nruns
      sb = (r-1)*(parms.detrend+1) + 1;
      if sb <= size(vol,2)
        baseline = baseline + single(vol(:,sb));
      end
    end
  else
    % non-REML: polort-0 Coef for run r is at 1-indexed column 2+(r-1)*(detrend+1)*2
    tmp = load(fname_mat,'vol');
    vol = tmp.vol;
    baseline = zeros(nroi,1,'single');
    for r = 1:nruns
      col = 2 + (r-1)*(parms.detrend+1)*2;
      if col <= size(vol,2)
        baseline = baseline + single(vol(:,col));
      end
    end
  end
  baseline = baseline / nruns;
  save(fname_out,'baseline','roinames','-v7');
  fprintf('%s: wrote %s\n',mfilename,fname_out);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nodata = roi_nodata(parms,nroi)
  % true for an ROI that is NaN in every frame of any run of the fit, from
  %   the {instem}_roi_data.mat that run_deconv_roi reads
  nodata = false(nroi,1);
  for r = 1:length(parms.fnames_extract)
    [~,instem] = cmig_get_fileparts(parms.fnames_extract{r});
    fname_roi_mat = sprintf('%s/%s_roi_data.mat',parms.tmpdir,instem);
    if ~exist(fname_roi_mat,'file'), continue; end
    tmp = load(fname_roi_mat,'data');  % [ntpoints, nroi]
    nd = mmil_colvec(all(isnan(tmp.data),1));
    n = min(nroi,length(nd));
    nodata(1:n) = nodata(1:n) | nd(1:n);
  end
return;
