# Changelog

Release notes for the CMIG BIDS pipeline, one section per version (X.Y.Z,
with its release date), newest first, added to as changes are made.

## 1.0.2 - 2026-11-25

Includes the changes of 1.0.1 (cmig_pipeline_dev only, host scripts,
2026-11-18).

- Post QC selection: the prediction layer (--predict_flag, default 1): a
  model per modality trained on the release 7.0 manual post QC, from the
  automatic post QC, raw QC __ok_, motion and Mahalanobis predictors
  (parms/autoqc_predictors.tsv, cmig_autoqc_predictors.py); rows above
  --predict_pthr are selected; the model and its scores are saved each round.
  Options --predictor_file, --predict_pthr, --predict_scanner_flag,
  --predict_max_frac_nan, --predict_n_draws, --predict_min_n (50) and
  --predict_min_nfail (5).
- cmig_autoqc_mahal.py: release 8.0 families, including a task correlation
  family per task; references fit on 10 random column draws (at most one
  column per 10 reference sessions), the distance the median over draws.
- Task fMRI motion (mr_y_qc__mot): a session mean FD over the frames the task
  model was fit to, __mot__fdsub__09_mean; and the censored and direction
  motion columns written for the tasks as for rest (57 new variables).
- cmig_check_batchdir.py: reads the per-cluster scriptlist_<cluster>.txt
  files the rounds write (it found no jobs in their batchdirs).
- run.py: group-level --manual_rawqc_dir and --manual_postqc_dir (default
  output_dir/cmig_manual_rawqc and cmig_manual_postqc, where they exist), as
  further roots of the tabulated compile; without --work_dir, SCRATCHROOT and
  MCR_CACHE_ROOT default per account (/tmp/cmig_<user>, .../mcr_cache) where
  unset or /tmp, so accounts no longer share one MATLAB Runtime cache.
- cmig_compile_bids_tabdat.py: --tables all|session|series.
- Rounds: job types (job_types in config.json: several pipelines from one
  run directory, each with its script, options, prefix, clusters, memory and
  CPUs); a retry at a higher memory limit after a job runs out of memory
  (retry_mem_factor); jobs.tsv gains job_type, oom and clustmem; ssh runs
  with no stdin; the submit log counts the jobs that got an id. SGE jobs
  with more than one CPU are submitted with -pe smp and -binding linear,
  with h_vmem per slot.
- Manual post QC selection: --tabulated_dir (manufacturer and model from
  the compiled mr_y_adm__info) replaces --qcinfo_file; no site column;
  --max_rows (a cap per modality, ranked by priority, ties by a seeded draw;
  --prev_selected_file keeps the rows still open); --filter
  new|old|any|pass|fail (default new).
- Protocol compliance: the image type is read from ImageTypeText where the
  sidecar has it (Siemens XA), else ImageType.
- Manual raw QC copy: --session_file.
- E-Prime conversion: an ASCII-only print (no crash in an ASCII locale).
- Error handling: every caught failure a module tolerates is appended to
  the status errmsg (no message overwritten or only logged); a resource
  error (out of memory, I/O, disk full, quota, stale file handle) sets
  errcode 1 even when other scans or analyses succeeded, so the job is
  rerun (with more memory after out of memory); a step that throws in
  analyze_func, and a missing or unreadable sphere map in the dMRI surface
  concat, now give errcode 1 with the message.
- Raw QC pipeline: errcode 1 when any modality fails (every modality still
  runs); mkjobs knows the per-modality raw QC status files, so a retry
  removes only the failed modality's outputs.
- Error catalog: a resource class (rules rs_memory, rs_disk, rs_io); the
  rounds raise the memory limit when a failed step's errmsg reports out of
  memory, not only when the scheduler does.
- Post QC fmap: dmri__fm__ap__pa__diff_rms is filled, from the AP dwi field
  map's first frame and the first b=0 frame of the processed PA dwi (proc
  dwi keeps no PA field map).
- Post QC selection: Python ports of the MMPS autoqc code, not yet used
  by the selection: cmig_autoqc_select.py (the multivariate path of
  autoqc_select_all.m; a model saved as JSON) and cmig_autoqc_mahal.py
  (autoqc_mahal.m with a saved reference); both match current MATLAB.
- cmig_autoqc_select.py, changed from the MATLAB: the score's histogram
  gate is replaced by a check that the Gaussian score can be formed (else
  the robust z, scaled by the SD where the MAD is 0); missing values stay
  unscored; the fallback median and MAD ignore missing values; only the
  training rows' labels fit the score; a perfectly separated fit is
  rejected with its own status.
- Raw QC (pcqc): ok and comp versions of the measures added in release 8.0,
  159 session columns: brain SNR (every series type), the auto__line
  measures and pepolar MI (dMRI and fMRI, ok and comp), and the overall FOV
  cutoff (ok max); their means divide by the scans that have the measure.
- Siemens XA rawdata: fs_load_mgh reads uint16 (MRI_USHRT) volumes and
  errors on an unknown data type; proc_func no longer crashes, without a
  status, on these sessions' field maps and bold.
- proc_func: a failed import, header read or field map header check gives
  errcode 1 and a status with the message (the import message shortened to
  its first line and mri_convert's error lines).
- Rounds: a job that ran and reported failure but left no status counts as
  failed (retried, then held), not lost.
- Tabulated compile: a file removed while the compile runs is skipped and
  listed in the report.
- Refusals: "not processed: " prefix on every run shorter than min_numTRs,
  task block durations that differ across runs, no processed fmap EPI data,
  and dMRI gradient directions rotated beyond the limit; the error rules
  class these messages as not_processed (held at once).
- Rounds: a job adopted from the queue gets a submit_date (the adoption
  time), and a done row with neither submit_date nor checked_date is
  compared with rawdata by its state_date, so reopen_changed reopens it when
  its rawdata change (213 adopted done rows were never reopened).
- Post QC fmap and func: every diff_rms (the AP-PA difference of the field
  maps, and each scan's difference from the inter-scan mean) is relative
  to the mean intensity in its mask, so scanners with different raw
  intensity scales are comparable (raw values were ~10x higher on Philips
  and Siemens XA).

## 1.0.0 - 2026-11-16

The first public release (cmig_pipeline).

- Licenses: the brain atlases (2010-2020) and the gradient nonlinearity
  maps (2006-2025, now from gradient_distortion_fields 1.0.0) are under UC
  San Diego's notice adapted for works, in place of CC-BY-NC 4.0; the paths
  in documentation/LICENSE.txt name where the image installs each part;
  the libraries of section C under UC San Diego's notice (2008-2010);
  libg2c listed; the atlas release asset carries LICENSE_atlases.txt.

- run.py (cmig_run_bids_app.py) as a BIDS app: bids_dir, output_dir,
  participant|group, --participant_label/--session_label, --module, root
  options; sessions from the step's input root for steps that do not read
  rawdata; --work_dir, --pc_rootdir, --dataset_name, --dry_flag; a FreeSurfer
  license is required for every participant-level run; the group level
  (tabdat, concat into tabulated/ and concatenated/) runs with forceflag 1 by
  default; --module cmig_compile_bids_tabdat and --module bids_validator.
- mkjobs and the rounds: image jobs run run.py in the CMIG BIDS app image
  (realpath binds, --image_binds added to them, dataset_name); the rounds
  write each device root's dataset_description.json (with realpaths in
  SourceDatasets and DatasetLinks) and atlas on the host;
  cmig_proc_bids_dwi as an image-only script; --matlab_addpath for testing
  MATLAB working copies with compiled_flag 0.
- BIDS naming: desc-preproc on proc outputs, desc-<cond> (was cond-),
  _desc-motion_timeseries.tsv, _corrmat.mat, .func.gii/.shape.gii, QC tables
  with the module as suffix and the table type in desc- (manual QC table
  names changed), n/a for empty cells, JSON sidecars for image and surface
  outputs, BIDS URIs (bids:<root>: across roots, DatasetLinks), ico
  SpatialReference fsaverage5/6 or fsaverage, GeneratedBy per package,
  CodeURL by tag, Container.Tag, .bidsignore per root.
- Versions: the image carries /opt/cmig/versions.json; every status file
  records app_version (the image tag), app_sha256 and one *_version column
  per package; pipeline_name is cmig.
- fMRI processing: the between-scan registration (bold_estimate_motion)
  starts from the header and from the brain-mask centroids and keeps the
  header unless they disagree by more than 3 mm (then the lower cost, with a
  warning); the slab-center start, which misplaced runs whose slab was moved
  more than 20 mm, is removed. The bold sidecar records
  B0UnwarpFieldMapShimMatch/Displacement/Coverage and
  BetweenScanStart/Distance, carried into the post QC series table.
- Analysis: baselines to the tmpdir, NaN for ROIs with baseline 0 or no
  data, motion regressors demeaned over the kept frames, ico5 extracted for
  rest too, one-run sessions' _deconv/_fstat under the session columns.
- dMRI: the intermediate folder (dwi/proc_dwi/, or --proc_dir) removed in
  full by cleanup;
  sidecar TORTOISEVersion, VolumeCount/ShellCount,
  ReversePEVolumeCount/ShellCount; the dwi status through the shared version
  columns.
- Job logs copied into {root}/logs/sub-X/ses-Y/ by the rounds.
- The rounds' ready_file (config.json) and cmig_bids_devices.py --ready_file:
  a participant_id,session_id list (e.g. a rawdata readiness check's) that
  limits the candidate sessions; a ready_file that is not found stops the
  round.
- Build: cmig_build_run_bids_module.py and cmig_compile_matlab take the
  repos folder from CMIG_REPOS_DIR, else the folder holding the
  cmig_pipeline tree, and the local TORTOISE folder from CMIG_TORTOISE_DIR
  or --tortoise_dir, in place of fixed site paths.
- Sidecars: the proc T1w, T2w, bold, epi and dwi sidecars have Description
  and Sources (bids:rawdata: URIs); a field map's IntendedFor lists the proc
  bold runs unwarped with it, and a bold's B0FieldSource names the pair used
  (for a fallback run, the substitute); B0UnwarpFieldMapCoverage and the
  post QC cover_prop are written with 6 decimals; TORTOISEVersion as in the
  status file; the atlas sidecars have Sources
  (bids:cmig_pipeline_atlases:...), linked in DatasetLinks to the atlas
  release asset, in place of Source and SourceDate.
- Post QC: the session mean bold is desc-meanNorm_boldref (was
  desc-meanNorm_bold).
- .bidsignore: folder patterns without a trailing slash (honored by
  bids-validator 2.4.1); per-datatype status files, QC .jpg, the proc
  motion time series, the analysis surface masks and meanNorm_boldref;
  DatasetLinks omit roots that do not exist.
- Concat: a datatype with no data is skipped (WARNING, errcode 0); each
  output .mat has a <outstem>_status.tsv (errcode and the number of
  sessions found, missing, padded and truncated) and each run a
  concat_<datatype>[_<meas>_sm<smooth>]_status.tsv, which the python check
  reads; run-axis padding is one WARNING per output file; the trimmed
  fmri_corr outputs have nruns 1.
- Compiled tables: n/a for a column a session lacks.
- Status files: app_version n/a when compiled_flag is 0.
- mkjobs: cmig_proc_bids_anat, cmig_proc_bids_func and
  cmig_analyze_bids_anat jobs pass compiled_flag.
- Release 8.0 column names (data dictionary review). Released names that
  change: the raw QC fsub/fdsub totals at 0.2 mm become fdsub__02 (48);
  the motion vol__fdsub_count, vol__fdsub__contig_count, vol__censor_count
  and vol_count become frame__fdsub__02_count,
  frame__fdsub__02__contig_count, frame__censor_count and
  frame__nodummy_count, and the task fdsub_count (0.9 mm) becomes
  frame__fdsub__09_count; the raw QC rep__ok_mean and rep__comp_mean
  become frame__ok_mean and frame__comp_mean (10); the detrended fMRI
  FWHM __ok/__comp columns take dtr__; the post QC fs_version, fs_build
  and fs__complete_indicator become fsurf_ver, fsurf_build and
  fsurf__complete_indicator; the E-Prime
  eprime__expt_version and eprime__runtime_version become eprime__expt_ver
  and eprime__runtime_ver. New-column renames: auto__line2__ to
  auto__line__fft__, _frac to _prop, rep__bN_count to frame__bN_count,
  diffrms to diff__rms, imgdist3/4 to imgdist1/2, fvox__miss_ratio to
  vox__miss_prop, inter__corr to inter__corr_score, series__<descriptor>
  to series_<descriptor>, series__vol_count to series__frame_count,
  multiple_indicator to multi_indicator, eprime__behav_indicator,
  trial__valid__resp_count; the nback condition c becomes cue. The dMRI
  raw QC fdsub columns are removed.
- FWHM: the 12 released fMRI mr_y_qc__raw__{rsfmri,tfmri__<task>}__
  {x,y,z}dir_fwhm now hold the smoothness of the data as imported, not the
  detrended value they held in release 7.0, which is now
  __dtr__{x,y,z}dir_fwhm; the raw__ and proc__ markers are dropped
  (unmarked is the data as it is at that step). Release 7.0 and 8.0 tables
  are not to be merged on these names.
- Compiled tables: new 8.0 table mr_y_qc__raw__fmri (the fMRI field map raw
  QC, was mr_y_qc__raw__fmri__fm); the tfMRI residual and task-predicted
  var and gpnet corr columns in the plain tables.
- Concat: the outputs in the release 7.0 folders (vertexwise/{smri,surf,
  tfmri,dti,isdti,rsi}, corrmat/{rsfmri,tfmri}), one vol_info.mat per
  folder, shared by fmri_corr and fmri_tseries; the run status
  concat_<datatype>[_<meas>_sm<smooth>]_status.tsv in the output root;
  smri_measures measlist default adds t1_gm, t2_gm, t1_wm and t2_wm.
- imgincl reads frame__ok_mean and frame__censor_count.

## 0.4.0 - 2026-11-11

- cmig_manual_postqc_bids.py: an artifact score of 3 fails the review (QC
  0) without asking for the pass/fail, and a reviewer's QC 1 with a score of
  3 counts as 0 in the session consensus (_score, __rev__disp_indicator),
  as for manual raw QC. Notes are required when QC is 0 with no score of 3.
- cmig_manual_rawqc_bids.py: notes are required when QC is 0 with no
  artifact scored 3 (they were required only when no artifact was scored).
- cmig_edit_manual_rawqc_abcd.py (new): writes the edited copy of the
  DAL_ABCD_QC merged_pcqcinfo csv, QC and revdisp set to the consensus of
  cmig_copy_manual_rawqc_bids_abcd.py (the severe-score override and the GE
  field-map pairs); --verify_flag 1 checks it against the original.
- cmig_edit_manual_postqc_abcd.py (new): writes edited copies of the
  release 7.0 and _desc-legacy manual post QC tables, _score and
  __rev__disp_indicator recomputed from the DAL_ABCD qcinfo_indiv csvs with
  the severe-score override; rows the csvs do not reproduce are listed, not
  edited.
- cmig_catalog_pipeline_errors.py: weekly copies of errors_log.tsv and
  errors_current.tsv in history/ (errors_log_<YYMMDD>.tsv,
  errors_current_<YYMMDD>.tsv), made when the newest copy is 7 or more days
  old; nothing in history/ is deleted.
- mmil_check_affine_trans.m: the warp score (WarpNorm; the post QC
  warp_score of cmig_postqc_bids_dwi and cmig_postqc_bids_fmap) is multiplied
  by 1/(1+exp(-(ImgDist3-0.13)/0.02)) instead of set to 0 when ImgDist3 <
  0.13; unchanged when ImgDist3 is NaN.
- cmig_check_event_pcqc.py: new session columns
  mr_y_qc__raw__{dmri,rsfmri,tfmri__<task>}__fov__cutoff__{dorsal,ventral}__ok_max,
  the max of the series FOV cutoff over the ok scans of each type (n/a when
  no ok scan has a value).
- cmig_set_imgincl_abcd.py (new): the 7 ABCD imaging inclusion flags
  (mr_y_qc__incl.tsv) from the compiled tables, with the release 8.0
  variable names; also mr_y_qc__incl_desc-criteria.tsv (one 0/1 column per
  criterion per session) and, with --report_flag 1, imgincl_report.txt
  (exclusion counts per criterion). The E-Prime override is read from the
  frozen _desc-legacy tables; the regT1 rigid and FOV cutoff criteria of
  release 7.0 are left out for now. A missing input table or column is an
  error.
- cmig_mkjobs_bids_pipeline.py: stale steps are found and rerun
  (--stale_flag, default 1), for the pipeline and the proc, fsurf, analyze
  and postqc module entries. A step the job runs is stale where an input
  step with errcode 0 ran after it started; for analyze_func, where an
  events copy in proc func/ was made after it started; for proc_anat,
  proc_func and proc_dwi, where the scans selected from the rawdata
  _scans.tsv (file name and SeriesInstanceUID) differ from those in the proc
  output sidecars (compared only where the _scans.tsv changed after the step
  started); or where an input step the job runs is stale. The job deletes
  the stale steps' outputs (analyze_func: the whole analysis func/ for a
  pipeline job without --task, else the task's analysis files, keeping the
  behavioral step's) and reruns them; each is listed in stale.tsv in the
  batch directory.
- cmig_mkjobs_bids_pipeline.py: errflag 1 deletes the outputs of each failed
  step the job runs and of the steps the job runs after it (by step, e.g.
  only the acq-func files of proc fmap/), in place of the failing module's
  directory; a failed behavioral step, and scripts without the step table,
  as before.
- cmig_select_bids_func_runs.py: with exclude_noncompliant_flag 1, the
  behavioral step is out of date where the compliant 0 bolds of a task in the
  rawdata _scans.tsv differ from the runs its _desc-select_runs.tsv notes as
  not compliant.
- cmig_fsurf_bids.py: a skip on a complete recon keeps an existing status
  with errcode 0 (procdate unchanged); otherwise it writes a new errcode 0
  status.
- cmig_proc_bids_func.m: a run that only copies new events files leaves the
  status file unchanged (procflag 0), so the post QC steps that do not read
  events are not made stale.
- cmig_write_status.py: inputs_newer returns (input path, its time, the
  output's time) or None, and takes ok_only (inputs whose errcode is not 0
  are skipped).
- cmig_run_bids_rounds.py: a done session returns to the candidates (state
  changed) where its rawdata _scans.tsv or an events file changed after its
  last job was submitted or it was last found done (jobs.tsv checked_date);
  a job that ends with steps stale gets state stale (jobs.tsv stale count,
  held after max_stale, default 2) with a WARNING in the round log, in place
  of lost; rounds.tsv gains the columns changed, stale and reopened.
- cmig_check_batchdir.py: with --raw_rootdir and the output roots
  (--proc_rootdir, --fsurf_rootdir, --analysis_rootdir, --postqc_rootdir),
  the batch's sessions are checked for steps left stale (as a job that runs
  only some steps leaves the steps after them) and listed as warnings in the
  summary.
- cmig_select_manual_postqc_bids.py: the rawqc layer, on by default
  (--rawqc_flag 1) with --rawqc_series_file, the compiled manual raw QC
  series table (mr_y_qc__raw__man__series.tsv): a session and modality is
  selected where a scan its proc outputs used (by the SeriesInstanceUID of
  the proc T1w, T2w, bold and dwi sidecars) passed manual raw QC with a QU
  score of 2 or more (rawqc_borderline) or with the reviewers split
  (rawqc_disp_pass); T1w feeds the fsurf review, T2w t2, bold fmri, dwi
  dmri. The scans selected are listed in
  manual_postqc_selected_<date>_rawqc_detail.tsv.

## 0.3.0 - 2026-11-06

- cmig_catalog_pipeline_errors.py (new): a cumulative catalog of processing
  errors by session and step (errors_log.tsv, every occurrence;
  errors_current.tsv, the latest per step with whether a later run
  succeeded, and hand notes from errors_notes.tsv), from mkjobs batchdirs or
  listed sessions; each error classed (not_processed, missing_input,
  failure, file_read, infrastructure, warning, unknown) by the rules in
  parms/errors_rules.tsv (new) or a file given with --rules.
- cmig_mkjobs_bids_pipeline.py: status_rows() lists every status file a
  session is expected to have for a script, with its row.
- cmig_run_bids_rounds.py: a session whose failing steps are all refusals
  (class not_processed) is held at once instead of retried with errflag 1
  and 2; config error_rules names the rules file.
- cmig_pipeline_bids_rawqc.py: the overall status is
  {pid}_{sid}_desc-rawqc_status.tsv, the file cmig_mkjobs_bids_pipeline.py
  checks, written by cmig_write_status with the version columns, instead of
  {pid}_{sid}_rawqc_status.json, so mkjobs no longer sees every rawqc
  session as pending.
- A MATLAB module that declines a session for its inputs now starts its
  status errmsg with "not processed: ", as cmig_proc_bids_dwi.py already
  did: no usable T1w/T2w (cmig_proc_bids_anat, cmig_analyze_bids_anat), no
  data found (cmig_rawqc_bids), no usable func scans or, with
  B0unwarp_flag 2, no usable field maps (cmig_proc_bids_func), no usable
  scans to analyze and a volume-only recon (cmig_analyze_bids_func,
  cmig_analyze_bids_dwi).
- cmig_manual_rawqc_bids.py (new): manual raw QC review of BIDS data with
  the cmig_rawqc_bids images: a 0-3 severity per artifact (per-modality
  lists) and a pass/fail per scan at the terminal; s at any prompt skips
  the scan; each scan is confirmed, repeated or skipped, and saved when
  confirmed, into the same per-session files as
  cmig_copy_manual_rawqc_bids_abcd.py; review and save locks for several
  reviewers.
- cmig_copy_manual_rawqc_bids_abcd.py: the severe-score override (a
  reviewer's pass with an artifact scored severe counts as a fail) applies
  to the consensus only; the per-reviewer rows keep the QC as given. The
  desc-manualrawqc_scans.tsv and desc-manualrawqc_qc.tsv headers use
  release-style names (tables mr_y_qc__raw__man__series and
  mr_y_qc__raw__man__series__review); files with the old names are still
  read, also by cmig_check_event_pcqc.py.
- cmig_compile_bids_tabdat.py: --series_outdir (default --outdir) for the
  per-series tables; a table is named from its first column up to
  __series; the per-reviewer manual raw and post QC files are compiled into
  mr_y_qc__raw__man__series__review and
  mr_y_qc__post__man__<modality>__review in the series outdir, and
  per-reviewer files with the old column names are skipped and reported.
- E-Prime events sidecars (task-<task>_events.json, written by
  cmig_abcd_eprime_to_events.py --rootdir): shorter descriptions of onset,
  task_name, start_delay, start_delay_diff, init_delay_diff, start_time,
  timing_ok, response_early, response_antic, feedback_window,
  response_time_stop_adjusted, response_time_ssd, the REC block_type and
  the SST early_stop trial type; response_accuracy is described as whether
  the response to the stimulus was correct.
- E-Prime events columns renamed in the _events.tsv files and their
  task-<task>_events.json sidecars (cmig_abcd_eprime_to_events.py): run to
  run_num, experiment_name to task_name, experiment_version to
  task_version, acc to response_correct, early_response to response_early,
  antic_resp to response_antic, stim_acc to response_accuracy,
  correct_go_response to response_correct_go, ssd_dur to duration_ssd,
  go_rt_adjusted to response_time_go_adjusted, stop_rt_adjusted to
  response_time_stop_adjusted, ssd_rt to response_time_ssd. The readers
  (cmig_abcd_events_to_behavioral.py, cmig_abcd_events_to_timing.py,
  cmig_select_bids_func_runs.py, cmig_pc_bids_abcd.py) take the new names
  only, so events files written before this change must be regenerated.
- cmig_select_manual_postqc_bids.py (new): selects the sessions for manual
  post-processing QC, per session and modality (t2, fmri, dmri, fsurf), as
  the union of layers: previous failures (release 7.0 and legacy manual post
  QC tables) and a random sample per modality stratified by event and
  manufacturer; rows without review inputs, or not in --id_events_file, are
  written to a _dropped file. The borderline raw QC and multivariate
  prediction layers are not yet implemented. --outdir, --review_rootdir,
  --prev_dir and --id_events_file are required (no ABCD path defaults).
- cmig_manual_postqc_bids.py: --selected_file may carry a modality column
  (each session reviewed for its listed modalities only) and a
  pipeline_rootdir column (the cmig_pipeline root of the session's inputs);
  --random_order (default 1) and --random_seed; --output_rootdir is
  required. The per-reviewer rows go to one file per modality,
  {pid}_{sid}_desc-manualpostqc<modality>_qc.tsv, with release-style column
  names (table mr_y_qc__post__man__<modality>__review); an old single
  desc-manualpostqc_qc.tsv is read and converted at the next save.
- dwi brain mask: cmig_proc_bids_dwi no longer writes a desc-brain_mask (the
  T1w mri_synthstrip mask resampled onto the dwi grid) into cmig_proc/dwi;
  cmig_check_bids_dwi_outputs.py no longer requires it. The shared dwi mask
  is the one cmig_analyze_bids_dwi.m writes beside the DTI and RSI
  coefficients (_mask.nii.gz), now always mri_synthstrip of the average
  b=0 of the fitted volume, never the processing's mask. The fits are not
  masked with it.
- cmig_classify_abcd_series.py (new): classifies an ABCD imaging series
  from the DICOM header of one file (series type, phase encoding
  direction, task, NORM and single-band variants, a description fallback),
  ported from the dispatch blocks of MMPS abcd_pc_GE.m, abcd_pc_Philips.m
  and abcd_pc_Siemens.m, for sourcedata generation; get_metadata_value
  reads values the enhanced (XA) object moves into its functional groups.
  Needs pydicom.
- cmig_calc_fmri_roi_corr.m: the ROI correlations and variances (the
  per-session corrmat and var tables and the concatenated fmri_corr files)
  leave out the ROIs whose mean is a WM or CSF regressor, those in wm_names
  (if wm_flag) and csf_names (if csf_flag), by default left and right
  cerebral white matter and lateral ventricle: with the mean of left and
  right in the design, their residuals are exact negatives of each other;
  cmig_prep_fmri_roi_tseries.m records them as regressor_roinames. The
  tfMRI betas keep them.
- cmig_compile_bids_tabdat.py: --id_events_file, an id_events release csv;
  only the visits it lists (participant_id and session_id) are compiled, and
  the report lists the processed visits left out.
- cmig_concat_bids_sessions.m: analysis_rootdir may be several roots, whose
  sessions are combined in participant_id, session_id order (a session
  under more than one root is an error); option id_events_file, as above;
  cmig_mkjobs_bids_concat.py and cmig_run_bids_module.py pass several
  --analysis_rootdir values and --id_events_file.
- cmig_cluster_utils.py: job_usage, the memory, CPU and disk I/O of
  finished jobs from the scheduler accounting (sacct on slurm, qacct on
  SGE): elapsed and CPU hours, cores, memory (resident MaxRSS on slurm,
  virtual maxvmem on SGE), GB written (both) and read (slurm only).
- cmig_check_batchdir.py: reports these per job (job ID and host from each
  .out header), in its jobs tsv, a usage section of the summary and
  histograms by cluster; --usageflag (default 1).
- cmig_run_bids_rounds.py: records them for each finished job in jobs.tsv,
  and per round and cluster in rounds.tsv (jobs finished, median memory and
  cores, total GB written and read).
- cmig_run_bids_module.m: runs the module with maxNumCompThreads set to
  CMIG_MATLAB_NTHREADS, when set, and prints the thread count;
  cmig_sing_utils.py sets CMIG_MATLAB_NTHREADS=1 with --reproducible_flag 1
  (the default), so MATLAB's results no longer depend on the cores of the
  host (on the SGE clusters a job used all of a node's cores: the fMRI
  analysis ROI time series and correlations, meanNorm bold and
  inter__diff_rms differed by up to about 1e-4); cmig_run_bids_module.py
  runs compiled_flag 0 through cmig_run_bids_module.m as well.
- mmil_check_affine_trans.m: its registration-to-T1 measures (e.g. the
  fmap rigid_score; used by cmig_postqc_bids_anat, cmig_postqc_bids_dwi,
  cmig_postqc_bids_fmap and cmig_analyze_bids_dwi), from imregtform after an
  8-bit histeq, changed with the thread count by up to 9e-03; the fixed
  thread count above covers them (a one-thread setting added in this
  release was removed again).
- cmig_pipeline_bids.py: --dwi_n_cpus (default 8), the thread count of the
  dwi processing tools (it followed the CPUs the job may use, all of a node
  on the SGE clusters); cmig_mkjobs_bids_pipeline.py passes it through.
- cmig_pipeline_bids.py: runs the dMRI steps, with --dwi_flag (default 1):
  cmig_proc_bids_dwi after cmig_proc_bids_func, for a session with a
  processed T1w and usable dwi data in rawdata; cmig_analyze_bids_dwi with a
  complete recon; cmig_postqc_bids_dwi and cmig_postqc_bids_fmap for the
  acq-dwi field maps. A session whose dwi cannot be processed returns
  errcode 1, as for any failed step. The steps need the image's dwi tools.
- cmig_mkjobs_bids_pipeline.py: --dwi_flag, passed to the pipeline; a
  session's completeness includes the dwi status files when its dwi steps
  run.
- cmig_analyze_bids_func.m: the completeness check expects the fstat and
  deconv tsvs for aseg, dsk and dst only, the atlases they are written for;
  it expected gp too, so each resubmit redid the mid, nback and sst
  analyses.
- README.txt and LICENSE.txt: the dMRI processing's authors and outside
  code (Diliana Pecheva in the Authors; cmig_proc_dwi.bash, based on M.
  Okan Irfanoglu's script, a United States Government work; the dwi python
  files' notice; cmig_synb0_disco.sh, based on Synb0-DISCO's pipeline.sh).
- Dockerfile: SynthStrip's infant model (models/synthstrip.infant.1.pt,
  by sha256) added to the image's FreeSurfer, which cmig_proc_dwi.bash needs
  for HBCD sessions.
- Dockerfile: TORTOISE pinned to release TORTOISEV4.3.0 (asset
  TORTOISEV4.3.0.tar.gz, by sha256), replacing the dwi code's bundled
  TORTOISE programs; four of the programs the dwi step runs
  (ApplyTransformationToScalar2, DisplacementFieldOperations,
  SetFieldToIsocenter, ZeroFieldDims) are in its bin/DRTAMAS/, which is
  added to the PATH.
- cmig_postqc_bids_fmap.m: calc_fmap_fwhm uses a working folder per field
  map group (fmap_fwhm_<label>); cmig_3dFWHM keeps an existing converted
  mask, so the func PA and dwi AP and PA FWHM were measured in the func AP
  mask.
- cmig_detect_brain_cutoff.m: passes fspath to mmil_brain_cutoff; with
  forceflag 1, cmig_postqc_bids_fmap stopped at the brain cutoff ("file
  /mri/aseg.mgz not found"), since mmil_brain_cutoff rebuilt the aseg brain
  mask without the FreeSurfer path.
- cmig_proc_bids_dwi.py: dMRI processing for one session (the dwi code
  from abcd_preproc: bash/cmig_proc_dwi.bash, which runs Synb0-DISCO
  through bash/cmig_synb0_disco.sh, and TORTOISE, with the cmig_*dwi* and
  related python helpers). Options --participant_id, --session_id,
  --raw_rootdir, --output_rootdir, --proc_anat_rootdir (the processed T1w;
  default output_rootdir), --fs_license, --n_cpus, --output_res,
  --output_size, --output_orientation, --max_rotation, --cleanup_flag and
  --forceflag; writes dwi/{pid}_{sid}_desc-proc_status.tsv, with errcode 1
  and the reason for a session it cannot process. The tools come from the
  environment (FSLDIR, ANTSPATH, FREESURFER_HOME, SYNB0_DISCO_DIR,
  CMIG_PIPELINE_ATLASDIR, GRADIENT_DISTORTION_FIELDS_DIR); HBCD sessions
  need FreeSurfer's models/synthstrip.infant.1.pt.
- cmig_build_run_bids_module.py: the image gets bash/.
- cmig_concat_fmri_ico_tseries.m (renamed from cmig_concat_tfmri_ico_tseries.m,
  since it is used for task and rest): each run's ico time series is NaN
  at vertices outside that run's coverage mask (-ico{N}-mask-{hemi}.mgz),
  where it was 0. The prep, residual and taskpred series and their
  correlations are NaN there; before, the taskpred correlations at those
  vertices were values computed from rounding noise. A missing mask file is
  an error. cmig_concat_tfmri_roi_tseries.m is renamed
  cmig_concat_fmri_roi_tseries.m.
- cmig_abcd_eprime_to_events.py: the events files no longer have a task
  column, which repeated the filename's task entity;
  cmig_abcd_events_to_behavioral.py and cmig_abcd_events_to_timing.py take
  the task from the filename (task-<label>), also for older files that
  have the column. The events sidecars give Levels for switch_flag,
  timing_ok, acc, early_response, antic_resp, feedback_window, stim_acc and
  correct_stop, Units "s" for onset and duration, and describe
  experiment_name as the task program run, named for its version and
  variant, with an example for each task.
- cmig_abcd_eprime_to_events.py: write_root_sidecars() (and the command
  line --rootdir DIR) writes task-mid_events.json, task-nback_events.json,
  task-sst_events.json and task-rec_events.json, the column descriptions of
  the E-Prime events files, at the rawdata root; eprime_to_events() no longer
  writes a per-session _events.json. trial_type Levels list every trial
  type; NBACK block_type Levels are 0_back/2_back, the values written.
- cmig_reg_T1w_to_atlas.m: AtlasRegistrationCorrelation is computed in
  double precision; in single it differed in the 7th digit with the number
  of MATLAB threads (e.g. between hosts with different CPU counts).
- mmil_synthstrip_brainmask.m: option nthreads (default 8), passed to
  mri_synthstrip as -t; with PyTorch's default (all CPUs) the mask could
  differ by a voxel at the border between hosts with different CPU counts.
- docker/Dockerfile: the FSL, AFNI, MATLAB Runtime and FreeSurfer stages
  install findutils, which prune_keep.sh needs (find, xargs).
- setup/SetUpCMIG.csh with no version argument sets up the highest version
  installed, so one copy (e.g. in /usr/pubsw/bin) works across versions.
- docker/files/keep_freesurfer.txt keeps the 66 FreeSurfer programs named
  in the FreeSurfer scripts the trace ran (e.g. AntsN4BiasFieldCorrectionFs,
  run by mri_nu_correct.mni), whose runs the trace did not record; without
  them recon-all stopped in the pruned image.
- python/cmig_list_bids_sessions.py (new) lists the sessions in a BIDS
  rawdata tree, with the modality directories each holds and its
  _scans.tsv mtime, into a TSV (sessions.tsv), merged with an existing one
  (sessions no longer found kept with present 0).
- --reproducible_flag (default 1) in cmig_run_bids_module.py and the two
  pipelines: the MATLAB modules run with OPENBLAS_CORETYPE=HASWELL,
  OPENBLAS_NUM_THREADS=1 and MKL_CBWR=AVX2 (each where not already set), so
  proc func gives the same results on every CPU type (EPYC 7763 and 7401P
  tested); 0 keeps the inherited environment, for a CPU without AVX2. Only
  the module's own call gets them, not the pipeline's environment.
- cmig_mkjobs_bids_pipeline.py: public check_dirs, write_root_files
  (dataset_description.json and the atlas, once per set of roots),
  make_opts, load_delay_flags, check_session and make_session_job (one
  session's status check and job script), used by its batch loop; job
  scripts unchanged; it passes --reproducible_flag when given.
- New, for running a pipeline over many sessions across storage devices
  and clusters: cmig_bids_devices.py (devices.tsv, assignments.tsv, df,
  projected use, device assignment), cmig_cluster_utils.py (query, submit
  and cancel jobs on slurm and SGE head nodes over ssh) and
  cmig_run_bids_rounds.py (one round: settle jobs, top each cluster up,
  assign devices, write and submit job scripts).
- cmig_manual_postqc_bids.py: a save lock around each review's write and
  a review lock on a visit and modality while it is reviewed, so two
  reviewers at once lose no rows and max_nrev holds; options --lock_hours
  (default 4) and --override_lock.
- docker/Dockerfile downloads the compiled app and the atlases from a
  cmig_pipeline release, and the gradient maps from the
  gradient_distortion_fields tag, each pinned by tag and sha256 in its
  ARGs (empty until the first release: the piece then comes from the
  build context's local/); the image labels record the source of each.
  cmig_build_run_bids_module.py reads the pins (--release_tag,
  --atlas_tag, --gdf_tag override them) and, with --assets_dir, writes the
  two release assets and prints their sha256.
- cmig_pipeline_bids.py --analyze_anat_flag (default 1): 0 skips
  cmig_analyze_bids_anat; cmig_mkjobs_bids_pipeline.py passes it when
  given.
- setup/versions.txt: cmig_pipeline_atlases 0.1.0 (0.0.0 plus NIfTI
  versions of T1T2 atlasT1, atlasT2 and maskT2, for the dwi step's atlas
  registration). cmig_build_run_bids_module.py takes the atlases for the
  image and the atlas release asset from that install, in place of
  atlases.2022.10.17.
- docker/Dockerfile: every stage FROM BASE_IMAGE, rockylinux:8.8 pinned by
  digest. cmig_build_run_bids_module.py: --pipeline_dir defaults to the
  tree it runs from and --docker_srcdir to its docker/; the build context
  is assembled from empty, without python caches; a docker build whose app
  is pinned by tag does not compile or write the help json.
- docker/Dockerfile for the dwi step: Python 3.12 (rpm-pinned) in place of
  python36, with numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, torch 2.14.1+cpu
  and the rest pinned with hashes (docker/files/requirements.txt), each
  package's license copied into the image; base image rockylinux 8.10 by
  digest; new stages for ANTs 2.6.5 (7 programs), Convert3D 1.1.0
  (c3d_affine_tool), Synb0-DISCO (commit d813f94 with
  docker/files/synb0_disco.patch: CPU inference and two pipeline edits) and
  TORTOISE (pinned by TORTOISE_TAG/ASSET/SHA256 once released, until then
  from a local bundle pruned to docker/files/keep_tortoise.txt); 12 more FSL
  programs kept; fsl_abspath replaced by a shell script so FSL's python
  need not be kept; the final stage built FROM a base stage (Rocky,
  system packages, Python and its packages); LICENSE.txt and licenses/ for
  the new tools.
- cmig_build_run_bids_module.py: --tortoise_tag, --tortoise_asset,
  --tortoise_sha256 and --tortoise_dir; --python_sing_flag and
  --python_image build the base stage alone and save it as a python
  singularity image ({appdir}/images/{python_image}-{appver}.sif).
- csh/cmig_python3 (new): runs python3 from that image with /space, /home
  and /usr/pubsw bound and the environment passed (CMIG_PYTHON_IMAGE,
  CMIG_PYTHON_SING).
- cmig_mkjobs_bids_pipeline.py --python_cmd (default python3): the
  command the job scripts run the python with, e.g. cmig_python3.
- docker/documentation/LICENSE.txt sections A and B: the code written at
  UC San Diego (2006 to April 2025) under The Regents of the University of
  California's notice (educational, research and non-profit use); the code
  written since, at JCVI, under the PolyForm Noncommercial License 1.0.0
  (licenses/PolyForm-Noncommercial-1.0.0.md), in place of the MIT notice.
  README.txt says the same and names both institutions in its
  Acknowledgments.

## 0.2.0 - 2026-10-27

- mmil_aseg_analysis and mmil_aseg_roi take fname_excolorlut (the extra
  color lookup table), passed to mmil_colorlut; empty keeps its default.

## 0.1.0 - 2026-10-27

- setup/SetUpCMIG.csh sets up the environment for an installed version
  (CMIG variables, SURFACE_MORPHOMETRY_TOOLBOX_DIR, PATH, PYTHONPATH), by
  itself or after SetUpMMPS.csh; setup/startup_CMIG.m adds the toolbox, the
  gradient maps and the pipeline's matlab/ to the MATLAB path.

## 0.0.0 - 2026-10-27

- The pipeline no longer needs MMPS: it finds its files through the
  environment variables CMIG_PIPELINE_DIR, CMIG_PIPELINE_VERSION,
  CMIG_PIPELINE_ATLASDIR and GRADIENT_DISTORTION_FIELDS_DIR, and stops with
  an error naming any that is not set.
- Status files and dataset_description.json record software_name
  cmig_pipeline and the version from CMIG_PIPELINE_VERSION.
- compiled_flag 2 (the interim singularity image) is removed; the choices
  are 0 (MATLAB code) and 1 (the compiled app).
- The compiled app is built by cmig_compile_matlab.m from cmig_pipeline,
  the surface_morphometry_toolbox and gradient_distortion_fields only.
- The image sets CMIG_DOCKER (was MMPS_DOCKER) and the CMIG variables, and
  installs cmig_pipeline, its atlases and the gradient maps under
  /usr/pubsw/packages.
- The build script's --prune_flag (default 1) prunes FSL, AFNI, the MATLAB
  Runtime and FreeSurfer in the image to keep-lists (docker/files).
- The extra colour lookup table is CMIG_ColorLUT.txt; MMIL_ColorLUT.txt is
  removed.
- Gradient distortion maps are found by gradient_distortion_files.m
  (gradient_distortion_fields); mmil_grad_files.m is removed.
