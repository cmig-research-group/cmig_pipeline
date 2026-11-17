function cmig_concat_fmri_roi_tseries(fnames_data,fnames_motion,fnames_events,indir,varargin)
%function cmig_concat_fmri_roi_tseries(fnames_data,fnames_motion,fnames_events,indir,[options])
%
% Purpose: concatenate surface-based and volumetric ROI
%          task or rest fMRI time series data
%          from output of cmig_extract_fmri_tseries
%
% Required Input:
%   fnames_data: cell array of processed data files
%   fnames_motion: cell array of _desc-motion_timeseries.tsv file names
%   fnames_events: cell array of events tsv file names
%   indir: input directory containing output of cmig_extract_fmri_tseries
%
% Optional Input ('key',value):
%   'fname_out': output file name
%     {default = 'concat_fmri_roi_tseries.mat'}
%   'skipTRs': number of dummy TRs already removed from the time series data,
%     used to align the motion tsv and the events onsets to it
%     {default = 10}
%   'TR': scan repetition time (seconds); required for onset/duration conversion
%     {default = []}
%   'fnames_design': cell array of per-run HRF-convolved design matrix files
%     (xmat.1D) from prep_design_per_run; if provided, stored as datamat_design
%     {default = {}}
%   'nmotion': number of motion regressor columns in design matrix
%     {default = 12}
%   'stim_labels': cell array of stimulus condition labels (from prep_contrasts)
%     if provided, used to order onset/duration columns; saved in output
%     {default = {}}
%   'nreg_per_stim': number of regressors per stimulus condition in design matrix
%     used with stimregs to compute number of conditions
%     {default = 2}
%   'detrend': polynomial detrending order used in 3dDeconvolve
%     used to compute nbase = detrend+1 for stimregs
%     {default = 3}
%   'censor_thresh': framewise displacement threshold for motion censoring (mm)
%     {default = 0.2}
%   'contig_thresh': minimum number of contiguous uncensored time points
%     {default = 5}
%   'censor_thresh_deconv': framewise displacement threshold for censvec_deconv,
%     which the task refit in cmig_task_predicted uses instead of censvec
%     {default = 0.9}
%   'contig_thresh_deconv': contig_thresh for censvec_deconv. the
%     deconvolution uses 0 where the correlation uses 5, so this is
%     deliberately not contig_thresh
%     {default = 0}
%   'min_ntpoints': minimum number of time points per run
%     {default = 100}
%   'max_ntpoints': minimum size of time point dimension in output matrix
%     {default = 200}
%   'roitypes': cell array of atlas types to include in output
%     allowed: 'gp', 'dsk', 'dst', 'aseg'
%     if empty, include all ROIs (backward compatible)
%     {default = {}}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  05/07/2026 by Don Hagler
% Prev Mod: 05/14/2026 by Don Hagler
% Prev Mod: 05/20/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 09/07/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/16/2026 by Don Hagler
% Prev Mod: 09/17/2026 by Don Hagler
% Prev Mod: 10/16/2026 by Don Hagler
% Prev Mod: 11/02/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on:
%   cmig_concat_rsfmri_roi_tseries
%     Created:  09/08/25 by Don Hagler
%     Last Mod: 09/09/25 by Don Hagler
%   ABCD_concat_roi_tfMRI_tseries
%     Created:  02/25/22 by Don Hagler
%     Last Mod: 12/29/25 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,4), return; end

parms = mmil_args2parms(varargin,{...
  'fname_out','concat_fmri_roi_tseries.mat',[],...
  'skipTRs',10,[],...
  'TR',[],[],...
  'fnames_design',{},[],...
  'nmotion',12,[],...
  'stim_labels',{},[],...
  'nreg_per_stim',2,[],...
  'detrend',3,[],...
  'censor_thresh',0.2,[],...
  'contig_thresh',5,[],...
  'censor_thresh_deconv',0.9,[],...
  'contig_thresh_deconv',0,[],...
  'min_ntpoints',100,[],...
  'max_ntpoints',200,[],...
  'roitypes',{},[],...
  'forceflag',false,[],...
});

ndirs = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input directory
if ~exist(indir,'dir'), error('indir %s not found',indir); end

% create output directory
outdir = fileparts(parms.fname_out);
if ~isempty(outdir), mmil_mkdir(outdir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set nruns_in based on number of data files
nruns_in = length(fnames_data);

if ~exist(parms.fname_out,'file') || parms.forceflag

  % load all fMRI ROI time series, motion, and events data
  array_tsdata       = {};
  array_motion       = {};
  array_onset        = {};
  array_duration     = {};
  array_stim_types   = {};
  array_design       = {};
  array_ntpoints_all = zeros(1,nruns_in);  % per-run ntpoints (including skipped runs)
  seri_per_run       = zeros(1,nruns_in);  % seri index for each run (0 if skipped)
  nframes_run        = zeros(1,nruns_in);  % per-run frame count before skipTRs removal
  stim_types = {};  % global stim types (union across all runs)
  roinames = [];
  nroi = [];
  nmot = 0;
  nreg = 0;
  diri = 1;
  seri = 0;
  for r = 1:nruns_in
    % get instem from fnames_data
    fname_data = fnames_data{r};
    [~,instem] = cmig_get_fileparts(fname_data);
    % set ROI timeseries data file name
    fname_tsdata = sprintf('%s/%s_roi_data.mat',indir,instem);
    if ~exist(fname_tsdata,'file')
      error('time series data file %s not found',fname_tsdata);
    end
    % get fname_motion from fnames_motion
    fname_motion = fnames_motion{r};
    if ~exist(fname_motion,'file')
      error('motion file %s not found',fname_motion);
    end
    % load ROI time series data
    tmp_tsdata = load(fname_tsdata);
    nframes_run(r) = size(tmp_tsdata.data,1) + parms.skipTRs;
    ntpoints_run = size(tmp_tsdata.data,1);
    array_ntpoints_all(r) = ntpoints_run;
    % exclude runs that are too short
    if ntpoints_run < parms.min_ntpoints
      fprintf('%s: WARNING: skipping run %d because ntpoints = %d\n',mfilename,r,ntpoints_run);
      continue;
    end
    parms.max_ntpoints = max(parms.max_ntpoints,ntpoints_run);
    roinames = tmp_tsdata.roinames;
    nroi = length(roinames);
    % load motion
    tmp_motion = mmil_load_motion_tsv(fname_motion);
    tmp_motion = tmp_motion(parms.skipTRs+1:end,:);
    [ntpoints_motion,nmot] = size(tmp_motion);
    if ntpoints_motion < ntpoints_run
      tmp_motion = cat(1,tmp_motion,zeros(ntpoints_run-ntpoints_motion,nmot));
    end
    % the data may stop short of the scan, capped at a frame count
    tmp_motion = tmp_motion(1:ntpoints_run,:);
    % load events and convert to per-TR onset/duration matrices
    fname_events_run = [];
    if ~isempty(fnames_events) && r <= length(fnames_events)
      fname_events_run = fnames_events{r};
    end
    [tmp_onset,tmp_duration,run_stim_types] = cmig_load_events_tsv(fname_events_run,...
      ntpoints_run,parms.TR,parms.stim_labels,parms.skipTRs*parms.TR);
    % update global stim_types with any new types found in this run
    new_types = setdiff(run_stim_types,stim_types,'stable');
    stim_types = [stim_types, new_types];
    % store in arrays (seri indexes only valid runs)
    seri = seri + 1;
    seri_per_run(r) = seri;
    array_tsdata{seri}     = tmp_tsdata.data;
    array_motion{seri}     = tmp_motion;
    array_onset{seri}      = tmp_onset;
    array_duration{seri}   = tmp_duration;
    array_stim_types{seri} = run_stim_types;
    % load per-run design matrix if available
    if ~isempty(parms.fnames_design) && r <= length(parms.fnames_design) && ...
        ~isempty(parms.fnames_design{r}) && exist(parms.fnames_design{r},'file')
      X_run = cell2mat(mmil_readtext(parms.fnames_design{r},' ','#'));
      array_design{seri} = single(X_run);
      nreg = size(X_run,2);
    end
  end
  ntpoints = parms.max_ntpoints;
  nruns = seri;
  if ~nruns
    fprintf('%s: WARNING: no valid runs\n',mfilename);
    return;
  end

  % reduce nframes_run to valid runs only, in seri order; nframes is the
  % session total number of frames before skipTRs removal, summed across
  % valid runs (NOT vol_count/ntpoints + skipTRs, which would undercount
  % for nruns > 1 since each run is trimmed individually)
  nframes_run = nframes_run(seri_per_run > 0);
  nframes = sum(nframes_run);
  TR = parms.TR;

  % resolve stim_labels: use provided labels if available, else events-derived types
  if ~isempty(parms.stim_labels)
    stim_labels = parms.stim_labels;
  else
    stim_labels = stim_types;
  end
  nstim = length(stim_labels);

  % compute stimregs from per-run design matrix width
  nbase = parms.detrend + 1;
  nmotion = parms.nmotion;
  stimregs = [];
  if nreg > 0 && ~isempty(stim_labels)
    nstim_cols = nreg - nbase - nmotion;
    if nstim_cols > 0
      stimregs = nbase+1 : nbase+nstim_cols;
    end
  end

  % initialize output matrices
  censvec          = false(nruns,ndirs,ntpoints);
  % the refit in cmig_task_predicted belongs on the deconvolution's frames,
  %   not the correlation's, so a second vector is built here where the
  %   motion and the censoring call already are
  censvec_deconv   = false(nruns,ndirs,ntpoints);
  datamat_tsdata   = NaN([nruns,ndirs,ntpoints,nroi],'single');
  datamat_motion   = NaN([nruns,ndirs,ntpoints,nmot],'single');
  datamat_onset    = NaN([nruns,ndirs,ntpoints,nstim],'single');
  datamat_duration = NaN([nruns,ndirs,ntpoints,nstim],'single');
  if nreg > 0
    datamat_design = NaN([nruns,ndirs,ntpoints,nreg],'single');
  else
    datamat_design = [];
  end

  % concatenate task fMRI ROI time series data
  for seri = 1:nruns
    tmp_data       = array_tsdata{seri};
    tmp_motion     = array_motion{seri};
    tmp_onset      = array_onset{seri};
    tmp_duration   = array_duration{seri};
    run_stim_types = array_stim_types{seri};
    % calculate censor vector from motion
    tmp_censor = mmil_censor_motion(tmp_motion,parms.censor_thresh,parms.contig_thresh);
    % available time points this run
    idx_time = 1:size(tmp_data,1);
    % add data to output matrices
    datamat_tsdata(seri,diri,idx_time,:) = tmp_data;
    datamat_motion(seri,diri,idx_time,:) = tmp_motion;
    censvec(seri,diri,idx_time)          = tmp_censor;
    censvec_deconv(seri,diri,idx_time)   = mmil_censor_motion(tmp_motion,...
      parms.censor_thresh_deconv,parms.contig_thresh_deconv);
    % align per-run onset/duration to stim_labels order
    for c = 1:length(run_stim_types)
      idx_stim = find(strcmp(run_stim_types{c},stim_labels));
      if ~isempty(idx_stim)
        datamat_onset(seri,diri,idx_time,idx_stim)    = tmp_onset(:,c);
        datamat_duration(seri,diri,idx_time,idx_stim) = tmp_duration(:,c);
      end
    end
    % design matrix slice for this run
    if nreg > 0 && seri <= length(array_design) && ~isempty(array_design{seri})
      nt_design = min(size(array_design{seri},1),ntpoints);
      datamat_design(seri,diri,1:nt_design,:) = array_design{seri}(1:nt_design,:);
    end
  end

  % filter roinames by roitypes if specified
  if ~isempty(parms.roitypes)
    keep_idx = get_roitype_idx(roinames, parms.roitypes);
    roinames = roinames(keep_idx);
    nroi = length(roinames);
    datamat_tsdata = datamat_tsdata(:,:,:,keep_idx);
  end

  % save as mat file
  fprintf('%s: saving results to %s...\n',mfilename,parms.fname_out);
  save(parms.fname_out,'datamat_tsdata','datamat_motion','datamat_onset','datamat_duration',...
                 'datamat_design','censvec','censvec_deconv','roinames','stim_types','stim_labels',...
                 'ntpoints','nroi','nruns','ndirs','nstim','nreg','stimregs',...
                 'nbase','nmotion','nframes','nframes_run','TR',...
                 '-v7');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function keep_idx = get_roitype_idx(roinames, roitypes)
  keep_flag = false(1, length(roinames));
  for i = 1:length(roinames)
    rn = roinames{i};
    is_gp  = ~isempty(regexp(rn, '^ctx-[lr]h-gp\d+$', 'once'));
    is_dsk = ~isempty(regexp(rn, '^ctx-[lr]h-', 'once')) && ~is_gp;
    is_dst = ~isempty(regexp(rn, '^ctx_[lr]h_', 'once'));
    is_aseg = ~is_gp && ~is_dsk && ~is_dst;
    if (any(strcmp(roitypes,'gp'))   && is_gp)   || ...
       (any(strcmp(roitypes,'dsk'))  && is_dsk)  || ...
       (any(strcmp(roitypes,'dst'))  && is_dst)  || ...
       (any(strcmp(roitypes,'aseg')) && is_aseg)
      keep_flag(i) = true;
    end
  end
  keep_idx = find(keep_flag);
return;
