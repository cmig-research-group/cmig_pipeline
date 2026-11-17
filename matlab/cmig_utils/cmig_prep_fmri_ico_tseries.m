function cmig_prep_fmri_ico_tseries(varargin)
%function cmig_prep_fmri_ico_tseries([options])
%
% Purpose: preprocess icosahedral surface-sampled
%          fMRI vertex time series data
%
% Optional Parameters ('key',value pairs):
%   'indir': input directory
%     {default = [pwd '/input']}
%   'instem': input file stem of concatenated ico time series data
%     {default = 'ico_tseries'}
%   'roi_instem': input file stem of concatenated ROI time series data
%     used only when wm_flag, csf_flag, or brain_flag are set
%     {default = 'tseries'}
%   'outdir': output directory
%     {default = [pwd '/output']}
%   'outstem': output file stem
%     {default = 'ico_tseries_prep'}
%   'task_residual_flag': [0|1|2] whether to calculate task-predicted signal or residual
%     0: use actual signal
%     1: calculate task-predicted signal
%     2: calculate task residual
%     if 1 or 2, requires input data include datamat_design (see cmig_concat_fmri_ico_tseries)
%     {default = 0}
%   'preproc_flag': [0|1] whether to preprocess data (e.g., raw tseries data as input)
%     {default = 1}
%   'prefilt_flag': [0|1] whether to apply bandpass filter
%     to time series and motion data before regression
%     unused if preproc_flag = 0
%     {default = 0}
%   'censor_flag': [0|1] whether to censor using censor_thresh or use existing censvec
%     {default = 0}
%   'motion_radius': approximate radius of typical head (mm)
%     for calculating distance from an angle of rotation
%     {default = 50}
%   'motion_absflag': calculate motion as sum of absolute values of each component
%     otherwise, calculate square root of sum of squares
%     {default = 1}
%   'motion_nodyflag': calculate motion ignoring dy (translation in y)
%     {default = 0}
%   'motion_gradflag': calculate framewise displacements using gradient instead of diff
%     {default = 0}
%   'mc_resp_flag': [0|1] notch filter motion estimates for aliased respiration signals
%     unused if preproc_flag=0 and censor_flag=0
%     {default = 0}
%   'resid_censor_thresh': framewise displacement censoring threshold (mm)
%     for pre-residualization
%     {default = 0.3}
%   'censor_thresh': framewise displacement censoring threshold (mm)
%     unused if preproc_flag=0 and censor_flag=0
%     {default = 0.2}
%   'contig_thresh': mininum number of contiguous frames for censoring
%     {default = 5}
%   'outlier_flag': [0|1] exclude frames that are outliers for spatial variation
%     criterion is standard deviation across ROIs
%     threshold is 3 x MAD (median absolute deviation) from median
%     {default = 1}
%   'TR': repetition time (s)
%     used only if preproc_flag = 1
%     {default = 0.8}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Created:  06/22/2026 by Don Hagler
% Prev Mod: 07/02/2026 by Don Hagler
% Prev Mod: 08/10/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 11/02/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on ABCD_concat_roi_rsfMRI_tseries_prep and ABCD_concat_ico_rsfMRI_corr
%   Created:  01/27/2026 by Don Hagler
%   Last Mod: 04/28/2026 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',[pwd '/input'],[],...
  'instem','ico_tseries',[],...
  'roi_instem','tseries',[],...
  'outdir',[pwd '/output'],[],...
  'outstem','ico_tseries_prep',[],...
  'task_residual_flag',0,[],... % [0:2],...
  'preproc_flag',1,[],...
  'prefilt_flag',0,[],...
  'censor_flag',0,[],...
  'motion_radius',50,[],...
  'motion_absflag',true,[false true],...
  'motion_nodyflag',false,[false true],...
  'motion_gradflag',false,[],... % [false true],...
  'mc_resp_flag',false,[],... % [false true],...
  'resid_censor_thresh',0.3,[],...
  'censor_thresh',0.2,[],... % [0.001,Inf],...
  'contig_thresh',5,[],... % [0,20],...
  'outlier_flag',true,[],... % [false true],...
  'TR',0.8,[],...
  'forceflag',false,[],... % [false true],...
... % residualization parameters
  'detrend',2,[],... % [0:2],...
  'mc_flag',true,[],... % [false true],...
  'wm_flag',true,[],... % [false true],...
  'csf_flag',true,[],... % [false true],...
  'brain_flag',true,[],... % [false true],...
  'deriv_flag',true,[],... % [false true],...
  'mc_sq_flag',true,[],... % [false true],...
... % additional parameters
  'interp_flag',true,[false true],...
  'bandpass_flag',true,[false true],...
  'bandpass_low',0.009,[],...
  'bandpass_high',0.08,[],...
  'bandpass_tb',0.001,[],...
  'bandpass_zpad_flag',true,[false true],...
  'demean_flag',true,[false true],...
});

% parse string inputs
%   as number, string, cell array, or numeric array
if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

% check input file
if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end
fname_in = sprintf('%s/%s.mat',parms.indir,parms.instem);
if ~exist(fname_in,'file')
  error('input file %s not found',fname_in);
end

% whether to load ROI data for nuisance regression
roi_flag = parms.wm_flag || parms.csf_flag || parms.brain_flag;

% check roi input file
if roi_flag
  fname_roi = sprintf('%s/%s.mat',parms.indir,parms.roi_instem);
  if ~exist(fname_roi,'file')
    error('roi input file %s not found',fname_roi);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = sprintf('%s/%s.mat',parms.outdir,parms.outstem);
if ~exist(fname_out,'file') || parms.forceflag

  % load concatenated ico fMRI time series data
  fprintf('%s: loading %s...\n',mfilename,fname_in);
  tic; data = load(fname_in); toc;
  nverts = data.nverts;
  nverts_lh = data.nverts_lh;
  nverts_rh = data.nverts_rh;

  % build vertex names
  vtxnames_lh = arrayfun(@(v) sprintf('lhv%d',v), data.verts_lh(:)', 'UniformOutput', false);
  vtxnames_rh = arrayfun(@(v) sprintf('rhv%d',v), data.verts_rh(:)', 'UniformOutput', false);
  vtxnames = [vtxnames_lh, vtxnames_rh];

  if isfield(data,'ndirs')
    ndirs = data.ndirs;
  else
    ndirs = size(data.datamat_tsdata,2);
  end

  % load ROI data for nuisance regression if needed
  if roi_flag
    fprintf('%s: loading %s...\n',mfilename,fname_roi);
    tic; roi_data = load(fname_roi); toc;
    roinames = roi_data.roinames;
  end

  % which vertices have no NaN timepoints
  defmat = isfinite(sum(data.datamat_tsdata,3));
  % NOTE: datamat_tsdata has size = [nruns,ndirs,ntpoints,nverts]

  % loop over visits
  for d=1:ndirs

    % periodic status update
    if d==1 || mod(d,100)==0
      fprintf('%s: diri = %d of %d (now=%s)\n',...
        mfilename,d,ndirs,datestr(now,'HH:MM:SS'));
    end

    % skip subject if no valid data points (all NaNs)
    if ~nnz(mmil_rowvec(defmat(:,d,:,:)))
      fprintf('%s: WARNING: %d/%d : skipping because no valid data...\n',...
        mfilename,d,ndirs);
      continue;
    end

    % compile data into struct for this visit
    tdata = [];
    % loop over runs
    r = 1;
    for s=1:data.nruns
      % skip run if no valid data (all NaNs)
      if ~nnz(squeeze(defmat(s,d,:,:)))
        fprintf('%s: WARNING: %d/%d : skipping run %d because no valid data...\n',...
          mfilename,d,ndirs,s);
        continue;
      end
      % extract ico vertex data for this run
      ico_data = squeeze(data.datamat_tsdata(s,d,:,:));
      % concatenate roi data for nuisance regression
      if roi_flag
        run_roi = squeeze(roi_data.datamat_tsdata(s,d,:,:));
        tdata(r).roi_data = cat(2, ico_data, run_roi);
        tdata(r).roinames = [vtxnames, roinames];
      else
        tdata(r).roi_data = ico_data;
        tdata(r).roinames = vtxnames;
      end
      if parms.task_residual_flag
        if ~isfield(data,'datamat_design')
          error('task_residual_flag = %d, but datamat_design is missing',...
            parms.task_residual_flag);
        end
        tdata(r).X = squeeze(data.datamat_design(s,d,:,:));
        tdata(r).ind_stim = data.stimregs;
      end
      tdata(r).motion_data = squeeze(data.datamat_motion(s,d,:,:));
      tdata(r).TR = parms.TR;
      tdata(r).numTRs = data.ntpoints;
      tdata(r).ntpoints = size(ico_data,1);
      tdata(r).time = [0:(tdata(r).ntpoints-1)]*tdata(r).TR;
      if isempty(ico_data)
        tdata(r).ind_valid = [];
        tdata(r).roi_data_censored = [];
      else
        % pre-computed censoring
        tdata(r).ind_valid = find(~squeeze(data.censvec(s,d,:)));
        % the refit's frames: the deconvolution censors looser than the
        %   correlation and only the refit should use them. ind_valid stays the
        %   correlation's, four other consumers reading it
        tdata(r).ind_valid_deconv = find(~squeeze(data.censvec_deconv(s,d,:)));
        tdata(r).roi_data_censored = tdata(r).roi_data(tdata(r).ind_valid,:);
        [tdata(r).censor_tseries,tdata(r).motion_fd,tdata(r).motion_fdmat,tdata(r).motion_stats] = ...
          mmil_censor_motion(tdata(r).motion_data,...
            parms.censor_thresh,parms.contig_thresh,...
            parms.motion_radius,parms.motion_absflag,...
            parms.motion_nodyflag,parms.motion_gradflag);
      end
      r = r + 1;
    end

    % option to calculate predicted signal or residual (task fMRI only)
    if parms.task_residual_flag
      for r=1:length(tdata)
        % the whole array, nuisance columns included, as the ico original
        %   ABCD_concat_ico_rsfMRI_corr does and cmig_calc_fmri_ico_corr
        %   already did. the ROIs are dropped further down, when the data goes
        %   back into datamat_tsdata, rather than held out of the
        %   substitution: otherwise the white matter, CSF and global signals
        %   abcd_prep_rsfMRI_roi_data regresses out would still carry the task
        % predicted time series from the stimulus regressors alone
        Ys = cmig_task_predicted(tdata(r).roi_data,tdata(r).X,...
          tdata(r).ind_valid_deconv,tdata(r).ind_stim);
        % replace data with predicted data or residual
        if parms.task_residual_flag==1
          tdata(r).roi_data = 100 + Ys;
        elseif parms.task_residual_flag==2
          tdata(r).roi_data = tdata(r).roi_data - Ys;
        end
      end
    end

    % option to preprocess and/or censor data
    if parms.preproc_flag || parms.censor_flag
      tp = [];
      tp.mc_resp_flag = parms.mc_resp_flag;
      tp.detrend = parms.detrend;
      tp.mc_flag = parms.mc_flag;
      tp.wm_flag = parms.wm_flag;
      tp.csf_flag = parms.csf_flag;
      tp.brain_flag = parms.brain_flag;
      tp.deriv_flag = parms.deriv_flag;
      tp.mc_sq_flag = parms.mc_sq_flag;
      tp.interp_flag = parms.interp_flag;
      tp.prefilt_flag = parms.prefilt_flag;
      tp.bandpass_flag = parms.bandpass_flag;
      tp.bandpass_low = parms.bandpass_low;
      tp.bandpass_high = parms.bandpass_high;
      tp.bandpass_tb = parms.bandpass_tb;
      tp.bandpass_zpad_flag = parms.bandpass_zpad_flag;
      tp.resid_censor_thresh = max(parms.resid_censor_thresh,parms.censor_thresh);
      tp.censor_thresh = parms.censor_thresh;
      tp.contig_thresh = parms.contig_thresh;
      tp.outlier_flag = parms.outlier_flag;
      tp.motion_radius = parms.motion_radius;
      tp.motion_absflag = parms.motion_absflag;
      tp.motion_nodyflag = parms.motion_nodyflag;
      tp.motion_gradflag = parms.motion_gradflag;
      if ~parms.preproc_flag
        tp.mc_resp_flag = 0; tp.norm_flag = 0; tp.detrend = 0;
        tp.mc_flag = 0; tp.wm_flag = 0; tp.brain_flag = 0; tp.csf_flag = 0;
        tp.interp_flag = 0; tp.prefilt_flag = 0;
        tp.bandpass_flag = 0; tp.demean_flag = 0;
      end
      ar = mmil_parms2args(tp);
      tdata = abcd_prep_rsfMRI_roi_data(tdata,ar{:});
    end

    % put pre-processed time series data back into data.datamat_tsdata
    % strip roi columns, keep only ico vertex columns (first nverts)
    r = 1;
    for s=1:data.nruns
      % skip run if no valid data
      if ~nnz(mmil_rowvec(defmat(s,d,:,:))), continue; end
      % restore NaN for vertices that were invalid for this run
      % (abcd_prep_rsfMRI_roi_data zeroes NaN before regressing and does not restore it)
      bad_idx = find(~squeeze(defmat(s,d,:,:)));
      tdata(r).roi_data(:,bad_idx) = NaN;
      data.datamat_tsdata(s,d,:,:) = tdata(r).roi_data(:,1:nverts);
      data.datamat_motion(s,d,:,:) = tdata(r).motion_data;
      censor = ones(1,1,data.ntpoints);
      censor(tdata(r).ind_valid) = 0;
      data.censvec(s,d,:) = censor;
      r = r + 1;
    end
  end

  % save as mat file (ico fields only, not roi_data fields from roi tseries)
  fprintf('%s: saving results to %s...\n',mfilename,fname_out);
  save(fname_out,'-struct','data','-v7');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
