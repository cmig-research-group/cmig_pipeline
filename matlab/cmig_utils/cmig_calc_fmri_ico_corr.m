function cmig_calc_fmri_ico_corr(varargin)
%function cmig_calc_fmri_ico_corr([options])
%
% Purpose: calculate icosahedral surface vertex-to-vertex resting-state fMRI
%          correlations from concatenated time series data
%
% Optional Parameters ('key',value pairs):
%   'indir': input directory
%     {default = [pwd '/input']}
%   'instem': input file fstem of concatenated time series data
%     {default = 'ico_tseries'}
%   'infostem': input file stem of info file
%     {default = 'vol_info'}
%   'roi_instem': input file fstem of concatenated ROI time series data
%     {default = 'roi_tseries'}
%   'outdir': output directory
%     {default = [pwd '/output']}
%   'outstem': output file stem of concatenated correlation matrices
%     {default = 'ico_corrmat'}
%   'task_residual_flag': [0|1|2] whether to calculate task-predicted signal or residual
%     0: use actual signal
%     1: calculate task-predicted signal
%     2: calculate task residual
%     if 1 or 2, requires input data include datamat_design (see ABCD_concat_roi_tfMRI_tseries)
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
%   'min_numTRs': minimum number of TRs (before censoring) required to analyze scan
%     {default = 100}
%   'min_ntpoints': minimum number of TRs after censoring required to include scan
%     {default = 64}
%   'trim_ntpoints': maximum number of TRs after censoring allowed (total across scans)
%     time series with more time points will be trimmed before calculating correlations
%     results for visits with fewer time points will be set to NaN
%     {default = Inf}
%   'trim_rand_flag': [0|1] randomly select trim_ntpoints frames; if 0, use first N
%     {default = 0}
%   'trim_seed': seed of the random number generator for trim_rand_flag=1
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Created:  06/22/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 08/10/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 10/13/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on ABCD_concat_ico_rsfMRI_corr
%   Created:  10/22/2025 by Don Hagler
%   Last Mod: 04/24/2026 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',[pwd '/input'],[],...
  'instem','ico_tseries',[],...
  'infostem','vol_info',[],...
  'roi_instem','roi_tseries',[],...
  'outdir',[pwd '/output'],[],...
  'outstem','ico_corrmat',[],...
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
  'min_numTRs',100,[],... % [1,1000],...
  'min_ntpoints',64,[],... % [1,1000],...
  'trim_ntpoints',Inf,[],... % [1,Inf],...
  'trim_rand_flag',false,[false true],...
  'trim_seed',0,[0,2^32-1],...
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

% check roi data
if parms.preproc_flag && (parms.wm_flag || parms.csf_flag || parms.brain_flag)
  parms.roi_flag = 1;
  fname_roi = sprintf('%s/%s.mat',parms.indir,parms.roi_instem);
  if ~exist(fname_roi,'file')
    error('ROI file %s not found',fname_roi);
  end
else
  parms.roi_flag = 0;
end

% check info file
fname_info = sprintf('%s/%s.mat',parms.indir,parms.infostem);
%if ~exist(fname_info,'file')
%  error('info file %s not found',fname_info);
%end

% seed the random number generator for reproducible random frame selection
if parms.trim_ntpoints<Inf && parms.trim_rand_flag, rng(parms.trim_seed,'twister'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = sprintf('%s/%s.mat',parms.outdir,parms.outstem);
if ~exist(fname_out,'file') || parms.forceflag

  % load concatenated rsfMRI info file
  if exist(fname_info,'file')
    fprintf('%s: loading %s...\n',mfilename,fname_info);
    tic; vol_info = load(fname_info); toc;
  else
    vol_info = [];
  end

  % load concatenated rsfMRI ico time series data
  fprintf('%s: loading %s...\n',mfilename,fname_in);
  tic; data = load(fname_in); toc;
  % create array of vertex names
  vtxnames = cat(2,...
    cellfun(@(x) sprintf('vtx-lh-%d',x),num2cell(data.verts_lh),'UniformOutput',false)',...
    cellfun(@(x) sprintf('vtx-rh-%d',x),num2cell(data.verts_rh),'UniformOutput',false)');
  nverts = data.nverts;
  ndirs = data.ndirs;

  % load concatenated rsfMRI roi time series data
  if parms.roi_flag
    fprintf('%s: loading %s...\n',mfilename,fname_roi);
    tic; roi_data = load(fname_roi); toc;
    roinames = roi_data.roinames;
    nroi = length(roi_data.roinames);
    roi_ndirs = roi_data.ndirs;
    if roi_ndirs ~= ndirs
      error('ROI data has %d rows but ico data has %d rows',roi_ndirs,ndirs);
    end
  end

  % upper triangle of nverts x nverts correlation matrix
  idx_triu = find(triu(ones(nverts,nverts)));
  nnodes = length(idx_triu);

  % vertex indices for each node
  vert1mat = repmat([1:nverts]',[1,nverts]);
  vert2mat = repmat([1:nverts],[nverts,1]);
  vert1vec = vert1mat(idx_triu);
  vert2vec = vert2mat(idx_triu);

  % initialize outputs
  corrmat = nan(ndirs,nnodes,'single');
  varmat = nan(ndirs,nverts,'single');
  nvolvec = nan(ndirs,1,'single');
  nvolfdsubvec = nan(ndirs,1,'single');
  nvolfdsubcontvec = nan(ndirs,1,'single');
  ntpointvec = nan(ndirs,1,'single');
  meanfdvec = nan(ndirs,1,'single');
  meanfdmat = nan(ndirs,6,'single');
  meanfdvec_uncensored = nan(ndirs,1,'single');
  meanfdmat_uncensored = nan(ndirs,6,'single');
  meanrotvec = nan(ndirs,1,'single');
  meantransvec = nan(ndirs,1,'single');
  maxrotvec = nan(ndirs,1,'single');
  maxtransvec = nan(ndirs,1,'single');
  maxfdvec = nan(ndirs,1,'single');

  % which vertices have no NaN timepoints
  defmat = isfinite(sum(data.datamat_tsdata,3));

  % loop over visits
  for d=1:ndirs

    % periodic status update
    if d==1 || mod(d,100)==0
      fprintf('%s: diri = %d of %d (now=%s)\n',...
        mfilename,d,ndirs,datestr(now,'HH:MM:SS'));
    end

    if isempty(vol_info)
      VisitID = [];
    else
      VisitID = vol_info.visitidvec{d};
    end

    % skip subject if no valid data points
    if ~nnz(mmil_rowvec(defmat(:,d,:,:)))
      fprintf('%s: WARNING: %d/%d : skipping %s because no valid data...\n',...
        mfilename,d,ndirs,VisitID);
      continue;
    end

    % compile data into struct for this visit
    tdata = [];
    % loop over runs
    r = 1;
    for s=1:data.nruns
      % skip run if no valid data
      if ~nnz(squeeze(defmat(s,d,:,:)))
        fprintf('%s: WARNING: %d/%d : skipping run %d for %s because no valid data...\n',...
          mfilename,d,ndirs,s,VisitID);
        continue;
      end
      % extract data for this run
      tdata(r).roi_data = squeeze(data.datamat_tsdata(s,d,:,:));
      % concatenate ROI data with ico data for GSR
      if parms.roi_flag
        tdata(r).roi_data = cat(2,tdata(r).roi_data,squeeze(roi_data.datamat_tsdata(s,d,:,:)));
        % combine vtxnames and roinames
        tdata(r).roinames = cat(2,vtxnames,roinames);
      else
        tdata(r).roinames = vtxnames;
      end
      % residualize task events
      if parms.task_residual_flag
        if ~isfield(data,'datamat_design')
          error('task_residual_flag = %d, but datamat_design is missing',...
            parms.task_residual_flag);
        end
        tdata(r).X = squeeze(data.datamat_design(s,d,:,:));
        tdata(r).ind_stim = data.stimregs;
      end
      tdata(r).motion_data = squeeze(data.datamat_motion(s,d,:,:));
      tdata(r).TR = data.TR;
      tdata(r).numTRs = data.ntpoints;
      tdata(r).ntpoints = size(tdata(r).roi_data,1);
      tdata(r).time = [0:(tdata(r).ntpoints-1)]*tdata(r).TR;
      if isempty(tdata(r).roi_data)
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
        [tdata(r).motion_stats,tdata(r).motion_fd,~,tdata(r).motion_fdmat] = mmil_motion_stats(tdata(r).motion_data,...
          parms.motion_radius,parms.motion_absflag,parms.censor_thresh,parms.motion_nodyflag,[],parms.motion_gradflag);
      end
      r = r + 1;
    end

    % option to calculate predicted signal or residual (task fMRI only)
    if parms.task_residual_flag
      for r=1:length(tdata)
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

    % concatenate time series
    D = [];
    fd_uncensored = []; fd_censored = [];
    fdmat_uncensored = []; fdmat_censored = [];
    mask_uncensored = []; mask_censored = [];
    nvols = 0; nvolfdsub = 0; nvolfdsubcont = 0;
    for r=1:length(tdata)
      tD = tdata(r).roi_data_censored;
      if isempty(tD), continue; end;
      % remove ROIs from processed time series data
      if parms.roi_flag, tD = tD(:,1:nverts); end
      %% NOTE: first frame FD always 0 if gradflag = 0
      %%       exclude from calculation of mean FD for compatibility with mmil_motion_stats
      % create masking variable to exclude first frame
      tm_uncensored = ~ismember([1:size(tdata(r).roi_data,1)]',1);
      tm_censored = ~ismember(tdata(r).ind_valid,1);
      % get uncensored and censored fd
      tfd_uncensored = tdata(r).motion_fd;
      tfd_censored = tdata(r).motion_fd(tdata(r).ind_valid);
      % get uncensored and censored motion_fdmat
      tfdmat_uncensored = tdata(r).motion_fdmat;
      tfdmat_censored = tdata(r).motion_fdmat(tdata(r).ind_valid,:);
      if isfield(tdata,'motion_stats')
        % number of frames
        nvols = nvols + tdata(r).motion_stats.nvols;
        % number of frames with FD < threshold
        nvolfdsub = nvolfdsub + tdata(r).motion_stats.subthresh_nvols;
      end
      if isfield(tdata,'censor_tseries')
        % number of frames with fd < censor_thresh and
        %   at least contig_thresh contiguous uncensored frames
        nvolfdsubcont = nvolfdsubcont + length(find(tdata(r).censor_tseries==0));
      end
      % exclude very short scans
      if tdata(r).numTRs<parms.min_numTRs, continue; end
      % exclude scans with very few valid time points
      if length(tdata(r).ind_valid)<parms.min_ntpoints, continue; end
      D = cat(1,D,tD);
      % concatenate uncensored and censored motion_fd
      fd_uncensored = cat(1,fd_uncensored,tfd_uncensored);
      fd_censored = cat(1,fd_censored,tfd_censored);
      % concatenate uncensored and censored motion_fdmat
      fdmat_uncensored = cat(1,fdmat_uncensored,tfdmat_uncensored);
      fdmat_censored = cat(1,fdmat_censored,tfdmat_censored);
      % concatenate uncensored and censored motion masks
      mask_uncensored = cat(1,mask_uncensored,tm_uncensored);
      mask_censored = cat(1,mask_censored,tm_censored);
    end
    ntpoints = size(D,1);

    if ~ntpoints
      fprintf('%s: WARNING: %d/%d : skipping %s because not enough valid data points...\n',...
        mfilename,d,ndirs,VisitID);
      continue;
    end
    if parms.trim_ntpoints<Inf
      if ntpoints < parms.trim_ntpoints
        % exclude if not enough time points
        fprintf('%s: WARNING: %d/%d : skipping %s because fewer than %d valid data points...\n',...
          mfilename,d,ndirs,VisitID,parms.trim_ntpoints);
        continue;
      else
        % select uniform number of time points (random or first N)
        if parms.trim_rand_flag
          idx_trim = randperm(ntpoints);
        else
          idx_trim = 1:ntpoints;
        end
        D = D(idx_trim(1:parms.trim_ntpoints),:);
        fd_censored = fd_censored(idx_trim(1:parms.trim_ntpoints),:);
        fdmat_censored = fdmat_censored(idx_trim(1:parms.trim_ntpoints),:);
        mask_censored = mask_censored(idx_trim(1:parms.trim_ntpoints),:);
        ntpoints = parms.trim_ntpoints;
      end
    end

    % calculate correlation matrix and variance
    R = corr(D);
    V = var(D,1,1);
    % compile output
    rvec = R(idx_triu);
    corrmat(d,:) = reshape(rvec,[1,nnodes]);
    varmat(d,:) = V;
    % number of frames after skipTR
    nvolvec(d) = nvols;
    % number of frames with fd < censor_thresh
    nvolfdsubvec(d) = nvolfdsub;
    % number of frames with fd < censor_thresh and at least 5 contiguous uncensored frames
    nvolfdsubcontvec(d) = nvolfdsubcont;
    % number of frames after censoring (fd, contig, and outlier)
    ntpointvec(d) = ntpoints;
    % calculate mean motion (framewise displacement) of censored frames
    meanfdvec(d) = mean(fd_censored(find(mask_censored)));
    meanfdvec_uncensored(d) = mean(fd_uncensored(find(mask_uncensored)));
    % calculate components of mean motion (dx,dy,dz,rx,ry,rz)
    meanfdmat(d,:) = mean(fdmat_censored(find(mask_censored),:),1);
    meanfdmat_uncensored(d,:) = mean(fdmat_uncensored(find(mask_uncensored),:),1);
    % calculate mean rotation (uncensored)
    meanrotvec(d) = mean(sum(fdmat_uncensored(find(mask_uncensored),4:6),2),1);
    % calculate mean translation
    meantransvec(d) = mean(sum(fdmat_uncensored(find(mask_uncensored),1:3),2),1);
    % calculate max rotation
    maxrotvec(d) = max(sum(fdmat_uncensored(find(mask_uncensored),4:6),2),[],1);
    % calculate max translation
    maxtransvec(d) = max(sum(fdmat_uncensored(find(mask_uncensored),1:3),2),[],1);
    % calculate max motion
    maxfdvec(d) = max(fd_uncensored(find(mask_uncensored)));
  end

  % save output file
  fprintf('%s: saving results to %s...\n',mfilename,fname_out);
  save(fname_out,'corrmat','varmat','ntpointvec',...
                 'nvolvec','nvolfdsubvec','nvolfdsubcontvec',...
                 'meanfdvec','meanfdmat',...
                 'meanfdvec_uncensored','meanfdmat_uncensored',...
                 'meanrotvec','meantransvec','maxrotvec','maxtransvec','maxfdvec',...
                 'nverts','nnodes','vtxnames','vert1vec','vert2vec','ndirs',...
                 '-v7');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

