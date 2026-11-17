function results = abcd_prep_rsfMRI_roi_data(results,varargin)
%function results = abcd_prep_rsfMRI_roi_data(results,[options])
%
% Required Input:
%   results: struct array with fields
%     roi_data, roinames, and motion_data
%     and number of elements equal to number of scans
%
% Optional Parameters:
%  'skipTRs': number of TRs at beginning of scan to be censored
%    {default = 0}
%  'mc_resp_flag': [0|1] notch filter motion estimates for aliased respiration signals
%    {default = 0}
%  'resp_low': respiration rate low frequency cut-off (respirations/minute)
%    {default = 18.6}
%  'resp_high': respiration rate high frequency cut-off (respirations/minute)
%    {default = 25.7}
%  'resid_censor_thresh': framewise displacement censoring threshold (mm)
%     for pre-residualization
%     {default = 0.3}
%  'resid_contig_thresh': mininum number of contiguous frames required
%     for motion censoring for pre-residualization
%     {default = 0}
%  'motion_radius': approximate radius of typical head (mm)
%    for calculating distance from an angle of rotation
%    {default = 50}
%  'motion_absflag': calculate frame-to-frame motion as sum of abs values
%    otherwise, calculate square root of sum of squares
%    {default = 1}
%  'motion_nodyflag': calculate motion ignoring dy (translation in y)
%    {default = 0}
%  'motion_gradflag': calculate framewise displacements using gradient instead of diff
%    {default = 0}
%  'norm_flag': [0|1] normalize input timeseries by mean for each voxel
%    NOTE: if detrend>0, new mean will be 0, but values will be % change
%    {default = 1}
%  'normval':  value to be assigned to mean
%    {defalut = 100}
%  'thresh': when normalize to mean, set ROI with original mean values
%    less than this to 0
%    {default = 1}
%  'detrend': [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    3: cubic detrend
%    {default = 2}
%  'mc_flag': [0|1] use motion time series as regressors
%    {default = 1}
%  'wm_flag': [0|1] use white matter time series as regressor
%    {default = 1}
%  'csf_flag': [0|1] use CSF time series as regressor
%    {default = 1}
%  'brain_flag': [0|1] use brain time series as regressor
%    {default = 1}
%  'deriv_flag': [0|1] include derivatives of motion, white matter, CSF,
%    and brain as regressors (if mc_flag, wm_flag, csf_flag, or brain_flag)
%    {default = 1}
%  'mc_sq_flag':  [0|1] include squares of motion
%    (if mc_flag, derivatives too if deriv_flag)
%    {default = 1}
%  'interp_flag': [0|1] interpolate censored frames
%    before applying bandpass filter
%    {default = 1}
%  'prefilt_flag': [0|1] whether to apply bandpass filter
%     to time series and motion data before regression
%    {default = 0}
%  'bandpass_flag': [0|1] bandpass filter time series
%    {default = 1}
%  'bandpass_low': low frequency cut-off (Hz)
%    {default = 0.009}
%  'bandpass_high': high frequency cut-off (Hz)
%    {default = 0.08}
%  'bandpass_tb': bandpass filter transition band (Hz)
%    {default = 0.001}
%  'bandpass_zpad_flag': [0|1] pad time series with zeros on each end to
%     supress edge artifacts
%     {default = 1}
%  'censor_thresh': framewise displacement censoring threshold (mm)
%     {default = 0.2}
%  'contig_thresh': mininum number of contiguous frames required
%     for final motion censoring
%     {default = 5}
%  'outlier_flag': [0|1] exclude frames that are outliers for spatial variation
%     criterion is standard deviation across ROIs
%     threshold is 3 x MAD (median absolute deviation) from median
%     {default = 1}
%  'demean_flag': [0|1] subtract mean for each ROI calculated without censored frames
%     {default = 1}
%
% Created:  01/07/21 by Don Hagler
% Prev Mod: 09/15/24 by Don Hagler
% Prev Mod: 05/05/25 by Don Hagler
% Prev Mod: 01/27/26 by Don Hagler
% Prev Mod: 01/30/26 by Don Hagler
% Last Mod: 04/28/26 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'skipTRs',0,[0 Inf],...
... % motion filtering
  'mc_resp_flag',false,[false true],...
  'resp_low',18.6,[],...
  'resp_high',25.7,[],...
... % pre-residualization censoring parameters
  'resid_censor_thresh',0.3,[],...
  'resid_contig_thresh',0,[],...
  'motion_radius',50,[],...
  'motion_absflag',true,[false true],...
  'motion_nodyflag',false,[false true],...
  'motion_gradflag',false,[false true],...
... % normalization parameters
  'norm_flag',true,[false true],...
  'normval',100,[1 1e10],...
  'thresh',1,[0 Inf],...
... % residualization parameters
  'detrend',2,[0:2],...
  'mc_flag',true,[false true],...
  'wm_flag',true,[false true],...
  'csf_flag',true,[false true],...
  'brain_flag',true,[false true],...
  'deriv_flag',true,[false true],...
  'mc_sq_flag',true,[false true],...
  'wm_names',{'Left-Cerebral-White-Matter',...
              'Right-Cerebral-White-Matter'},[],...
  'csf_names',{'Left-Lateral-Ventricle',...
               'Right-Lateral-Ventricle'},[],...
... % filtering parameters
  'interp_flag',true,[false true],...
  'prefilt_flag',false,[false true],...
  'bandpass_flag',true,[false true],...
  'bandpass_low',0.009,[],...
  'bandpass_high',0.08,[],...
  'bandpass_tb',0.001,[],...
  'bandpass_zpad_flag',true,[false true],...
  'TR',0.8,[],...
... % final censoring parameters
  'censor_thresh',0.2,[0.001,Inf],...
  'contig_thresh',5,[0,20],...
  'outlier_flag',true,[false true],...
... % demean option
  'demean_flag',true,[false true],...
... % hidden options
  'verbose',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over scans
for i=1:length(results)
  TR = mmil_getfield(results(i),'TR',parms.TR);
  % force motion_data to double (avoid problem with pinv)
  results(i).motion_data = double(results(i).motion_data);
  % filter respiratory signals from motion
  if parms.mc_resp_flag
    results(i).motion_data = filter_motion(results(i).motion_data,TR,parms);
  end
  % set valid frames for pre-residualization
  results(i).ind_valid_resid = ...
    censor_resid(results(i).motion_data,parms);
  % normalize data
  if parms.norm_flag
    results(i).roi_data = ...
      norm_data(results(i).roi_data,...
                results(i).ind_valid_resid,parms);
  end
  
  if parms.prefilt_flag
    results(i).motion_data_filt = results(i).motion_data;
    % interpolate censored frames before filtering
    if parms.interp_flag
      results(i).roi_data = interp_censored(results(i).roi_data,...
                                 results(i).ind_valid_resid);
      results(i).motion_data_filt = interp_censored(results(i).motion_data_filt,...
                                 results(i).ind_valid_resid);
    end
    % apply bandpass filter to roi_data and motion_data
    if parms.bandpass_flag
      results(i).roi_data = filter_data(results(i).roi_data,TR,parms);
      results(i).motion_data_filt = filter_data(results(i).motion_data_filt,TR,parms);
    end
    % residualize filtered data
    if parms.detrend || parms.mc_flag ||...
       parms.wm_flag || parms.brain_flag || parms.csf_flag
      results(i).roi_data = ...
        resid_data(results(i).roi_data,...
                  results(i).roinames,...
                  results(i).motion_data_filt,...
                  results(i).ind_valid_resid,parms);
    end
  else
    % residualize data
    if parms.detrend || parms.mc_flag ||...
       parms.wm_flag || parms.brain_flag || parms.csf_flag
      results(i).roi_data = ...
        resid_data(results(i).roi_data,...
                  results(i).roinames,...
                  results(i).motion_data,...
                  results(i).ind_valid_resid,parms);
    end
    % interpolate censored frames
    if parms.interp_flag
      results(i).roi_data = interp_censored(results(i).roi_data,...
                                 results(i).ind_valid_resid);
    end
    % apply bandpass filter to roi_data
    if parms.bandpass_flag
      results(i).roi_data = filter_data(results(i).roi_data,TR,parms);
    end
  end
  
  % censor data
  [results(i).roi_data_censored,results(i).ind_valid,...
   results(i).motion_fd,results(i).motion_fdmat,...
   results(i).motion_stats,results(i).censor_tseries] = ...
    censor_data(results(i).roi_data,...
                results(i).motion_data,parms);

  % subtract mean from timeseries for each ROI
  if parms.demean_flag
    results(i).roi_data = ...
      demean_data(results(i).roi_data,results(i).ind_valid);
    results(i).roi_data_censored = results(i).roi_data(results(i).ind_valid,:);
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function motion_data = filter_motion(motion_data,TR,parms)
  motion_data = abcd_resp_filter_motion(motion_data,TR,...
    parms.resp_low,parms.resp_high);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind_valid = censor_resid(motion_data,parms)
  censor_tseries_resid = mmil_censor_motion(motion_data,...
    parms.resid_censor_thresh,...
    parms.resid_contig_thresh,...
    parms.motion_radius,...
    parms.motion_absflag,...
    parms.motion_nodyflag,...
    parms.motion_gradflag);
  if parms.skipTRs, censor_tseries_resid(1:parms.skipTRs) = 1; end;
  ind_valid = find(~censor_tseries_resid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = norm_data(roi_data,ind_valid,parms)
  nframes = size(roi_data,1);
  roi_mean = repmat(mean(roi_data(ind_valid,:),1),[nframes,1]);
  roi_data = parms.normval*roi_data./(roi_mean+eps);
  roi_data(roi_mean<parms.thresh) = 0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = resid_data(roi_data,roinames,motion_data,ind_valid,parms)
  % force data type to double
  roi_data = double(roi_data);
  % set NaNs to 0
  roi_data(~isfinite(roi_data)) = 0;  
  % construct design matrix
  X = design_matrix(roi_data,roinames,motion_data,parms);
  % calculate beta coefficients
  V = roi_data';
  if isempty(ind_valid)
    beta = (pinv(X)*V')';
  else
    beta = (pinv(X(ind_valid,:))*V(:,ind_valid)')';
  end;
  % calculate model fit
  Y = beta*X';
  % subtract fit from data
  V = V - Y;
  % transpose to original shape
  roi_data = V';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = design_matrix(roi_data,roinames,motion_data,parms)
  X = [];

  % quadratic detrending
  nframes = size(roi_data,1);
  t = [0:nframes-1]';
  X = ones(size(t));
  if parms.detrend
    for d=1:parms.detrend
      X = [X t.^d];
    end;
  end;

  % motion
  if parms.mc_flag
    m = motion_data;
    X = cat(2,X,m);
    if parms.deriv_flag
      % motion derivatives
      if parms.motion_gradflag
        dm = gradient(m);
      else
        dm = cat(1,zeros(1,6),diff(m));
      end
      X = cat(2,X,dm);
      % squard motion and squared motion derivatives
      if parms.mc_sq_flag
        X = cat(2,X,m.^2,dm.^2);
      end
    elseif parms.mc_sq_flag
      % squared motion
      X = cat(2,X,m.^2);
    end
  end

  % white matter
  if parms.wm_flag
    % get wm roi indices by name
    [~,i_wm] = intersect(roinames,parms.wm_names);
    % add wm time series and derivative to design matrix
    X = add_roi_to_design(X,i_wm,roi_data,'wm',parms);
  end
  
  % brain
  if parms.brain_flag
    % use all rois for global "brain" regressor
    i_bm = [1:length(roinames)];
    % add brain time series and derivative to design matrix
    X = add_roi_to_design(X,i_bm,roi_data,'brain',parms);
  end
  
  % CSF
  if parms.csf_flag
    % get csf roi indices by name
    [~,i_csf] = intersect(roinames,parms.csf_names);
    % add brain time series and derivative to design matrix
    X = add_roi_to_design(X,i_csf,roi_data,'csf',parms);
  end  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = interp_censored(roi_data,ind_valid)
  int_data = roi_data;
  nframes = size(roi_data,1);
  censor_tseries = ones(nframes,1);
  censor_tseries(ind_valid) = 0;
  % identify censored and uncensored frames
  ind_censor = find(censor_tseries);
  if isempty(ind_censor), return; end;
  if length(ind_valid)>1
    % interpolate censored frames
    tmp_data = interp1(ind_valid,roi_data(ind_valid,:),ind_censor);
    int_data(ind_censor,:) = tmp_data;
    % if last frame is censored, use mean across frames
    roi_mean = mean(roi_data(ind_valid,:),1);
    if censor_tseries(end), int_data(end,:) = roi_mean; end;
    % relace censored frames with NaNs
    %   with mean across uncensored frames
    ind_valid_rois = find(~isnan(roi_mean));
    ind_nans = find(isnan(mean(tmp_data(:,ind_valid_rois),2)));
    if ~isempty(ind_nans)
      int_data(ind_censor(ind_nans),ind_valid_rois) = repmat(roi_mean(ind_valid_rois),numel(ind_nans),1);
    end;
  elseif length(ind_valid)==1
    int_data(ind_censor,:) = repmat(roi_data(ind_valid,:),length(ind_censor),1);
  else
    int_data(ind_censor,:) = 0;
  end;
  roi_data = int_data;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = filter_data(roi_data,TR,parms)
  roi_data = ts_freq_filt(roi_data,1/TR,...
    [parms.bandpass_low,parms.bandpass_high],...
    [parms.bandpass_tb,parms.bandpass_tb],'bandpass',...
    parms.bandpass_zpad_flag);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [roi_data,ind_valid,motion_fd,motion_fdmat,motion_stats,censor_tseries] = censor_data(roi_data,motion_data,parms)
  ind_valid = [];
  % censor frames based on FD
  [censor_tseries,motion_fd,motion_fdmat,motion_stats] = mmil_censor_motion(motion_data,...
    parms.censor_thresh,parms.contig_thresh,...
    parms.motion_radius,parms.motion_absflag,...
    parms.motion_nodyflag,parms.motion_gradflag);
  % censor frames based on outliers
  ind_valid = find(~censor_tseries);
  if parms.outlier_flag && ~isempty(ind_valid)
    % exclude outliers based on standard deviation across ROIs  
    tmp_data = roi_data(ind_valid,:);
    ind_valid_roi = find(~isnan(sum(tmp_data,1)));
    std_vals = std(tmp_data(:,ind_valid_roi),[],2);
    [mad_val,med_val] = mmil_wtd_mad(std_vals);
    % threshold is 3 x MAD (median absolute deviation) from median
    ind_keep = find(std_vals < med_val + 3*mad_val &...
                    std_vals > med_val - 3*mad_val);
    ind_valid = ind_valid(ind_keep);
  end;
  roi_data = roi_data(ind_valid,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = demean_data(roi_data,ind_valid)
  roi_data = roi_data - nanmean(roi_data(ind_valid,:),1);
  roi_data(~isfinite(roi_data)) = 0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = add_roi_to_design(X,i_roi,roi_data,roi_name,parms)
  % exclude ROIs with NaNs
  i_roi = setdiff(i_roi,find(isnan(mean(roi_data,1))));
  if ~isempty(i_roi)
    % average time series
    x_roi = mean(roi_data(:,i_roi),2);
    X = cat(2,X,x_roi);
    % calculate derivative
    if parms.deriv_flag
      if parms.motion_gradflag
        d_roi = gradient(x_roi);
      else
        d_roi = [0;diff(x_roi)];
      end
      X = cat(2,X,d_roi);
    end
  elseif parms.verbose
    fprintf('%s: NOTE: excluded %s ROI time series containing NaNs from design matrix\n',...
      mfilename,roi_name);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


