function [Ys,D] = cmig_task_predicted(D,X,ind_valid,ind_stim)
%function [Ys,D] = cmig_task_predicted(D,X,ind_valid,ind_stim)
%
% Purpose: fit a deconvolution design on the uncensored frames and return the
%   stimulus-only predicted time series, for building task-predicted or
%   task-residual data
%
% Required Input:
%   D: time series data, [ntpoints,ncols]
%   X: design matrix, [ntpoints,nreg]
%   ind_valid: indices of the frames to fit on
%     the caller chooses the frame set; these are the deconvolution's
%     uncensored frames, not the correlation's
%   ind_stim: indices of the stimulus columns of X
%
% Output:
%   Ys: predicted time series from the stimulus regressors alone,
%     [ntpoints,ncols]
%   D: the input demeaned on ind_valid, as the fit saw it, returned so a
%     caller that wants the residual of the demeaned data need not repeat it
%
% Extracted from four copies in cmig_calc_fmri_{ico,roi}_corr and
%   cmig_prep_fmri_{ico,roi}_tseries, which had drifted apart. The frame set
%   is an argument so that the choice of threshold lives in one place per
%   caller rather than inside this function.
%
% Created:  09/09/2026 by Don Hagler
% Last Mod: 09/09/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % subtract mean of non-censored frames
  D = D - nanmean(D(ind_valid,:),1);
  % censor design matrix
  Xc = X(ind_valid,:);
  % calculate inverse matrix
  Wc = pinv(Xc);
  % get design matrix for stimulus parameters only
  Xs = X(:,ind_stim);
  % calculate beta weights using censored design matrix
  B = Wc*D(ind_valid,:);
  Bs = B(ind_stim,:);
  % predicted data (stimulus only)
  Ys = Xs*Bs;
return;
