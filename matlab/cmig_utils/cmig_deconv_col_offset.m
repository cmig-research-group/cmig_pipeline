function col_offset = cmig_deconv_col_offset(parms,nruns)
%function col_offset = cmig_deconv_col_offset(parms,nruns)
%
% Purpose: count the columns of a 3dDeconvolve bucket that precede the GLT
%   (contrast) block, so callers can index the contrast beta and sem
%
% Required Input:
%   parms: struct with these fields
%     baseline_flag: [0|1] whether the deconvolution kept the raw regressor
%       coefficients, i.e. whether -nocout was dropped
%     stim_fnames: cell array of stim timing file names, one per condition
%     nreg_per_stim: regressors per condition: 1 for SPMG1, GAM and BLOCK,
%       2 for SPMG and SPMG2, stim_times_nbasis for TENT and CSPLIN
%     detrend: polort order passed to 3dDeconvolve
%     fname_motion: motion regressor file name; [] if motion was not modeled
%     motion_labels: cell array of motion regressor labels
%     deriv_flag: [0|1] whether motion derivatives were included
%   nruns: number of runs in the deconvolution
%
% Output:
%   col_offset: number of columns before the GLT block; 0 unless
%     baseline_flag is set
%
% NOTE: with baseline_flag set, 3dDeconvolve drops -nocout and the bucket
%   carries a Coef and a Tstat for every regressor ahead of the contrasts,
%   hence each term below doubled. The bucket order is Full_Fstat, then
%   Run#r Pol#p per run and polort, then nreg_per_stim per condition, then
%   motion and its derivatives, then the GLT block. The sub-brick labels in
%   the stats_info file cannot be matched instead: AFNI truncates them at 31
%   usable characters and the nback contrast names already reach 31.
%
% Created:  09/07/2026 by Don Hagler
% Last Mod: 09/07/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col_offset = 0;
if ~mmil_check_nargs(nargin,2), return; end

if ~parms.baseline_flag, return; end

nconds = length(parms.stim_fnames);
nmotion = 0;
if ~isempty(parms.fname_motion), nmotion = length(parms.motion_labels); end
col_offset = nruns * (parms.detrend+1) * 2 ...
           + nconds * parms.nreg_per_stim * 2 ...
           + nmotion * (1+parms.deriv_flag) * 2;

return;
