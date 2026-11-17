function [line_score,line_nslices,line_extent] = autoqc_detect_lines_fft(vol,outdir,outstem,verbose,plotflag,brainmask)
%function [line_score,line_nslices,line_extent] = autoqc_detect_lines_fft(vol,outdir,outstem,verbose,plotflag,brainmask)
%
% Purpose: Determine if scan has line artifacts (horizontal banding)
%
% The detector scans all in-brain slices along dim1.
%
% Required Input:
%   vol: 3D scan volume
%
% Optional Input:
%   outdir: output directory
%     {default = pwd}
%   outstem: output file stem
%     {default = 'vol'}
%   verbose: [0|1] display messages about lines detected
%     {default = 0}
%   plotflag: [0|1] plot image of the most-affected slice if lines detected
%     {default = 1}
%
% Output:
%   line_score: scalar measure of line artifact (90th percentile of the
%     per-slice Nyquist-band score across in-brain slices)
%   line_nslices: number of in-brain slices with detected banding
%   line_extent: column-coherence of the banding on the most-affected slice -
%     fraction of in-brain columns whose Nyquist-band ripple is in-phase
%     (correlated) with the row-mean ripple. A coherent full-width line -> ~1;
%     incoherent speckle/noise -> ~0. Distinguishes real banding from the
%     row-to-row noise that can also raise line_score.
%   line_score/line_nslices/line_extent are NaN when the scan is UNEVALUABLE (no
%     scorable slice, even after the intensity-mask fallback) - distinct from 0,
%     which means "evaluated, no line artifact".
%
% If a supplied brainmask scores no slice (e.g. SynthStrip failed on a severe
% artifact), the detector automatically falls back to its intensity mask.
%
% Method:
%   Line artifacts are regular, evenly-spaced horizontal lines that appear at
%   the Nyquist frequency (period ~2 rows) in the phase-encode/readout
%   reformats. For each dim1 slice, the volume is reduced to a row-mean profile
%   over brain voxels (largest connected component); the smooth anatomical
%   envelope is removed (movmean high-pass); and the RMS amplitude of the
%   ripple in the Nyquist band (period 2-3.5 rows) is measured. Coherent
%   banding survives the row-mean while incoherent anatomy averages out.
%   Non-brain / speckle slices are excluded so out-of-brain noise cannot score.
%
% Created:   08/29/21 by Payam Saisan
% Prev Mod:  11/29/21 by Payam Saisan
% Prev Mod:  06/26/2026 by Don Hagler
% Prev Mod:  07/10/2026 by Don Hagler
% Last Mod:  07/13/2026 by Don Hagler
%

line_score = 0; line_nslices = 0; line_extent = 0;

if ~mmil_check_nargs(nargin,1), return; end
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end
if ~exist('outstem','var') || isempty(outstem), outstem = 'vol'; end
if ~exist('verbose','var') || isempty(verbose), verbose = 0; end
if ~exist('plotflag','var') || isempty(plotflag), plotflag = 1; end
if ~exist('brainmask','var'), brainmask = []; end   % 3D SynthStrip mask; [] -> intensity mask

% detection parameters
SPAN        = 9;          % dome-removal moving-average span (rows)
MIN_AREA    = 0.08;       % min brain fraction of a slice to be scored
PMIN        = 2.0;        % Nyquist band, min period (rows)
PMAX        = 3.5;        % Nyquist band, max period (rows)
SLICE_THRESH = 8.0;       % per-slice score above which a slice counts as banded
CTHRESH     = 0.5;        % min column-vs-row-mean ripple correlation to count as coherent

% allow single 2D slice as input
if length(size(vol))<=2
  tmpvol(1,:,:) = squeeze(vol);
  vol = tmpvol;
end

% if multi-frame, use middle frame (callers normally pass one frame at a time)
if length(size(vol))==4
  vol = vol(:,:,:,round(end/2));
end

% score every dim1 slice; non-brain slices come back NaN and are dropped below
slice_scores = all_slice_scores(vol,brainmask,SPAN,MIN_AREA,PMIN,PMAX);
valid = slice_scores(~isnan(slice_scores));

% AUTO-FALLBACK: a supplied (SynthStrip) mask can fail on severe-artifact scans -
% ghosting confuses SynthStrip into a tiny wrong mask that gates out every slice.
% If nothing scored, retry with the intensity mask so the artifact is not missed.
if isempty(valid) && ~isempty(brainmask)
  if verbose
    fprintf('%s: supplied mask scored no slice; falling back to intensity mask\n',mfilename);
  end
  brainmask = [];
  slice_scores = all_slice_scores(vol,brainmask,SPAN,MIN_AREA,PMIN,PMAX);
  valid = slice_scores(~isnan(slice_scores));
end

% still nothing scorable -> the scan is UNEVALUABLE. Return NaN (not 0) so the
% caller can tell "could not measure" apart from "measured, no line artifact".
if isempty(valid)
  if verbose, fprintf('%s: no scorable slice; returning NaN\n',mfilename); end
  line_score = NaN; line_nslices = NaN; line_extent = NaN;
  return;
end

% reduce across in-brain slices: p90 is robust to a single outlier slice;
% nslices = how many slices show banding above SLICE_THRESH
line_score = prctile(valid,90);
line_nslices = sum(valid > SLICE_THRESH);

% most-affected slice; measure the column-coherence extent of its banding
[~,smax] = max(slice_scores);
if isempty(brainmask), msx=[]; else msx=squeeze(brainmask(smax,:,:)); end
line_extent = coherence_extent(squeeze(vol(smax,:,:)),msx,SPAN,MIN_AREA,PMIN,PMAX,CTHRESH);

if verbose
  fprintf('%s: line_score=%.3f line_extent=%.2f (n in-brain slices=%d, n banded=%d)\n',...
    mfilename,line_score,line_extent,numel(valid),line_nslices);
end

% plot the most-affected slice, with the detected banded rows marked
if plotflag && line_score > 0
  mmil_mkdir(outdir);
  img = local_normalize(rot90(squeeze(vol(smax,:,:))));
  figure; clf; colormap(bone);
  imagesc(img); axis image off; hold on;

  % overlay the line locations: rows where the Nyquist-band ripple peaks (yellow
  % ticks at the row ends so the brain stays visible)
  rows = banded_rows(squeeze(vol(smax,:,:)),msx,SPAN,MIN_AREA,PMIN,PMAX,10);
  W = size(img,2);
  for rr = rows(:)'
    plot([1:8, W-7:W], rr*ones(1,16),'y.','MarkerSize',6);
  end

  title(sprintf('%s: slice=%d line\\_score=%0.2f nslices=%d extent=%0.2f',...
    regexprep(outstem,'_',' '),smax,line_score,line_nslices,line_extent),'FontSize',8);
  outfile = sprintf('%s/%s_lines_slice%03d.jpg',outdir,outstem,smax);
  mmil_save_fig(gcf,outfile);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function score = score_slice(slice2d,maskslice,span,min_area,pmin,pmax)

  % returns NaN for a non-brain slice (so the caller can drop it), else the
  % Nyquist-band amplitude of the dome-removed row-mean ripple
  score = NaN;

  % orient + normalize the slice to [0,1]
  imgo = local_normalize(rot90(squeeze(slice2d)));
  if max(imgo(:))<=0, return; end

  if isempty(maskslice)
    % fallback (no brain mask supplied): intensity threshold, erode, keep the
    % LARGEST connected component (drops speckle/edge blobs), fill holes
    bw = medfilt2(imgo,[3 3]) > 0.15*max(imgo(:));
    bw = imerode(bw,strel('disk',2));
    if ~any(bw(:)), return; end
    bw = keep_largest_cc(bw);
    bw = imfill(bw,'holes');
  else
    % use the supplied (SynthStrip) brain mask, rot90'd to align with imgo
    bw = rot90(squeeze(maskslice)) > 0.5;
    if ~any(bw(:)), return; end
  end

  % gate: skip slices without enough brain (out-of-brain noise can't score)
  if mean(bw(:)) < min_area, return; end

  % row-mean over brain voxels only (NaN elsewhere so nanmean ignores it)
  m = imgo; m(~bw) = NaN;
  r = nanmean(m,2);

  % keep only interior rows (>25% of the widest brain row) to avoid thin edges
  brainw = sum(bw,2);
  idx = find(brainw > 0.25*max(brainw));
  if numel(idx) < 3*span, return; end
  seg = r(idx);

  % normalize by median brain intensity -> scale-invariant fractional modulation
  med = median(seg,'omitnan');
  if med<=0 || isnan(med), return; end

  % HIGH-PASS: subtract the smooth anatomical dome (movmean) -> banding ripple only
  rd = (seg - movmean(seg,span,'omitnan')) / med;
  rd(isnan(rd)) = 0;

  % power spectrum of the (windowed, demeaned) ripple along the row direction
  L = numel(rd);
  x = (rd - mean(rd)) .* hann(L);   % Hann window reduces spectral leakage
  P = abs(fft(x)).^2;
  half = floor(L/2);
  k = (2:half)';               % one-sided; drop DC and lowest bin (dome residual)
  Pk = P(k+1);
  period = L ./ k;             % rows-per-cycle for each frequency bin

  % SCORE = RMS amplitude in the Nyquist band (period ~2-3.5 rows). Line
  % artifacts are regular every-other-line structure (period ~2); anatomy/noise
  % is low-frequency and falls outside this band.
  innyq = period>=pmin & period<=pmax;
  score = 100 * sqrt(sum(Pk(innyq))/L);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bw = keep_largest_cc(bw)
  cc = bwconncomp(bw);
  if cc.NumObjects<=1, return; end
  n = cellfun(@numel,cc.PixelIdxList);
  [~,mi] = max(n);
  bw = false(size(bw));
  bw(cc.PixelIdxList{mi}) = true;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ext = coherence_extent(slice2d,maskslice,span,min_area,pmin,pmax,cthresh)

  % fraction of in-brain columns whose Nyquist-band ripple is in-phase with the
  % row-mean ripple. Real full-width banding is coherent across columns (-> ~1);
  % row-to-row noise has random per-column phase (-> ~0).
  ext = 0;

  imgo = local_normalize(rot90(squeeze(slice2d)));
  if max(imgo(:))<=0, return; end

  % same brain mask as score_slice (supplied mask, else intensity fallback)
  if isempty(maskslice)
    bw = medfilt2(imgo,[3 3]) > 0.15*max(imgo(:));
    bw = imerode(bw,strel('disk',2));
    if ~any(bw(:)), return; end
    bw = keep_largest_cc(bw);
    bw = imfill(bw,'holes');
  else
    bw = rot90(squeeze(maskslice)) > 0.5;
    if ~any(bw(:)), return; end
  end
  if mean(bw(:)) < min_area, return; end

  brainw = sum(bw,2);
  idx = find(brainw > 0.25*max(brainw));
  if numel(idx) < 3*span, return; end

  m = imgo; m(~bw) = NaN;

  % reference ripple = Nyquist-bandpass of the (dome-removed) row-mean profile
  seg = nanmean(m,2); seg = seg(idx);
  med = median(seg,'omitnan');
  if med<=0 || isnan(med), return; end
  ref = nyq_bandpass((seg - movmean(seg,span,'omitnan'))/med,pmin,pmax);
  if norm(ref)<=0, return; end

  % per in-brain column: correlate its Nyquist ripple with the reference
  cols = find(sum(bw(idx,:),1) > 0.25*numel(idx));
  nc = 0; ncoh = 0;
  for c=cols(:)'
    cseg = m(idx,c);
    if sum(~isnan(cseg)) < 3*span, continue; end
    cs = (cseg - movmean(cseg,span,'omitnan'))/med; cs(isnan(cs)) = 0;
    rc = nyq_bandpass(cs,pmin,pmax);
    if norm(rc)<=0, continue; end
    nc = nc + 1;
    if corr_coef(ref,rc) >= cthresh, ncoh = ncoh + 1; end
  end
  if nc>0, ext = ncoh/nc; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = nyq_bandpass(x,pmin,pmax)
  % keep only Nyquist-band frequencies (period pmin-pmax rows), inverse-transform
  x = x(:); L = numel(x);
  X = fft(x - mean(x));
  k = (0:L-1)';
  period = L ./ max(k,0.5);
  keep = (period>=pmin & period<=pmax);
  keep = keep | flipud(keep);   % keep conjugate-symmetric partner
  X(~keep) = 0;
  y = real(ifft(X));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = corr_coef(a,b)
  a = a(:) - mean(a); b = b(:) - mean(b);
  d = norm(a)*norm(b);
  if d<=0, c = 0; else c = (a'*b)/d; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slice_scores = all_slice_scores(vol,brainmask,span,min_area,pmin,pmax)
  % per dim1 slice Nyquist-band score; NaN for non-brain slices
  dslices = size(vol,1);
  slice_scores = nan(dslices,1);
  for s=1:dslices
    if isempty(brainmask), ms=[]; else ms=squeeze(brainmask(s,:,:)); end
    slice_scores(s) = score_slice(squeeze(vol(s,:,:)),ms,span,min_area,pmin,pmax);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rows = banded_rows(slice2d,maskslice,span,min_area,pmin,pmax,topN)

  % row positions (in the rot90 plot image) of the strongest Nyquist-band ripple
  % peaks on this slice, for overlaying the detected line locations. [] if the
  % slice is not scorable. Mirrors score_slice's preprocessing.
  rows = [];
  imgo = local_normalize(rot90(squeeze(slice2d)));
  if max(imgo(:))<=0, return; end
  if isempty(maskslice)
    bw = medfilt2(imgo,[3 3]) > 0.15*max(imgo(:));
    bw = imerode(bw,strel('disk',2));
    if ~any(bw(:)), return; end
    bw = keep_largest_cc(bw);
    bw = imfill(bw,'holes');
  else
    bw = rot90(squeeze(maskslice)) > 0.5;
    if ~any(bw(:)), return; end
  end
  if mean(bw(:)) < min_area, return; end

  m = imgo; m(~bw) = NaN;
  r = nanmean(m,2);
  brainw = sum(bw,2);
  idx = find(brainw > 0.25*max(brainw));
  if numel(idx) < 3*span, return; end
  seg = r(idx);
  med = median(seg,'omitnan');
  if med<=0 || isnan(med), return; end
  rd = (seg - movmean(seg,span,'omitnan')) / med; rd(isnan(rd)) = 0;

  % peaks of the Nyquist-bandpassed ripple above half its max amplitude
  ripple = nyq_bandpass(rd,pmin,pmax);
  a = abs(ripple);
  if max(a)<=0, return; end
  [pk,loc] = findpeaks(a);
  keep = pk > 0.5*max(a);
  pk = pk(keep); loc = loc(keep);
  [~,ord] = sort(pk,'descend'); loc = loc(ord);
  if numel(loc) > topN, loc = loc(1:topN); end
  rows = idx(loc);   % map back to absolute rot90-image rows
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = local_normalize(x)
  x = double(x);
  y = x - min(x(:));
  mx = max(y(:));
  if mx>0, y = y/mx; end
return;
