function [line_scores,line_counts,line_scores_fft,line_nslices_fft,line_extents_fft] = mmil_detect_lines(fname_in,varargin)
%function [line_scores,line_counts,line_scores_fft,line_nslices_fft,line_extents_fft] = mmil_detect_lines(fname_in,[options])
%
% Purpose: Calculates a measure of line artifacts
%
% Loads the input volume once and runs the original detector, plus - ONLY when the
% fft outputs are requested (nargout>2) - the redesigned fft detector:
%   autoqc_detect_lines     (original; line_*)   -> line_scores / line_counts
%   autoqc_detect_lines_fft (Nyquist-band; line2_*) -> line_scores_fft /
%                                                      line_nslices_fft / line_extents_fft
% With <=2 outputs the fft branch and its SynthStrip brain mask are skipped, so the
% line_* path is unchanged in cost/behavior (ABCD line_* stays byte-identical). The
% fft detector is the redesign that recovers the HBCD fmap false negatives.
%
% Usage:  [ls,lc,ls_fft,ns_fft,le_fft] = mmil_detect_lines(fname_in,'key1',value1,...);
%
% Required Input:
%   fname_in: full path of input volume file
%
% Optional Input:
%   'outdir': output directory
%     {default = pwd}
%   'plotflag': [0|1|2] plot images for lines detected
%     0: no plots
%     1: plot for one frame with highest line score
%     2: plots for all frames with line artifact detected
%     {default = 1}
%   'verbose': [0|1] display messagess
%     {default = 0}
%   'nslices': passed to the original detector (the fft detector always scans all
%     in-brain slices, ignoring this); if 0, analyze all. Default 1 = the single
%     mid-plane slice, matching production so line_* stays byte-identical (nslices=0
%     makes the original fire on - and plot - every out-of-brain noise slice)
%     {default = 1}
%   'frames': vector of frame numers to be analyzed
%     if empty, will analyze all frames
%     {default = []}
%   'fname_mask': precomputed SynthStrip brain mask file for the fft detector; if
%     empty, it is computed here with mmil_synthstrip_brainmask
%     {default = []}
%
% Output (one entry per analyzed frame):
%   line_scores / line_counts:                      original detector (always)
%   line_scores_fft / line_nslices_fft / line_extents_fft: Nyquist-band redesign
%     (computed only when nargout>2; [] otherwise)
%
% Created:  08/29/21 by Payam Saisan
% Prev Mod: 02/01/22 by Don Hagler
% Prev Mod: 07/12/2026 by Don Hagler
% Last Mod: 07/13/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

line_scores = []; line_counts = [];
line_scores_fft = []; line_nslices_fft = []; line_extents_fft = [];
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin,{...
  'outdir',[],[],...
  'plotflag',1,[0:2],...
  'verbose',false,[false,true],...
  'nslices',1,[],...
  'frames',[],[],...
  'fname_mask',[],[],...
});

[fpath,fstem] = fileparts(fname_in);

% run the fft/line2 branch only when its outputs are requested (nargout>2); a
% caller wanting only line_* ([scores,counts]) skips fft AND the SynthStrip mask
do_fft = (nargout > 2);

[vol,M] = fs_load_mgh(fname_in);   % load/uncompress once, reused by both detectors

% brain mask for the fft detector (it falls back to an intensity mask if empty).
% Reuse a caller-supplied fname_mask if given (MMIL_autoQC_Exam already builds
% one); else compute here with a unique name per input container dir so series
% sharing the BOLD1_*.mgz stem don't collide.
brainmask = [];
if do_fft
  if isempty(parms.outdir), maskdir = tempdir; else maskdir = parms.outdir; end
  [~,pdir] = fileparts(fpath);
  try
    fname_mask = parms.fname_mask;
    if isempty(fname_mask) || ~exist(fname_mask,'file')
      if isempty(fname_mask)
        fname_mask = fullfile(maskdir,sprintf('%s_%s_brainmask.mgz',pdir,fstem));
      end
      fname_mask = mmil_synthstrip_brainmask(fname_in,'fname_mask',fname_mask,...
        'binarize_flag',1,'verbose',parms.verbose);
    end
    bm = fs_load_mgh(fname_mask);
    brainmask = bm(:,:,:,1) > 0.5;
  catch ME
    fprintf('%s: WARNING: SynthStrip brain mask failed (%s); using intensity mask\n',...
      mfilename,ME.message);
  end
end

if isempty(parms.frames)
  parms.frames = 1:size(vol,4);
else
  parms.frames = intersect(parms.frames,1:size(vol,4));
end
nframes = length(parms.frames);

if parms.verbose
  if do_fft, dstr = 'orig + fft'; else dstr = 'orig'; end
  fprintf('%s: checking for lines in %d frames (%s)...\n',mfilename,nframes,dstr);
  tic
end
line_scores = zeros(nframes,1); line_counts = zeros(nframes,1);
if do_fft
  line_scores_fft  = zeros(nframes,1); line_nslices_fft = zeros(nframes,1);
  line_extents_fft = zeros(nframes,1);
end
for f=1:nframes
  fnum = parms.frames(f);
  V = vol(:,:,:,fnum);
  os = sprintf('%s_f%03d',fstem,fnum);
  [line_scores(f),line_counts(f)] = autoqc_detect_lines(...
    V,parms.nslices,parms.outdir,os,parms.verbose,parms.plotflag==2);
  if do_fft
    [line_scores_fft(f),line_nslices_fft(f),line_extents_fft(f)] = autoqc_detect_lines_fft(...
      V,parms.outdir,[os '_fft'],parms.verbose,parms.plotflag==2,brainmask);
  end
end

% plotflag==1: re-plot the single highest-scoring frame for each detector run
if parms.plotflag==1
  fns = {@(V,os) autoqc_detect_lines(V,parms.nslices,parms.outdir,os,parms.verbose,1)};
  scoresets = {line_scores};
  suffixes  = {''};
  if do_fft
    fns{end+1}       = @(V,os) autoqc_detect_lines_fft(V,parms.outdir,os,parms.verbose,1,brainmask);
    scoresets{end+1} = line_scores_fft;
    suffixes{end+1}  = '_fft';
  end
  for d=1:numel(fns)
    sc = scoresets{d};
    if any(sc>0)
      [~,fbest] = max(sc);
      fnum = parms.frames(fbest);
      os = sprintf('%s_f%03d%s',fstem,fnum,suffixes{d});
      fns{d}(vol(:,:,:,fnum),os);
    end
  end
end

if parms.verbose, toc; end

return;
