function [outdir,errcode,errmsg] = cmig_concat_tfmri_stats(participant_id,session_id,output_rootdir,task,varargin)
%function [outdir,errcode,errmsg] = cmig_concat_tfmri_stats(participant_id,session_id,output_rootdir,task,[options])
%
% Purpose: Save per-subject task fMRI deconvolution stats in BIDS .mat format
%   Reads the ico surface deconvolution bucket produced by
%   cmig_analyze_bids_func (run_deconv_ico).
%   Outputs one .mat per condition x hemisphere, for beta and sem, and for
%   amp and lat when iresp_flag is set.
%
% Required Input:
%   participant_id: participant identifier string (e.g. 'sub-XXXXXXXX')
%   session_id: session identifier string (e.g. 'ses-00A')
%   output_rootdir: root directory for analysis outputs
%   task: task name (e.g. 'nback'), naming the analyze_func_{task} subdir
%
% Optional Parameters ('key',value pairs):
%   'fname_info': full path of the vol_info mat file written by
%     cmig_save_vol_info, supplying the cortex vertex masks
%     {default = []}
%   'nruns': number of runs in the deconvolution
%     {default = 1}
%   'hemilist': cell array of hemisphere labels
%     {default = {'lh','rh'}}
%   'ico_order': icosahedral order of the deconvolution output
%     {default = 5}
%   'run_tag': prefix identifying a per-run deconvolution; empty for the
%     session-level one
%     {default = []}
%   'iresp_flag': [0|1] also save the impulse response amp and lat outputs
%     {default = 0}
%   'baseline_flag': [0|1] the deconvolution kept the raw regressor
%     coefficients, so the bucket carries columns ahead of the GLT block, and
%     the beta and sem are post-normalized by the baseline,
%     {pid}_{sid}_task-{task}[_run-{r}]_hemi-{L|R}_atlas-ico5_baseline.mat in
%     the analyze_func_{task} subdir
%     {default = 0}
%   'mask_flag': [0|1] set the beta, sem, amp and lat to NaN at vertices outside the
%     coverage mask, {pid}_{sid}_task-{task}[_run-{r}]_hemi-{L|R}_atlas-ico5_mask.mat
%     in outdir (measmat [1, nverts], 1 covered, 0 not)
%     {default = 0}
%   'outliers_flag': [0|1] apply the per-condition censoring recorded in
%     deconv_outliers.mat, leaving those conditions all-NaN
%     {default = 1}
%   'stim_fnames': cell array of stim timing file names, one per condition;
%     needed for the bucket column offset when baseline_flag is set
%     {default = []}
%   'nreg_per_stim': regressors per condition; needed for the offset
%     {default = 1}
%   'detrend': polort order used in the deconvolution; needed for the offset
%     {default = 3}
%   'fname_motion': motion regressor file name; needed for the offset
%     {default = []}
%   'motion_labels': cell array of motion regressor labels
%     {default = {'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'}}
%   'deriv_flag': [0|1] motion derivatives were included
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output files {default = 0}
%
% Output:
%   outdir: output directory
%   errcode: 0 if successful, 1 if error
%   errmsg: error message string if errcode=1, otherwise []
%
% Input naming (from run_deconv_ico, in {outdir}/analyze_func_{task}):
%   {run_tag_}ico{N}_deconv_{hemi}_3dDeconv.mgh
%   {run_tag_}ico{N}_deconv_{hemi}_3dDeconv_stats_info.mat
%   deconv_condnames.mat, {run_tag_}deconv_outliers.mat
%   {pid}_{sid}_task-{task}[_run-{r}]_hemi-{L|R}_atlas-ico5_baseline.mat
%
% Output BIDS .mat naming:
%   {outdir}/{pid}_{sid}_task-{task}[_run-{r}]_hemi-{L|R}_atlas-ico5_desc-{cond}_{suffix}.mat
%   suffix: beta, sem, and amp and lat when iresp_flag
%
% NOTE: unlike cmig_concat_smri_measures and cmig_concat_dwi_measures, each
%   output carries dofvec alongside measmat; downstream readers use the
%   degrees of freedom to weight across runs
%
% NOTE: ico_order selects the input deconvolution bucket, but the atlas entity
%   in the output names is the literal ico5, as it is throughout
%   cmig_analyze_bids_func, whose check_output looks for that literal. Only
%   matters if ico_order is ever changed from 5, and it would need changing in
%   both places together.
%
% Created:  09/07/2026 by Don Hagler
% Prev Mod: 09/07/2026 by Don Hagler
% Prev Mod: 09/08/2026 by Don Hagler
% Prev Mod: 09/16/2026 by Don Hagler
% Prev Mod: 09/17/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,4), errcode = 1; return; end

parms = mmil_args2parms(varargin,{...
  'participant_id',participant_id,[],...
  'session_id',session_id,[],...
  'output_rootdir',output_rootdir,[],...
  'task',task,[],...
  'fname_info',[],[],...
  'nruns',1,[1,Inf],...
  'hemilist',{'lh','rh'},[],...
  'ico_order',5,[0 7],...
  'run_tag',[],[],...
  'iresp_flag',false,[false true],...
  'baseline_flag',false,[false true],...
  'mask_flag',false,[false true],...
  'outliers_flag',true,[false true],...
  'stim_fnames',[],[],...
  'nreg_per_stim',1,[1,Inf],...
  'detrend',3,[0:3],...
  'fname_motion',[],[],...
  'motion_labels',{'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'},[],...
  'deriv_flag',true,[false true],...
  'forceflag',false,[false true],...
...
  'modality','func',[],...
});

outdir = sprintf('%s/%s/%s/%s',...
  parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);
tmpdir = sprintf('%s/analyze_%s_%s',outdir,parms.modality,parms.task);
mmil_mkdir(outdir);

% get condnames (conditions + contrasts) from deconvolution output
fname_condnames = sprintf('%s/deconv_condnames.mat',tmpdir);
if ~exist(fname_condnames,'file')
  errcode = 1;
  errmsg = sprintf('condnames file not found: %s',fname_condnames);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  return;
end
tmp = load(fname_condnames,'condnames');
condnames = tmp.condnames;
% conditions with no event, written by prep_stims; per run when a run_tag
%   names one, otherwise absent from every run
empty_conds = cmig_get_empty_conds(fname_condnames,parms.run_tag);
if length(empty_conds) ~= length(condnames)
  empty_conds = false(1,length(condnames));
end
if isempty(condnames)
  errcode = 1;
  errmsg = sprintf('no condition names found in %s',fname_condnames);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  return;
end
npars = length(condnames);

% baseline_flag drops -nocout (see mmil_run_3dDeconv), so the bucket also
% contains raw condition and motion/deriv coefficients before the GLT block
% that condnames indexes into; cmig_deconv_col_offset counts them
col_offset = cmig_deconv_col_offset(parms,parms.nruns);

% load outlier censoring determined from roi sem rms (see postnorm_deconv_roi);
% propagates the same per-condition/contrast censoring to ico outputs
excluded = false(1,npars);
if parms.outliers_flag
  if ~isempty(parms.run_tag)
    fname_outliers = sprintf('%s/%s_deconv_outliers.mat',tmpdir,parms.run_tag);
  else
    fname_outliers = sprintf('%s/deconv_outliers.mat',tmpdir);
  end
  if exist(fname_outliers,'file')
    tmp_out = load(fname_outliers,'excluded');
    if length(tmp_out.excluded) == npars
      excluded = tmp_out.excluded;
    end
  end
end

bids_base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
if ~isempty(parms.run_tag)
  run_entity = sprintf('_run-%s',strrep(parms.run_tag,'run',''));
else
  run_entity = '';
end

% load vol_info
if isempty(parms.fname_info) || ~exist(parms.fname_info,'file')
  errcode = 1;
  errmsg = sprintf('vol_info file not found: %s',parms.fname_info);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  return;
end
vol_info = load(parms.fname_info);

fprintf('%s: saving tfmri stats in concatenated format...\n',mfilename);
for hemii = 1:length(parms.hemilist)
  hemi = parms.hemilist{hemii};
  switch hemi
    case 'lh'
      nverts = length(vol_info.vol_mask_lh);
      verts  = find(vol_info.vol_mask_lh);
    case 'rh'
      nverts = length(vol_info.vol_mask_rh);
      verts  = find(vol_info.vol_mask_rh);
  end

  if ~isempty(parms.run_tag)
    outstem_deconv = sprintf('%s_ico%d_deconv_%s',parms.run_tag,parms.ico_order,hemi);
  else
    outstem_deconv = sprintf('ico%d_deconv_%s',parms.ico_order,hemi);
  end
  fname_info_deconv = sprintf('%s/%s_3dDeconv_stats_info.mat',tmpdir,outstem_deconv);

  % load DOF from stats_info file
  dofvec = 0;
  if ~exist(fname_info_deconv,'file')
    errcode = 1;
    errmsg = sprintf('deconv stats info file not found: %s',fname_info_deconv);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  tmp_info = load(fname_info_deconv,'stats_info');
  si = tmp_info.stats_info;
  for si_i = 1:length(si)
    if strcmp(si(si_i).type,'tstat') && ~isempty(si(si_i).dofs)
      dofvec = si(si_i).dofs;
      break;
    end
  end

  % load beta/sem deconv output; reshape to [nverts,nframes]
  fname_deconv_beta = sprintf('%s/%s_3dDeconv.mgh',tmpdir,outstem_deconv);
  if ~exist(fname_deconv_beta,'file')
    errcode = 1;
    errmsg = sprintf('deconv output not found for %s: %s',hemi,fname_deconv_beta);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  deconv_beta = fs_load_mgh(fname_deconv_beta);
  deconv_beta = reshape(deconv_beta,[nverts,size(deconv_beta,4)]);

  % load amp/lat deconv output if iresp_flag
  if parms.iresp_flag
    fname_deconv_amp = sprintf('%s/%s_3dDeconv_iresp_amp.mgh',tmpdir,outstem_deconv);
    fname_deconv_lat = sprintf('%s/%s_3dDeconv_iresp_lat.mgh',tmpdir,outstem_deconv);
    if ~exist(fname_deconv_amp,'file')
      errcode = 1;
      errmsg = sprintf('iresp amp output not found for %s: %s',hemi,fname_deconv_amp);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    if ~exist(fname_deconv_lat,'file')
      errcode = 1;
      errmsg = sprintf('iresp lat output not found for %s: %s',hemi,fname_deconv_lat);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    deconv_amp = fs_load_mgh(fname_deconv_amp);
    deconv_amp = reshape(deconv_amp,[nverts,size(deconv_amp,4)]);
    deconv_lat = fs_load_mgh(fname_deconv_lat);
    deconv_lat = reshape(deconv_lat,[nverts,size(deconv_lat,4)]);
  end

  % load baseline for post-normalization (per hemisphere)
  baseline_hemi = [];
  if parms.baseline_flag
    bids_hemi_bl = cmig_get_bids_hemi(hemi);
    fname_baseline_hemi = sprintf('%s/%s%s_hemi-%s_atlas-ico5_baseline.mat',...
      tmpdir,bids_base,run_entity,bids_hemi_bl);
    if exist(fname_baseline_hemi,'file')
      tmp_bl = load(fname_baseline_hemi,'measmat');
      baseline_hemi = tmp_bl.measmat;  % [1, nverts]
    else
      fprintf('%s: WARNING: baseline file not found, %s stats not post-normalized: %s\n',...
        mfilename,hemi,fname_baseline_hemi);
    end
  end

  % load coverage mask (per hemisphere); vertices outside it get NaN
  mask_hemi = [];
  if parms.mask_flag
    fname_mask_hemi = sprintf('%s/%s%s_hemi-%s_atlas-ico5_mask.mat',...
      outdir,bids_base,run_entity,cmig_get_bids_hemi(hemi));
    if ~exist(fname_mask_hemi,'file')
      errcode = 1;
      errmsg = sprintf('coverage mask file not found: %s',fname_mask_hemi);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    tmp_mask = load(fname_mask_hemi,'measmat');
    mask_hemi = tmp_mask.measmat > 0;  % [1, nverts]
    if numel(mask_hemi) ~= nverts
      errcode = 1;
      errmsg = sprintf('coverage mask %s has %d vertices, not %d',...
        fname_mask_hemi,numel(mask_hemi),nverts);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    mask_hemi = reshape(mask_hemi,[1,nverts]);
  end

  for c = 1:npars
    condname = condnames{c};
    bids_hemi = cmig_get_bids_hemi(hemi);
    fname_beta = sprintf('%s/%s%s_hemi-%s_atlas-ico5_desc-%s_beta.mat',outdir,bids_base,run_entity,bids_hemi,condname);
    fname_sem  = sprintf('%s/%s%s_hemi-%s_atlas-ico5_desc-%s_sem.mat', outdir,bids_base,run_entity,bids_hemi,condname);
    if parms.iresp_flag
      fname_amp = sprintf('%s/%s%s_hemi-%s_atlas-ico5_desc-%s_amp.mat',outdir,bids_base,run_entity,bids_hemi,condname);
      fname_lat = sprintf('%s/%s%s_hemi-%s_atlas-ico5_desc-%s_lat.mat',outdir,bids_base,run_entity,bids_hemi,condname);
      all_exist = exist(fname_beta,'file') && exist(fname_sem,'file') && ...
                  exist(fname_amp,'file') && exist(fname_lat,'file');
    else
      all_exist = exist(fname_beta,'file') && exist(fname_sem,'file');
    end
    if all_exist && ~parms.forceflag, continue; end

    % stats file: col 1=F-stat, then polort section (if baseline_flag), then interleaved [beta,sem] per condition
    if size(deconv_beta,2) < 1+col_offset+npars*2, continue; end
    vals_beta = deconv_beta(:,1+col_offset+(c-1)*2+1);
    vals_sem  = deconv_beta(:,1+col_offset+(c-1)*2+2);
    % censor: leave measmat all-NaN for this condition. A condition with no
    %   event was never measured, so it takes the same path rather than
    %   publishing the zero its all-zero regressor returns
    is_outlier = excluded(c) || empty_conds(c);

    % save beta
    measmat = NaN(1,nverts,'single');
    if ~is_outlier
      measmat(verts) = vals_beta(verts);
      if ~isempty(baseline_hemi)
        measmat = 100 * measmat ./ (baseline_hemi + eps);
      end
      if ~isempty(mask_hemi), measmat(~mask_hemi) = NaN; end
    end
    save(fname_beta,'measmat','dofvec','-v7');
    fprintf('%s: wrote %s\n',mfilename,fname_beta);
    % save sem
    measmat = NaN(1,nverts,'single');
    if ~is_outlier
      measmat(verts) = vals_sem(verts);
      if ~isempty(baseline_hemi)
        measmat = 100 * measmat ./ (baseline_hemi + eps);
      end
      if ~isempty(mask_hemi), measmat(~mask_hemi) = NaN; end
    end
    save(fname_sem,'measmat','dofvec','-v7');
    fprintf('%s: wrote %s\n',mfilename,fname_sem);

    if parms.iresp_flag
      if size(deconv_amp,2) < c, continue; end
      vals_amp = deconv_amp(:,c);
      vals_lat = deconv_lat(:,c);
      % save amp
      measmat = NaN(1,nverts,'single');
      if ~is_outlier
        measmat(verts) = vals_amp(verts);
      end
      if ~isempty(mask_hemi), measmat(~mask_hemi) = NaN; end
      save(fname_amp,'measmat','dofvec','-v7');
      fprintf('%s: wrote %s\n',mfilename,fname_amp);
      % save lat
      measmat = NaN(1,nverts,'single');
      if ~is_outlier
        measmat(verts) = vals_lat(verts);
      end
      if ~isempty(mask_hemi), measmat(~mask_hemi) = NaN; end
      save(fname_lat,'measmat','dofvec','-v7');
      fprintf('%s: wrote %s\n',mfilename,fname_lat);
    end
  end
end

return;
