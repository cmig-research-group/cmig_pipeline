function cmig_write_tabdat_tfmri(varargin)
%function cmig_write_tabdat_tfmri([options])
%
% Optional Input:
%   'indir': input directory containing roi_deconv mat file
%     {default = pwd}
%   'participant_id': participant ID; if empty, multi-subject mode
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'infostem': input file stem for info file (used in multi-subject mode)
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mr_y_tfmri'}
%   'colstem': prefix attached to column names
%     {default = 'mr_y_tfmri'}
%   'atlas': atlas name / ROI type
%     'dsk' = Desikan aparc
%     'dst' = Destrieux aparc
%     'aseg' = automated brain segmentation
%     {default = 'dsk'}
%   'fs_version': FreeSurfer version number; e.g., 450, 530, 711, 810
%     {default = 711}
%   'instem': input file stem for roi_deconv mat file
%     {default = 'roi_deconv_3dDeconv'}
%   'condidx': condition index (1-based) into deconv mat columns
%     {default = 0}
%   'condname': condition name for output file stem
%     {default = ''}
%   'task': task name for output file stem
%     {default = ''}
%   'run_tag': per-run tag (e.g. 'r01') for column names
%     {default = ''}
%   'fname_run_tag': per-run tag for the file name's run entity (e.g.
%     'run03' -> run-03), when it differs from the columns'
%     if empty, run_tag
%     {default = ''}
%   'tabdir': directory containing ROI abbreviation lookup table CSV files
%     {default = [getenv('MMPS_PARMS') '/tabulated']}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  05/15/2026 by Don Hagler
% Prev Mod: 05/20/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 09/16/2026 by Don Hagler
% Prev Mod: 09/17/2026 by Don Hagler
% Prev Mod: 10/03/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'participant_id',[],[],...
  'session_id',[],[],...
  'infostem','vol_info',[],...
...
  'outdir',pwd,[],...
  'outstem','mr_y_tfmri',[],...
  'colstem','mr_y_tfmri',[],...
  'atlas','dsk',{'dsk','dst','aseg'},...
  'fs_version',711,[],...
...
  'aparc_roilist',[1001:1003,1005:1035,2001:2003,2005:2035],[],...
  'aparc2_roilist',[11101:11141,11143:11175,12101:12141,12143:12175],[],...
  'aseg_roilist',[2,4:5,7,8,10:18,24,26,28,41,43:44,46,47,49:54,58,60,77,251:255],[1,Inf],...
...
  'tabdir',[],[],...
...
  'ext','.tsv',[],...
...
  'instem','roi_deconv_3dDeconv',[],...
  'condidx',0,[],...
  'condname','',[],...
  'task','',[],...
  'run_tag','',[],...
  'fname_run_tag','',[],...
...
  'forceflag',false,[false true],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

% the file name may carry the BIDS run while the columns carry the slot
fname_run_tag = parms.fname_run_tag;
if isempty(fname_run_tag), fname_run_tag = parms.run_tag; end
if ~isempty(parms.participant_id) && ~isempty(parms.task)
  fname_out = build_bids_fname_tfmri(parms.outdir,parms.participant_id,parms.session_id,...
    parms.task,parms.condname,fname_run_tag,parms.atlas,parms.ext);
else
  outstem = build_outstem(parms.outstem,parms.task,parms.condname,parms.run_tag,parms.atlas);
  fname_out = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
end
if exist(fname_out,'file') && ~parms.forceflag, return; end

parms = get_roi_info(parms);
colstem = build_outstem(parms.colstem,parms.task,parms.condname,parms.run_tag,parms.atlas);
mmil_mkdir(parms.outdir);

if isempty(parms.participant_id)
  % multi-subject mode: load participant IDs from info file
  fname_info = sprintf('%s/%s.mat',parms.indir,parms.infostem);
  if ~exist(fname_info,'file')
    fprintf('%s: WARNING: info file not found: %s\n',mfilename,fname_info);
    return;
  end
  vol_info = load(fname_info);
  participant_ids = vol_info.participant_id;
  session_ids     = vol_info.session_id;
  nsubs = length(participant_ids);
  fname_mat = sprintf('%s/%s.mat',parms.indir,parms.instem);
  if ~exist(fname_mat,'file')
    fprintf('%s: WARNING: deconv mat not found: %s\n',mfilename,fname_mat);
    return;
  end
  tmp = load(fname_mat,'vol','roinames');
  all_beta = NaN(nsubs,length(parms.roinames));
  all_sem  = NaN(nsubs,length(parms.roinames));
  for i = 1:nsubs
    [b,s] = extract_deconv(tmp,parms,i);
    if ~isempty(b)
      all_beta(i,:) = b;
      all_sem(i,:)  = s;
    end
  end
  write_tsv_multi(fname_out,all_beta,all_sem,parms.shortnames,colstem,...
    participant_ids,session_ids);
else
  % single-subject mode
  [beta,sem] = load_deconv_data(parms);
  if isempty(beta), return; end
  write_tsv_file(fname_out,beta,sem,parms.shortnames,colstem,parms);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = get_roi_info(parms)
  [roi_names, roi_abbrevs] = cmig_load_abbrev('indir', parms.tabdir);
  switch parms.atlas
    case {'dsk','dst'}
      parms = set_aparc_rois(parms);
    case 'aseg'
      parms = set_aseg_rois(parms);
  end
  parms.shortnames = cell(size(parms.roinames));
  for i = 1:length(parms.roinames)
    idx_roi = find(strcmp(parms.roinames{i},roi_names));
    if isempty(idx_roi)
      fprintf('%s: WARNING: no ROI in look-up table called %s\n',mfilename,parms.roinames{i});
      parms.shortnames{i} = parms.roinames{i};
    else
      parms.shortnames{i} = roi_abbrevs{idx_roi};
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_aparc_rois(parms)
  [all_roicodes,all_roinames] = mmil_colorlut('FS_version',parms.fs_version);
  switch parms.atlas
    case 'dsk'
      roilist = parms.aparc_roilist;
    case 'dst'
      roilist = parms.aparc2_roilist;
  end
  [~,ind_rois] = intersect(all_roicodes,roilist);
  parms.roicodes = mmil_rowvec(all_roicodes(ind_rois));
  parms.roinames = mmil_rowvec(all_roinames(ind_rois));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = set_aseg_rois(parms)
  [all_roicodes,all_roinames] = mmil_colorlut('FS_version',parms.fs_version);
  [~,ind_rois] = intersect(all_roicodes,parms.aseg_roilist);
  parms.roicodes = mmil_rowvec(all_roicodes(ind_rois));
  parms.roinames = mmil_rowvec(all_roinames(ind_rois));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stem = build_outstem(base,task,condname,run_tag,atlas)
  if ~isempty(run_tag)
    stem = sprintf('%s__%s__%s__%s__%s',base,task,condname,run_tag,atlas);
  else
    stem = sprintf('%s__%s__%s__%s',base,task,condname,atlas);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = build_bids_fname_tfmri(outdir,pid,sid,task,condname,run_tag,atlas,ext)
  base = sprintf('%s_%s_task-%s',pid,sid,task);
  if ~isempty(run_tag)
    run_num = strrep(run_tag,'run','');
    base = sprintf('%s_run-%s',base,run_num);
  end
  base = sprintf('%s_atlas-%s',base,atlas);
  if ~isempty(condname)
    base = sprintf('%s_desc-%s',base,condname);
  end
  fname = sprintf('%s/%s_deconv%s',outdir,base,ext);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta,sem] = load_deconv_data(parms)
  beta = []; sem = [];
  fname_mat = sprintf('%s/%s.mat',parms.indir,parms.instem);
  if ~exist(fname_mat,'file')
    fprintf('%s: WARNING: deconv mat not found: %s\n',mfilename,fname_mat);
    return;
  end
  tmp = load(fname_mat,'vol','roinames');
  [beta,sem] = extract_deconv(tmp,parms,[]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta,sem] = extract_deconv(tmp,parms,subj_idx)
  beta = []; sem = [];
  vol          = tmp.vol;
  roinames_all = tmp.roinames;
  c = parms.condidx;
  if ~isempty(subj_idx)
    % multi-subject: vol is [nrois x ncols x nsubj]
    if ndims(vol) < 3 || size(vol,3) < subj_idx
      return;
    end
    v = vol(:,:,subj_idx);
  else
    v = vol;
  end
  if size(v,2) < 2*c+1
    fprintf('%s: WARNING: deconv mat has %d cols, need %d for condidx %d\n',...
      mfilename,size(v,2),2*c+1,c);
    return;
  end
  % a condition with no event in any run was never measured; its all-zero
  %   regressor returns a beta of zero, which would publish as a value, so
  %   leave it NaN and let the writer mark it n/a
  fname_condnames = sprintf('%s/deconv_condnames.mat',parms.indir);
  empty_conds = cmig_get_empty_conds(fname_condnames,parms.run_tag);
  if c <= length(empty_conds) && empty_conds(c)
    beta = NaN(1,length(parms.roinames));
    sem  = NaN(1,length(parms.roinames));
    return;
  end
  beta_all = v(:,2*c);
  sem_all  = v(:,2*c+1);
  [~,idx_parms,idx_all] = intersect(parms.roinames,roinames_all,'stable');
  beta = NaN(1,length(parms.roinames));
  sem  = NaN(1,length(parms.roinames));
  beta(1,idx_parms) = beta_all(idx_all)';
  sem(1,idx_parms)  = sem_all(idx_all)';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tsv_file(fname_out,beta,sem,shortnames,colstem,parms)
  idx_lh  = find(contains(shortnames,'__lh'));
  idx_rh  = find(contains(shortnames,'__rh'));
  idx_mid = find(~contains(shortnames,'__lh') & ~contains(shortnames,'__rh'));
  idx_ord = [idx_lh, idx_rh, idx_mid];
  names_ord = shortnames(idx_ord);
  beta_ord  = beta(idx_ord);
  sem_ord   = sem(idx_ord);
  colnames  = [strcat(colstem,'__',names_ord,'_beta'), strcat(colstem,'__',names_ord,'_sem')];
  vals      = [beta_ord, sem_ord];
  prepvals  = num2cell(vals);
  prepvals(isnan(vals)) = {'n/a'};
  header   = cat(2,{'participant_id','session_id'},colnames);
  data_row = cat(2,{parms.participant_id,parms.session_id},prepvals);
  mmil_write_csv(fname_out,cat(1,header,data_row),'separator','\t','quote_flag',0);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tsv_multi(fname_out,all_beta,all_sem,shortnames,colstem,...
    participant_ids,session_ids)
  nsubs = length(participant_ids);
  idx_lh  = find(contains(shortnames,'__lh'));
  idx_rh  = find(contains(shortnames,'__rh'));
  idx_mid = find(~contains(shortnames,'__lh') & ~contains(shortnames,'__rh'));
  idx_ord = [idx_lh, idx_rh, idx_mid];
  names_ord = shortnames(idx_ord);
  colnames = [strcat(colstem,'__',names_ord,'_beta'), strcat(colstem,'__',names_ord,'_sem')];
  header = cat(2,{'participant_id','session_id'},colnames);
  rows = cell(nsubs,1);
  for i = 1:nsubs
    beta_ord = all_beta(i,idx_ord);
    sem_ord  = all_sem(i,idx_ord);
    vals = [beta_ord, sem_ord];
    pv = num2cell(vals);
    pv(isnan(vals)) = {'n/a'};
    rows{i} = cat(2,{participant_ids{i},session_ids{i}},pv);
  end
  data_mat = cat(1,header,rows{:});
  mmil_write_csv(fname_out,data_mat,'separator','\t','quote_flag',0);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

