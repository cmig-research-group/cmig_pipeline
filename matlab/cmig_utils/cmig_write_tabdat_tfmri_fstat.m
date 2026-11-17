function cmig_write_tabdat_tfmri_fstat(varargin)
%function cmig_write_tabdat_tfmri_fstat([options])
%
% Purpose: write tabulated F-statistics from tFMRI deconvolution
%
% Optional Input:
%   'indir': input directory containing roi_deconv mat file
%     {default = pwd}
%   'instem': input file stem for roi_deconv mat file (without run_tag prefix)
%     {default = 'roi_deconv_3dDeconv'}
%   'infostem': input file stem for info file (for future multi-subject use)
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': base output file stem (task/run_tag/atlas appended automatically)
%     {default = 'mr_y_tfmri'}
%   'colstem': base column label stem (task/run_tag/atlas appended automatically)
%     {default = 'mr_y_tfmri'}
%   'task': task name embedded in column and file names
%     {default = ''}
%   'run_tag': per-run tag in mat filename (e.g. 'run-01'); drives col_run_tag
%   'col_run_tag': per-run tag for the column names (e.g. 'r01'), when it
%     differs from the one run_tag gives
%     {default = ''}
%     {default = ''}
%   'tabdir': directory containing ROI abbreviation lookup table CSV files
%     {default = [getenv('MMPS_PARMS') '/tabulated']}
%   'participant_id': participant ID; if empty, multi-subject mode (not yet implemented)
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  05/20/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 10/03/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'instem','roi_deconv_3dDeconv',[],...
  'infostem','vol_info',[],...
  'outdir',pwd,[],...
  'outstem','mr_y_tfmri',[],...
  'colstem','mr_y_tfmri',[],...
  'task','',[],...
  'run_tag','',[],...
  'col_run_tag','',[],...
  'tabdir',[],[],...
  'participant_id',[],[],...
  'session_id',[],[],...
...
  'atlaslist',{'aseg','dsk','dst'},[],...
  'ext','.tsv',[],...
...
  'forceflag',false,[false true],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

% build column run tag from run_tag (e.g. 'run-01' -> 'r-01'; 'run01' -> 'r01')
if ~isempty(parms.run_tag)
  col_run_tag = strrep(parms.run_tag,'run','r');
  instem_full = sprintf('%s_%s',parms.run_tag,parms.instem);
else
  col_run_tag = '';
  instem_full = parms.instem;
end
% the columns may carry the release slot while the file names keep the run
if ~isempty(parms.col_run_tag), col_run_tag = parms.col_run_tag; end

% check if any output file is needed before loading data
any_needed = false;
for a = 1:length(parms.atlaslist)
  atlas = parms.atlaslist{a};
  fname_chk = build_bids_fname_fstat(parms,col_run_tag,atlas);
  if ~exist(fname_chk,'file') || parms.forceflag
    any_needed = true; break;
  end
end
if ~any_needed, return; end

fname_mat = sprintf('%s/%s.mat',parms.indir,instem_full);
if ~exist(fname_mat,'file')
  fprintf('%s: WARNING: deconv mat not found: %s\n',mfilename,fname_mat);
  return;
end

data = load(fname_mat,'vol','roinames');
fstat = data.vol(:,1);  % overall model F-stat [nroi, 1]

[roi_names, roi_abbrevs, roi_types] = cmig_load_abbrev('indir',parms.tabdir);
if isempty(roi_names), return; end

participant_id = parms.participant_id;
session_id     = parms.session_id;

mmil_mkdir(parms.outdir);

for a = 1:length(parms.atlaslist)
  atlas = parms.atlaslist{a};
  colstem  = build_colstem(parms.colstem,parms.task,col_run_tag,atlas);
  fname_out = build_bids_fname_fstat(parms,col_run_tag,atlas);
  if exist(fname_out,'file') && ~parms.forceflag, continue; end

  % select ROIs matching this atlas type
  idx_type = find(strcmp(atlas,roi_types));
  if isempty(idx_type), continue; end
  tmp_roinames = roi_names(idx_type);
  tmp_abbrevs  = roi_abbrevs(idx_type);
  [~,idx_type_rois,idx_rois] = intersect(tmp_roinames,data.roinames,'stable');
  if isempty(idx_rois), continue; end

  vals = fstat(idx_rois)';
  colnames = cell(1,length(idx_rois));
  for i = 1:length(idx_rois)
    colnames{i} = sprintf('%s__%s_fstat',colstem,tmp_abbrevs{idx_type_rois(i)});
  end

  prepvals = num2cell(vals);
  prepvals(isnan(vals)) = {'n/a'};
  header   = cat(2,{'participant_id','session_id'},colnames);
  data_row = cat(2,{participant_id,session_id},prepvals);
  mmil_write_csv(fname_out,cat(1,header,data_row),'separator','\t','quote_flag',0);
  fprintf('%s: wrote %s\n',mfilename,fname_out);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colstem = build_colstem(base,task,col_run_tag,atlas)
  if ~isempty(col_run_tag)
    colstem = sprintf('%s__%s__fstat__%s__%s',base,task,col_run_tag,atlas);
  elseif ~isempty(task)
    colstem = sprintf('%s__%s__fstat__%s',base,task,atlas);
  else
    colstem = sprintf('%s__fstat__%s',base,atlas);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outstem = build_outstem(base,task,col_run_tag,atlas)
  if ~isempty(col_run_tag)
    outstem = sprintf('%s__%s__fstat__%s__%s',base,task,col_run_tag,atlas);
  elseif ~isempty(task)
    outstem = sprintf('%s__%s__fstat__%s',base,task,atlas);
  else
    outstem = sprintf('%s__fstat__%s',base,atlas);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = build_bids_fname_fstat(parms,col_run_tag,atlas)
  if ~isempty(parms.participant_id) && ~isempty(parms.task)
    base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
    if ~isempty(parms.run_tag)
      run_num = strrep(parms.run_tag,'run','');
      base = sprintf('%s_run-%s',base,run_num);
    end
    fname = sprintf('%s/%s_atlas-%s_fstat%s',parms.outdir,base,atlas,parms.ext);
  else
    outstem = build_outstem(parms.outstem,parms.task,col_run_tag,atlas);
    fname = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

