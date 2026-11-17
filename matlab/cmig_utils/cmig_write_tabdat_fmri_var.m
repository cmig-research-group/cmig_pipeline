function cmig_write_tabdat_fmri_var(varargin)
%function cmig_write_tabdat_fmri_var([options])
%
% Purpose: write tabulated ROI variance from fMRI data
%          (rsfmri or tfmri); one output file per ROI type
%
% Optional Input:
%   'indir': input directory containing roi_corrmat .mat file
%     {default = pwd}
%   'instem': input file stem for roi_corrmat .mat file (tag appended automatically)
%     {default = 'roi_corrmat'}
%   'infostem': input file stem for info file (used in multi-subject mode)
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': base output file stem (tag and roitype appended automatically)
%     {default = 'mr_y_fmri__var'}
%   'colstem': base column label stem (tag and roitype appended automatically)
%     {default = 'mr_y_fmri__var'}
%   'roitypes': ROI types to process; one output file per type
%     {default = {'aseg','dsk','dst','gp'}}
%   'tag': suffix appended to instem/outstem/colstem for tfmri variants
%     '' = direct signal; '__task' = task-predicted; '__resid' = residual
%     {default = ''}
%   'tabdir': directory containing ROI abbreviation lookup table CSV files
%     {default = [getenv('MMPS_PARMS') '/tabulated']}
%   'participant_id': participant ID; if empty, multi-subject mode
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  05/20/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Prev Mod: 09/29/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'instem','roi_corrmat',[],...
  'infostem','vol_info',[],...
  'outdir',pwd,[],...
  'outstem','mr_y_fmri__var',[],...
  'colstem','mr_y_fmri__var',[],...
  'task',[],[],...
  'desc',[],[],...
  'roitypes',{'aseg','dsk','dst','gp'},[],...
  'tag','',[],...
  'tabdir',[],[],...
  'participant_id',[],[],...
  'session_id',[],[],...
  'forceflag',false,[false true],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

full_instem = sprintf('%s%s',parms.instem,parms.tag);
fname_in = sprintf('%s/%s.mat',parms.indir,full_instem);
if ~exist(fname_in,'file')
  fprintf('%s: WARNING: input file not found: %s\n',mfilename,fname_in);
  return;
end

% check if any output file is needed
any_needed = false;
for t = 1:length(parms.roitypes)
  fname_chk = build_bids_fname_var(parms,parms.roitypes{t});
  if ~exist(fname_chk,'file') || parms.forceflag
    any_needed = true; break;
  end
end
if ~any_needed, return; end

if isempty(parms.participant_id)
  % multi-subject mode
  fname_info = sprintf('%s/%s.mat',parms.indir,parms.infostem);
  if ~exist(fname_info,'file')
    fprintf('%s: WARNING: info file not found: %s\n',mfilename,fname_info);
    return;
  end
  vol_info = load(fname_info);
  participant_ids = vol_info.participant_id;
  session_ids     = vol_info.session_id;
else
  participant_ids = {parms.participant_id};
  session_ids     = {parms.session_id};
end

fprintf('%s: loading %s...\n',mfilename,fname_in);
data = load(fname_in,'varmat','roinames');

[roi_names, roi_abbrevs, roi_types] = cmig_load_abbrev('indir',parms.tabdir);
if isempty(roi_names), return; end

nsubs = length(participant_ids);
mmil_mkdir(parms.outdir);

for t = 1:length(parms.roitypes)
  roitype = parms.roitypes{t};
  colstem_t = sprintf('%s%s__%s',parms.colstem,parms.tag,roitype);
  fname_out = build_bids_fname_var(parms,roitype);
  if exist(fname_out,'file') && ~parms.forceflag, continue; end

  % get ROIs for this type
  idx_type = find(strcmp(roitype,roi_types));
  if isempty(idx_type), continue; end
  tmp_roinames = roi_names(idx_type);
  tmp_abbrevs  = roi_abbrevs(idx_type);
  [tmp_roinames,idx_type_rois,idx_rois] = intersect(tmp_roinames,data.roinames,'stable');
  if isempty(idx_rois), continue; end
  nroi = length(idx_rois);

  % build column names
  colnames = cell(1,nroi);
  for i = 1:nroi
    abbrev = tmp_abbrevs{idx_type_rois(i)};
    colnames{i} = sprintf('%s__%s_mean',colstem_t,abbrev);
  end

  % build data rows
  header = cat(2,{'participant_id','session_id'},colnames);
  rows = cell(nsubs,2+nroi);
  for i = 1:nsubs
    d = i;
    vals = data.varmat(d,idx_rois);
    pv = num2cell(vals);
    pv(isnan(vals)) = {'n/a'};
    rows(i,:) = cat(2,{participant_ids{i},session_ids{i}},pv);
  end

  fprintf('%s: writing %s...\n',mfilename,fname_out);
  mmil_write_csv(fname_out,cat(1,header,rows),'separator','\t','quote_flag',0);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = build_bids_fname_var(parms,roitype)
  if ~isempty(parms.participant_id) && ~isempty(parms.task)
    base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
    base = sprintf('%s_atlas-%s',base,roitype);
    if ~isempty(parms.desc)
      base = sprintf('%s_desc-%s',base,parms.desc);
    end
    fname = sprintf('%s/%s_var.tsv',parms.outdir,base);
  else
    outstem = sprintf('%s%s__%s',parms.outstem,parms.tag,roitype);
    fname = sprintf('%s/%s.tsv',parms.outdir,outstem);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

