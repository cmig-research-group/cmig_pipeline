function cmig_write_tabdat_fmri_network_aseg_corr(varargin)
%function cmig_write_tabdat_fmri_network_aseg_corr([options])
%
% Purpose: write tabulated network-to-subcortical ROI correlations
%          from fMRI data (rsfmri or tfmri)
%
% Optional Input:
%   'indir': input directory containing gpnet_aseg_corrmat .mat file
%     {default = pwd}
%   'instem': input file stem for gpnet_aseg_corrmat .mat file
%     {default = 'gpnet_aseg_corrmat'}
%   'infostem': input file stem for info file (used in multi-subject mode)
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mr_y_fmri__corr__gpnet__aseg'}
%   'colstem': column label stem
%     {default = 'mr_y_fmri__corr__gpnet__aseg'}
%   'tabdir': directory containing abbreviation lookup table CSV files
%     {default = [getenv('MMPS_PARMS') '/tabulated']}
%   'participant_id': participant ID; if empty, multi-subject mode
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  05/20/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'instem','gpnet_aseg_corrmat',[],...
  'infostem','vol_info',[],...
  'outdir',pwd,[],...
  'outstem','mr_y_fmri__corr__gpnet__aseg',[],...
  'colstem','mr_y_fmri__corr__gpnet__aseg',[],...
  'task',[],[],...
  'desc',[],[],...
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

if ~isempty(parms.participant_id) && ~isempty(parms.task)
  base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
  base = sprintf('%s_atlas-gpnetaseg',base);
  if ~isempty(parms.desc)
    base = sprintf('%s_desc-%s',base,parms.desc);
  end
  fname_out = sprintf('%s/%s_corr.tsv',parms.outdir,base);
else
  fname_out = sprintf('%s/%s.tsv',parms.outdir,parms.outstem);
end
if exist(fname_out,'file') && ~parms.forceflag, return; end

fname_in = sprintf('%s/%s.mat',parms.indir,parms.instem);
if ~exist(fname_in,'file')
  fprintf('%s: WARNING: input file not found: %s\n',mfilename,fname_in);
  return;
end

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
data = load(fname_in,'corrmat','network_names','roinames','net1vec','roi2vec');

% load network and aseg abbreviations
[net_names,  net_abbrevs]  = cmig_load_abbrev('indir',parms.tabdir,'instem','gpnet_abbrev');
[aseg_names, aseg_abbrevs] = cmig_load_abbrev('indir',parms.tabdir,'instem','aseg_abbrev');
if isempty(net_names) || isempty(aseg_names), return; end

% build column names
nnodes = size(data.corrmat,2);
colnames = cell(1,nnodes);
for i = 1:nnodes
  net_name  = data.network_names{data.net1vec(i)};
  net_abbrev  = net_abbrevs{find(strcmp(net_name,net_names))};
  roi_name  = data.roinames{data.roi2vec(i)};
  roi_abbrev  = aseg_abbrevs{find(strcmp(roi_name,aseg_names))};
  colnames{i} = sprintf('%s__%s__%s_mean',parms.colstem,net_abbrev,roi_abbrev);
end

nsubs = length(participant_ids);
mmil_mkdir(parms.outdir);

header = cat(2,{'participant_id','session_id'},colnames);
rows = cell(nsubs,2+nnodes);
for i = 1:nsubs
  d = i;
  vals = data.corrmat(d,:);
  pv = num2cell(vals);
  pv(isnan(vals)) = {'n/a'};
  rows(i,:) = cat(2,{participant_ids{i},session_ids{i}},pv);
end

fprintf('%s: writing %s...\n',mfilename,fname_out);
mmil_write_csv(fname_out,cat(1,header,rows),'separator','\t','quote_flag',0);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

