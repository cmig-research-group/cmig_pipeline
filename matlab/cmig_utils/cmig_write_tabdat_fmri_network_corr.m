function cmig_write_tabdat_fmri_network_corr(varargin)
%function cmig_write_tabdat_fmri_network_corr([options])
%
% Purpose: write tabulated network-to-network correlations from fMRI data
%          (rsfmri or tfmri)
%
% Optional Input:
%   'indir': input directory containing gpnet_corrmat .mat file
%     {default = pwd}
%   'instem': input file stem for gpnet_corrmat .mat file
%     {default = 'gpnet_corrmat'}
%   'infostem': input file stem for info file (used in multi-subject mode)
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mr_y_fmri__corr__gpnet'}
%   'colstem': column label stem
%     {default = 'mr_y_fmri__corr__gpnet'}
%   'bidirflag': [0|1] write columns for bidirectional correlations
%     {default = 1}
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
  'instem','gpnet_corrmat',[],...
  'infostem','vol_info',[],...
  'outdir',pwd,[],...
  'outstem','mr_y_fmri__corr__gpnet',[],...
  'colstem','mr_y_fmri__corr__gpnet',[],...
  'task',[],[],...
  'desc',[],[],...
  'bidirflag',true,[false true],...
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
  base = sprintf('%s_atlas-gpnet',base);
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
data = load(fname_in,'corrmat','ndirs','nnet','network_names','net1vec','net2vec');

% load network abbreviations
[net_names, net_abbrevs] = cmig_load_abbrev('indir',parms.tabdir,'instem','gpnet_abbrev');
if isempty(net_names), return; end

nsubs = length(participant_ids);
nnet = data.nnet;
ndirs = data.ndirs;
mmil_mkdir(parms.outdir);

if parms.bidirflag
  % restore full correlation matrix [ndirs x nnet*nnet]
  idx_triu = find(triu(ones(nnet,nnet)));
  idx_tril = find(tril(ones(nnet,nnet)));
  corrmat_full = nan(nnet*nnet,ndirs,'single');
  corrmat_full(idx_triu,:) = data.corrmat';
  corrmat_lower = reshape(permute(reshape(corrmat_full,[nnet,nnet,ndirs]),[2,1,3]),[nnet*nnet,ndirs]);
  corrmat_full(idx_tril,:) = corrmat_lower(idx_tril,:);
  clear corrmat_lower;
  corrmat_full = corrmat_full';
  % column names: all combinations in both orders
  colnames = cell(nnet,nnet);
  for i = 1:nnet
    for j = 1:nnet
      abbrev1 = net_abbrevs{find(strcmp(data.network_names{i},net_names))};
      abbrev2 = net_abbrevs{find(strcmp(data.network_names{j},net_names))};
      colnames{i,j} = sprintf('%s__%s__%s_mean',parms.colstem,abbrev1,abbrev2);
    end
  end
  colnames = reshape(colnames',[1,nnet*nnet]);
  corrmat_out = corrmat_full;
else
  % upper triangle only
  nnodes = size(data.corrmat,2);
  colnames = cell(1,nnodes);
  for i = 1:nnodes
    abbrev1 = net_abbrevs{find(strcmp(data.network_names{data.net1vec(i)},net_names))};
    abbrev2 = net_abbrevs{find(strcmp(data.network_names{data.net2vec(i)},net_names))};
    colnames{i} = sprintf('%s__%s__%s_mean',parms.colstem,abbrev1,abbrev2);
  end
  corrmat_out = data.corrmat;
end

% build full data matrix
header = cat(2,{'participant_id','session_id'},colnames);
rows = cell(nsubs,2+size(corrmat_out,2));
for i = 1:nsubs
  d = i;
  vals = corrmat_out(d,:);
  pv = num2cell(vals);
  pv(isnan(vals)) = {'n/a'};
  rows(i,:) = cat(2,{participant_ids{i},session_ids{i}},pv);
end

fprintf('%s: writing %s...\n',mfilename,fname_out);
mmil_write_csv(fname_out,cat(1,header,rows),'separator','\t','quote_flag',0);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

