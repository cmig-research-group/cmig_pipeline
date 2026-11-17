function cmig_write_fmri_corrmat(varargin)
%function cmig_write_fmri_corrmat([options])
%
% Purpose: write ROI-ROI correlation matrix as a BIDS derivative TSV file
%   output format: (nroi+1) x (nroi+1) matrix with ROI name labels on rows and columns
%   includes all ROIs present in the input corrmat .mat file
%
% Optional Input ('key',value):
%   'indir': input directory containing corrmat .mat file
%     {default = pwd}
%   'instem': input file stem for corrmat .mat file
%     {default = 'roi_corrmat'}
%   'outdir': output directory
%     {default = pwd}
%   'atlas': atlas label for BIDS filename; if empty, omitted from filename
%     {default = []}
%   'task': task name for BIDS filename
%     {default = []}
%   'desc': desc entity for BIDS filename; if empty, omitted from filename
%     {default = []}
%   'participant_id': participant ID; if empty, uses instem as output base name
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  05/29/2026 by Don Hagler
% Prev Mod: 05/29/2026 by Don Hagler
% Prev Mod: 06/05/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'instem','roi_corrmat',[],...
  'outdir',pwd,[],...
  'atlas',[],[],...
  'task',[],[],...
  'desc',[],[],...
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

fname_out = build_fname(parms);
if exist(fname_out,'file') && ~parms.forceflag, return; end

fname_in = sprintf('%s/%s.mat',parms.indir,parms.instem);
if ~exist(fname_in,'file')
  fprintf('%s: WARNING: input file not found: %s\n',mfilename,fname_in);
  return;
end

data = load(fname_in,'corrmat','roinames','roi1vec','roi2vec');
nroi = length(data.roinames);
if size(data.corrmat,1) < 1 || all(isnan(data.corrmat(1,:)))
  fprintf('%s: WARNING: no valid data in %s\n',mfilename,fname_in);
  return;
end

% reconstruct full symmetric nroi x nroi matrix from upper triangle
R = NaN(nroi,nroi,'single');
R(sub2ind([nroi,nroi],data.roi1vec(:),data.roi2vec(:))) = data.corrmat(1,:);
idx_lower = find(tril(ones(nroi,nroi),-1));
[row_l,col_l] = ind2sub([nroi,nroi],idx_lower);
R(idx_lower) = R(sub2ind([nroi,nroi],col_l,row_l));

% write matrix TSV
mmil_mkdir(parms.outdir);
fid = fopen(fname_out,'wt');
if fid < 0
  fprintf('%s: ERROR: failed to open %s for writing\n',mfilename,fname_out);
  return;
end

% header row: 'roi_name' then all ROI names
header = ['roi_name', data.roinames(:)'];
fprintf(fid,'%s\n',strjoin(header,'\t'));

% data rows: ROI name then correlation values
for i = 1:nroi
  vals = R(i,:);
  parts = cell(1,nroi+1);
  parts{1} = data.roinames{i};
  for j = 1:nroi
    if isnan(vals(j))
      parts{j+1} = 'n/a';
    else
      parts{j+1} = num2str(vals(j),'%.6g');
    end
  end
  fprintf(fid,'%s\n',strjoin(parts,'\t'));
end
fclose(fid);
fprintf('%s: wrote %s\n',mfilename,fname_out);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = build_fname(parms)
  if ~isempty(parms.participant_id) && ~isempty(parms.task)
    base = sprintf('%s_%s_task-%s',parms.participant_id,parms.session_id,parms.task);
    if ~isempty(parms.atlas)
      base = sprintf('%s_atlas-%s',base,parms.atlas);
    end
    if ~isempty(parms.desc)
      base = sprintf('%s_desc-%s',base,parms.desc);
    end
  else
    base = parms.instem;
  end
  fname = sprintf('%s/%s_corrmat.tsv',parms.outdir,base);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
