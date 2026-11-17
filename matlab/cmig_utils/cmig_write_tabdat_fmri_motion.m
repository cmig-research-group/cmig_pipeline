function cmig_write_tabdat_fmri_motion(varargin)
%function cmig_write_tabdat_fmri_motion([options])
%
% Purpose: write tabulated motion quality control data
%          from concatenated fMRI data (rsfmri or tfmri)
%
% Optional Input:
%   'indir': input directory containing roi_corrmat .mat file
%     {default = pwd}
%   'instem': input file stem for roi_corrmat .mat file
%     {default = 'roi_corrmat'}
%   'infostem': input file stem for info file (used in multi-subject mode)
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mr_y_fmri__qc__mot'}
%   'colstem': prefix attached to column names
%     {default = 'mr_y_fmri__qc__mot'}
%   'participant_id': participant ID; if empty, multi-subject mode
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'dofvec': session dof, in place of the one in the roi_corrmat file
%     single-subject mode only; NaN writes no dof columns
%     {default = []}
%   'dofvec_run': per-run dof, in place of the one in the roi_corrmat file
%     one value per run of the roi_corrmat file; NaN for a run with none
%     {default = []}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  05/20/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 09/24/2026 by Don Hagler
% Prev Mod: 10/03/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Last Mod: 11/15/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'instem','roi_corrmat',[],...
  'infostem','vol_info',[],...
  'outdir',pwd,[],...
  'outstem','mr_y_fmri__qc__mot',[],...
  'colstem','mr_y_fmri__qc__mot',[],...
  'task',[],[],...
  'participant_id',[],[],...
  'session_id',[],[],...
  'dofvec',[],[],...
  'dofvec_run',[],[],...
  'forceflag',false,[false true],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

if ~isempty(parms.participant_id) && ~isempty(parms.task)
  fname_out = sprintf('%s/%s_%s_task-%s_motionqc.tsv',...
    parms.outdir,parms.participant_id,parms.session_id,parms.task);
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
  % multi-subject mode: load participant IDs from info file
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

% load motion data from roi_corrmat.mat
fprintf('%s: loading %s...\n',mfilename,fname_in);
data = load(fname_in,'nvolvec','nvolfdsubvec','nvolfdsubcontvec','ntpointvec',...
  'meanfdvec','maxfdvec','meanrotvec','maxrotvec','meantransvec','maxtransvec',...
  'meanfdcvec','maxfdcvec','meanrotcvec','maxrotcvec','meantranscvec','maxtranscvec',...
  'meanfdmat','meanfdcmat',...
  'nvolfdsubdeconvvec','nvolfdsubdeconvvec_run',...
  'nvolvec_run','nvolfdsubvec_run','nvolfdsubcontvec_run','ntpointvec_run',...
  'meanfdvec_run','maxfdvec_run','meanrotvec_run','maxrotvec_run',...
  'meantransvec_run','maxtransvec_run','dofvec','dofvec_run',...
  'nframesvec','nframesvec_run','TR','nruns');

nsubs = length(participant_ids);
nruns = data.nruns;

% the deconvolution's frame count is absent from roi_corrmat files written before
%   it was added, and is then written as n/a
if ~isfield(data,'nvolfdsubdeconvvec')
  data.nvolfdsubdeconvvec = nan(size(data.nvolvec),'single');
  data.nvolfdsubdeconvvec_run = nan(size(data.nvolvec_run),'single');
end

% the caller's dof, from the 3dDeconvolve fits, replaces the roi_corrmat one
if ~isempty(parms.dofvec)
  if nsubs ~= 1
    fprintf('%s: ERROR: dofvec is for single-subject mode only\n',mfilename);
    return;
  end
  data.dofvec = single(parms.dofvec(1));
  dofvec_run = nan(1,nruns,'single');
  n = min(nruns,length(parms.dofvec_run));
  dofvec_run(1:n) = parms.dofvec_run(1:n);
  data.dofvec_run = dofvec_run;
end

% DOF columns are written only when dofvec is not all-NaN
dof_flag = any(isfinite(data.dofvec(:)));

% censored mean/max mot/rotat/trans columns are rsfmri-only (release 7.0
% data dictionary defines no equivalent tfmri fields)
rsfmri_flag = strcmp(parms.task,'rest');

mmil_mkdir(parms.outdir);

% build column names
colnames = build_colnames(parms.colstem, nruns, dof_flag, rsfmri_flag);

% build data matrix (one row per subject)
data_mat = cell(nsubs+1, 2+length(colnames));
data_mat(1,:) = cat(2,{'participant_id','session_id'},colnames);
for i = 1:nsubs
  d = i;
  row = build_row(data, d, nruns, dof_flag, rsfmri_flag);
  data_mat(i+1,:) = cat(2,{participant_ids{i},session_ids{i}},row);
end

fprintf('%s: writing %s...\n',mfilename,fname_out);
mmil_write_csv(fname_out,data_mat,'separator','\t','quote_flag',0);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function colnames = build_colnames(colstem, nruns, dof_flag, rsfmri_flag)
  colnames = {};
  colnames{end+1} = sprintf('%s__frame_count',           colstem);
  colnames{end+1} = sprintf('%s__tr_mean',               colstem);
  colnames{end+1} = sprintf('%s__frame__nodummy_count',  colstem);
  colnames{end+1} = sprintf('%s__frame__fdsub__02_count', colstem);
  colnames{end+1} = sprintf('%s__frame__fdsub__02__contig_count',colstem);
  colnames{end+1} = sprintf('%s__frame__censor_count',    colstem);
  if ~rsfmri_flag
    colnames{end+1} = sprintf('%s__frame__fdsub__09_count',colstem);
  end
  colnames{end+1} = sprintf('%s__mot_mean',               colstem);
  colnames{end+1} = sprintf('%s__mot_max',                colstem);
  if rsfmri_flag
    colnames{end+1} = sprintf('%s__mot__censor_mean',     colstem);
    colnames{end+1} = sprintf('%s__mot__censor_max',      colstem);
    colnames{end+1} = sprintf('%s__mot__xdir_mean',       colstem);
    colnames{end+1} = sprintf('%s__mot__ydir_mean',       colstem);
    colnames{end+1} = sprintf('%s__mot__zdir_mean',       colstem);
    colnames{end+1} = sprintf('%s__mot__xdir__censor_mean',colstem);
    colnames{end+1} = sprintf('%s__mot__ydir__censor_mean',colstem);
    colnames{end+1} = sprintf('%s__mot__zdir__censor_mean',colstem);
  end
  colnames{end+1} = sprintf('%s__rotat_mean',             colstem);
  colnames{end+1} = sprintf('%s__rotat_max',              colstem);
  if rsfmri_flag
    colnames{end+1} = sprintf('%s__rotat__censor_mean',   colstem);
    colnames{end+1} = sprintf('%s__rotat__censor_max',    colstem);
    colnames{end+1} = sprintf('%s__rotat__xdir_mean',     colstem);
    colnames{end+1} = sprintf('%s__rotat__ydir_mean',     colstem);
    colnames{end+1} = sprintf('%s__rotat__zdir_mean',     colstem);
    colnames{end+1} = sprintf('%s__rotat__xdir__censor_mean',colstem);
    colnames{end+1} = sprintf('%s__rotat__ydir__censor_mean',colstem);
    colnames{end+1} = sprintf('%s__rotat__zdir__censor_mean',colstem);
  end
  colnames{end+1} = sprintf('%s__trans_mean',             colstem);
  colnames{end+1} = sprintf('%s__trans_max',              colstem);
  if rsfmri_flag
    colnames{end+1} = sprintf('%s__trans__censor_mean',   colstem);
    colnames{end+1} = sprintf('%s__trans__censor_max',    colstem);
  end
  if dof_flag
    colnames{end+1} = sprintf('%s_dof',                   colstem);
  end
  % per-run columns go with the per-run betas, so tasks only; the correlation
  %   analysis's censoring counts are session-level only
  if rsfmri_flag, return; end
  for r = 1:nruns
    rtag = sprintf('r%02d',r);
    colnames{end+1} = sprintf('%s__%s__frame_count',            colstem,rtag);
    colnames{end+1} = sprintf('%s__%s_tr',                      colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__frame__nodummy_count',   colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__frame__fdsub__09_count',colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__mot_mean',               colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__mot_max',                colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__rotat_mean',             colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__rotat_max',              colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__trans_mean',             colstem,rtag);
    colnames{end+1} = sprintf('%s__%s__trans_max',              colstem,rtag);
    if dof_flag
      colnames{end+1} = sprintf('%s__%s_dof',                   colstem,rtag);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = build_row(data, d, nruns, dof_flag, rsfmri_flag)
  vals = [];
  vals(end+1) = data.nframesvec(d);
  vals(end+1) = data.TR;
  vals(end+1) = data.nvolvec(d);
  vals(end+1) = data.nvolfdsubvec(d);
  vals(end+1) = data.nvolfdsubcontvec(d);
  vals(end+1) = data.ntpointvec(d);
  if ~rsfmri_flag
    vals(end+1) = data.nvolfdsubdeconvvec(d);
  end
  vals(end+1) = data.meanfdvec(d);
  vals(end+1) = data.maxfdvec(d);
  if rsfmri_flag
    vals(end+1) = data.meanfdcvec(d);
    vals(end+1) = data.maxfdcvec(d);
    vals(end+1) = data.meanfdmat(d,1);
    vals(end+1) = data.meanfdmat(d,2);
    vals(end+1) = data.meanfdmat(d,3);
    vals(end+1) = data.meanfdcmat(d,1);
    vals(end+1) = data.meanfdcmat(d,2);
    vals(end+1) = data.meanfdcmat(d,3);
  end
  vals(end+1) = data.meanrotvec(d);
  vals(end+1) = data.maxrotvec(d);
  if rsfmri_flag
    vals(end+1) = data.meanrotcvec(d);
    vals(end+1) = data.maxrotcvec(d);
    vals(end+1) = data.meanfdmat(d,4);
    vals(end+1) = data.meanfdmat(d,5);
    vals(end+1) = data.meanfdmat(d,6);
    vals(end+1) = data.meanfdcmat(d,4);
    vals(end+1) = data.meanfdcmat(d,5);
    vals(end+1) = data.meanfdcmat(d,6);
  end
  vals(end+1) = data.meantransvec(d);
  vals(end+1) = data.maxtransvec(d);
  if rsfmri_flag
    vals(end+1) = data.meantranscvec(d);
    vals(end+1) = data.maxtranscvec(d);
  end
  if dof_flag
    vals(end+1) = data.dofvec(d);
  end
  % per-run values for tasks only, matching build_colnames
  if ~rsfmri_flag
    for r = 1:nruns
      vals(end+1) = data.nframesvec_run(d,r);
      vals(end+1) = data.TR;
      vals(end+1) = data.nvolvec_run(d,r);
      vals(end+1) = data.nvolfdsubdeconvvec_run(d,r);
      vals(end+1) = data.meanfdvec_run(d,r);
      vals(end+1) = data.maxfdvec_run(d,r);
      vals(end+1) = data.meanrotvec_run(d,r);
      vals(end+1) = data.maxrotvec_run(d,r);
      vals(end+1) = data.meantransvec_run(d,r);
      vals(end+1) = data.maxtransvec_run(d,r);
      if dof_flag
        vals(end+1) = data.dofvec_run(d,r);
      end
    end
  end
  row = num2cell(vals);
  row(isnan(vals)) = {'n/a'};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

