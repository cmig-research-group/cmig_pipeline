function [outdir,errcode,errmsg] = cmig_concat_bids_sessions(datatype,analysis_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_concat_bids_sessions(datatype,analysis_rootdir,output_rootdir,[options])
%
% Purpose: concatenate per-session concat-format mat files across
%   subject-sessions, producing whole-cohort files in the release layout
%
% Required Input:
%   datatype: which set of outputs to concatenate
%     'smri_measures'  vertex-wise sMRI measures      (measmat)
%     'smri_surfaces'  ico surface coordinates        (coordmat)
%     'dmri_dti'       vertex-wise DTI measures       (measmat)
%     'dmri_isdti'     isotropic-suppressed DTI       (measmat)
%     'dmri_rsi'       vertex-wise RSI measures       (measmat)
%     'tfmri_stats'    task fMRI beta and sem maps    (measmat)
%     'fmri_corr'      ROI correlation matrices       (corrmat, varmat)
%     'fmri_tseries'   ROI and ico time series        (datamat_tsdata,
%                                                      datamat_motion, censvec)
%   analysis_rootdir: root of the per-session cmig_analysis tree, or a cell
%     array of such roots, whose sessions are combined; a session found under
%     more than one root is an error
%   output_rootdir: root directory for the concatenated output
%
% Optional Parameters ('key',value pairs):
%   'sessions': struct array with fields participant_id, session_id and indir,
%     bypassing the tree scan
%     {default = []}
%   'seslist': cell array of session ids to keep, e.g. {'ses-00A'}
%     if empty, keep all
%     {default = []}
%   'id_events_file': id_events release csv; only the visits it lists
%     (its participant_id and session_id columns) are kept
%     if empty, keep all
%     {default = []}
%   'measlist': cell array of measure names; if empty, the datatype's default:
%     smri_measures: {'thk','area','sulc','t1_gwc','t2_gwc','t1_gm','t2_gm',
%                     't1_wm','t2_wm'}
%     dmri_dti, dmri_isdti: {'FA','MD','LD','TD'}
%     dmri_rsi: {'RNI','RND','RNT','HNT','FNI'}
%     tfmri_stats: {'beta','sem'}
%     not used by smri_surfaces, fmri_corr and fmri_tseries
%     {default = []}
%   'layerlist': cell array of surface layer labels, dmri datatypes only
%     {default = {'gm','wm'}}
%   'surflist': cell array of surface names, smri_surfaces only
%     {default = {'white','pial'}}
%   'smoothlist': sphere smoothing steps
%     {default = [0,16,256,1000]}
%   'hemilist': cell array of hemisphere labels
%     {default = {'lh','rh'}}
%   'tasklist': cell array of task labels, fmri datatypes only
%     {default = {'mid','nback','sst','rest'}}
%   'atlaslist': cell array of atlas labels, the fmri datatypes only.
%     fmri_tseries is the expensive one: ico5 rest is 81 MB per session
%     against 2.8 MB for the gpdstdskaseg union atlas
%     {default = {'gpaseg'}}
%   'ico_order': icosahedral order of the inputs
%     {default = 5}
%   'nruns_{task}': number of runs the output allows for each fmri task,
%     for rest, mid, nback and sst. Runs beyond this are not read. The trimmed
%     fmri_corr outputs (_trim5, _trim10) have nruns 1 regardless, since each
%     session's trimmed corrmat comes from one concatenated series
%     {default = 4, 2, 2, 2}
%   'ntpoints_{task}': number of time points the output allows for each fmri
%     task. A session with more is truncated
%     {default = 375, 403, 362, 437}
%   'vol_info_flag': [0|1] write vol_info.mat and return without
%     concatenating anything
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output files {default = 0}
%
% Output:
%   outdir: output directory, {output_rootdir}/{folder} (below); for
%     fmri_corr and fmri_tseries, {output_rootdir}/corrmat
%   errcode: 0 if successful, 1 if error; a datatype with no data (no
%     session with one of its per-session files, or no input file for any
%     target) is skipped with a warning, errcode 0, and nothing kept
%   errmsg: error message string if errcode=1, otherwise []
%
% Output folders follow the release 7.0 layout, under output_rootdir:
%   smri_measures  vertexwise/smri     smri_surfaces  vertexwise/surf
%   tfmri_stats    vertexwise/tfmri    dmri_dti       vertexwise/dti
%   dmri_isdti     vertexwise/isdti    dmri_rsi       vertexwise/rsi
%   fmri_corr and fmri_tseries: corrmat/rsfmri (rest) and corrmat/tfmri (the
%     other tasks)
% Output naming follows the release layout, except that the smoothing axis is
%   named by fwhm rather than by smoothing steps:
%   {folder}/{stem}_fwhm{F}_{hemi}.mat  and  {folder}/vol_info.mat
%
% Status files (tsv, an empty value n/a):
%   {folder}/{stem}_status.tsv, beside each .mat when it is written (or the
%     target fails): procdate, errcode, errmsg, startdate, nsess, nfound,
%     nmissing, npad_run, ntrunc_run, npad_time, ntrunc_time
%   {output_rootdir}/concat_{datatype}_status.tsv, or
%     concat_{datatype}_{meas}_sm{smooth}_status.tsv for smri_measures and the
%     dmri datatypes given one measure and one smoothing step, once per
%     invocation: procdate, errcode, errmsg, startdate, vol_info_flag, nsess
%     (the most of any of its folders), ntargets, nwritten,
%     nskipped_existing, nnodata. It is named as cmig_mkjobs_bids_concat.py
%     names the job, whose smoothlist value, in smoothing steps, is the
%     {smooth} of sm{smooth}, where the output files have the fwhm{F} it
%     gives. A no-data skip has errcode 0 and errmsg
%     "no {datatype} data in {roots}"
%
% NOTE: session i of vol_info is index i along each variable's session axis,
%   which is dim 1 except for the tseries arrays, where dim 1 is the run and
%   the session axis is dim 2. vol_info is written before any measure file and
%   read back by later calls, so every file in a folder shares one session
%   axis no matter how many jobs write them. Run once with vol_info_flag=1
%   before the jobs that do the concatenating
%
% NOTE: fmri_corr and fmri_tseries share the vol_info.mat of each corrmat
%   folder: the first to run writes it and the other reads it back, so both
%   are to be given the same sessions. With forceflag, a vol_info.mat with
%   the same rows is kept, not rewritten; one with other rows is rewritten,
%   and the other datatype is then to be rerun with forceflag too. A
%   vol_info.mat is never removed while other .mat files in its folder use it
%
% NOTE: tfmri_stats per-run outputs are by run slot, {stem}_r01 the first
%   selected run and so on to nruns_{task}: slot k is filled from the run with
%   slot_num k in the session's _desc-select_runs.tsv, or, without one, from
%   the k-th of the session's per-run files in run order
%
% Created:  08/28/2026 by Don Hagler
% Prev Mod: 08/28/2026 by Don Hagler
% Prev Mod: 08/31/2026 by Don Hagler
% Prev Mod: 09/02/2026 by Don Hagler
% Prev Mod: 09/14/2026 by Don Hagler
% Prev Mod: 10/12/2026 by Don Hagler
% Prev Mod: 11/04/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/13/2026 by Don Hagler
% Prev Mod: 11/14/2026 by Don Hagler
% Last Mod: 11/15/2026 by Don Hagler
%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,3), errcode = 1; return; end
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

% coerce char values from a compiled deployment before the typed ranges below
%   are checked against them
varargin = cmig_parse_deployed_args(varargin);
% several roots arrive from a compiled deployment as one string, {'a' 'b'}
if isdeployed && ischar(analysis_rootdir) && ~isempty(regexp(analysis_rootdir,'^{','once'))
  analysis_rootdir = cmig_parseval(analysis_rootdir);
end

parms = mmil_args2parms(varargin,{...
  'datatype',datatype,[],...
  'analysis_rootdir',analysis_rootdir,[],...
  'output_rootdir',output_rootdir,[],...
  'sessions',[],[],...
  'seslist',[],[],...
  'id_events_file',[],[],...
  'measlist',[],[],...
  'layerlist',{'gm','wm'},[],...
  'surflist',{'white','pial'},[],...
  'smoothlist',[0,16,256,1000],[],...
  'hemilist',{'lh','rh'},[],...
  'tasklist',{'mid','nback','sst','rest'},[],...
  'atlaslist',{'gpaseg'},[],...
  'ico_order',5,[0 7],...
  'nruns_rest',4,[],...
  'nruns_mid',2,[],...
  'nruns_nback',2,[],...
  'nruns_sst',2,[],...
  'ntpoints_rest',375,[],...
  'ntpoints_mid',403,[],...
  'ntpoints_nback',362,[],...
  'ntpoints_sst',437,[],...
  'vol_info_flag',false,[false true],...
  'forceflag',false,[false true],...
});

if ischar(parms.analysis_rootdir)
  parms.analysis_rootdir = {parms.analysis_rootdir};
end
if ischar(parms.tasklist), parms.tasklist = {parms.tasklist}; end

% one status file per invocation, in output_rootdir, since fmri_corr and
%   fmri_tseries span two folders; a status left by an earlier invocation is
%   removed first, so a run that stops without writing one is not read as done
fname_status = sprintf('%s/%s_status.tsv',parms.output_rootdir,status_stem(parms));
if exist(fname_status,'file'), delete(fname_status); end

[outdir,errcode,errmsg,counts] = concat_datatype(parms);

write_concat_status(fname_status,errcode,errmsg,startdate,counts);
% a no-data skip carries its message in the status file but is not an error
if ~errcode, errmsg = []; end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outdir,errcode,errmsg,counts] = concat_datatype(parms)
  % the work of one invocation, folder by folder. counts are the invocation
  %   status columns, empty (n/a) where the step that sets them was not
  %   reached. A datatype with no data in any folder returns errcode 0 with
  %   errmsg "no <datatype> data in ..."
  outdir = [];
  counts = struct('vol_info_flag',double(parms.vol_info_flag),'nsess',[],...
    'ntargets',[],'nwritten',[],'nskipped_existing',[],'nnodata',[]);
  [parms,errcode,errmsg] = check_datatype(parms);
  if errcode, return; end
  % folders are made only when vol_info is written, so a datatype with no
  %   data leaves no empty folder
  outdir = sprintf('%s/%s',parms.output_rootdir,parms.outsubdir);
  ndata = 0;
  for f = 1:length(parms.folders)
    fparms = parms;
    fparms.tasklist = parms.folders(f).tasklist;
    fparms.data_patterns = parms.folders(f).data_patterns;
    fdir = sprintf('%s/%s',parms.output_rootdir,parms.folders(f).subdir);
    [errcode,errmsg,fcounts,nodata,scan] = concat_folder(fparms,fdir);
    counts = add_counts(counts,fcounts);
    if errcode, return; end
    if ~nodata, ndata = ndata + 1; end
    % the next folder takes its rows from the same scan
    if isempty(parms.sessions), parms.sessions = scan; end
  end
  % a datatype with no data is skipped, not an error: nothing is written
  if ~ndata
    errmsg = sprintf('no %s data in %s',parms.datatype,...
      strjoin(parms.analysis_rootdir,', '));
    fprintf('%s: WARNING: %s\n',mfilename,errmsg);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,counts,nodata,scan] = concat_folder(parms,outdir)
  % the targets of one output folder. nodata is true where the folder has no
  %   data, and nothing of this call is kept in it. scan is the session list
  %   of a tree scan, if write_vol_info made one
  counts = struct('nsess',[],'ntargets',[],'nwritten',[],...
    'nskipped_existing',[],'nnodata',[]);
  nodata = true;

  % fix the rows before anything is concatenated. Jobs of one datatype, and
  %   fmri_corr and fmri_tseries, share a folder, so the session list must
  %   come from a single scan: the first call writes vol_info and the rest
  %   read it, keeping every file of the folder on the same session axis
  [parms,errcode,errmsg] = write_vol_info(parms,outdir);
  scan = parms.scan;
  if errcode, return; end
  nsess = length(parms.sessions);
  counts.nsess = nsess;
  if ~nsess
    fprintf('%s: WARNING: no %s data for %s\n',mfilename,parms.datatype,outdir);
    return;
  end
  fprintf('%s: %d sessions for datatype %s in %s\n',...
    mfilename,nsess,parms.datatype,outdir);
  if parms.vol_info_flag, nodata = false; return; end

  % concatenate each target
  if strcmp(parms.datatype,'tfmri_stats'), parms = get_run_slots(parms); end
  targets = get_targets(parms);
  ntargets = length(targets);
  fprintf('%s: %d targets\n',mfilename,ntargets);
  counts.ntargets = ntargets;
  counts.nwritten = 0; counts.nskipped_existing = 0; counts.nnodata = 0;
  for t = 1:ntargets
    [errcode,errmsg,result] = concat_target(targets(t),parms,outdir);
    if errcode, return; end
    switch result
      case 'written',  counts.nwritten = counts.nwritten + 1;
      case 'existing', counts.nskipped_existing = counts.nskipped_existing + 1;
      case 'nodata',   counts.nnodata = counts.nnodata + 1;
    end
  end
  % no input file for any target is no data too: a vol_info this run wrote
  %   is removed, unless other files of the folder use it
  if ~counts.nwritten && ~counts.nskipped_existing
    fprintf('%s: WARNING: no %s data for %s\n',mfilename,parms.datatype,outdir);
    if parms.vol_info_written, remove_vol_info(outdir); end
    return;
  end
  nodata = false;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function counts = add_counts(counts,fcounts)
  % a folder's counts into the invocation's: nsess is the most of any
  %   folder, the others are summed; an empty count is not reached
  names = fieldnames(fcounts);
  for i = 1:length(names)
    n = names{i};
    b = fcounts.(n);
    if isempty(b), continue; end
    a = counts.(n);
    if isempty(a)
      counts.(n) = b;
    elseif strcmp(n,'nsess')
      counts.(n) = max(a,b);
    else
      counts.(n) = a + b;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = status_stem(parms)
  % the invocation status is named as cmig_mkjobs_bids_concat.py names the
  %   job: concat_{datatype}, or concat_{datatype}_{meas}_sm{smooth} for a
  %   datatype expanded by measure and smoothing, given one of each
  expandable = {'smri_measures','dmri_dti','dmri_isdti','dmri_rsi'};
  measlist = parms.measlist;
  if ischar(measlist), measlist = {measlist}; end
  if ismember(parms.datatype,expandable) && iscell(measlist) &&...
     numel(measlist)==1 && numel(parms.smoothlist)==1
    fstem = sprintf('concat_%s_%s_sm%g',parms.datatype,measlist{1},...
      parms.smoothlist);
  else
    fstem = sprintf('concat_%s',parms.datatype);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_concat_status(fname,errcode,errmsg,startdate,extra)
  % a dataset-level status file with the columns of cmig_write_status but no
  %   participant_id or session_id: procdate, errcode, errmsg, startdate, then
  %   the fields of extra; an empty value is written as n/a
  if ~ischar(errmsg) || isempty(errmsg), errmsg = 'n/a'; end
  errmsg = strrep(strrep(errmsg,sprintf('\n'),' '),sprintf('\r'),' ');
  procdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');
  names = {}; vals = {};
  if ~isempty(extra)
    names = fieldnames(extra);
    vals = cell(size(names));
    for i = 1:length(names)
      val = extra.(names{i});
      if isnumeric(val) || islogical(val)
        vals{i} = num2str(val);
      elseif ischar(val)
        vals{i} = val;
      else
        vals{i} = '';
      end
      if isempty(vals{i}), vals{i} = 'n/a'; end
    end
  end
  mmil_mkdir(fileparts(fname));
  fid = fopen(fname,'wt');
  if fid < 0
    error('%s: failed to open %s for writing',mfilename,fname);
  end
  fprintf(fid,'procdate\terrcode\terrmsg\tstartdate');
  for i = 1:length(names), fprintf(fid,'\t%s',names{i}); end
  fprintf(fid,'\n');
  fprintf(fid,'%s\t%d\t%s\t%s',procdate,errcode,errmsg,startdate);
  for i = 1:length(vals), fprintf(fid,'\t%s',vals{i}); end
  fprintf(fid,'\n');
  fclose(fid);
  fprintf('%s: wrote file %s\n',mfilename,fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_datatype(parms)
  errcode = 0; errmsg = [];
  % modality is the BIDS derivative subdirectory the inputs sit in.
  %   data_patterns are the per-session files after {pid}_{sid}, any one of
  %   which is data (the same patterns as cmig_run_bids_module.py
  %   CONCAT_DATA). outsubdir is the output folder under output_rootdir, as
  %   in the release 7.0 layout
  switch parms.datatype
    case 'smri_measures'
      parms.modality = 'anat';
      parms.outsubdir = 'vertexwise/smri';
      parms.data_patterns = {'*_hemi-*_atlas-ico*_desc-fwhm*.mat'};
      if isempty(parms.measlist)
        parms.measlist = {'thk','area','sulc','t1_gwc','t2_gwc',...
                          't1_gm','t2_gm','t1_wm','t2_wm'};
      end
    case 'smri_surfaces'
      parms.modality = 'anat';
      parms.outsubdir = 'vertexwise/surf';
      parms.data_patterns = {'*_hemi-*_atlas-ico*_white.mat',...
                             '*_hemi-*_atlas-ico*_pial.mat'};
    case {'dmri_dti','dmri_isdti','dmri_rsi'}
      parms.modality = 'dwi';
      model = regexprep(parms.datatype,'^dmri_','');
      parms.outsubdir = ['vertexwise/' model];
      parms.data_patterns = {sprintf('*_desc-fwhm*_%s*.mat',model)};
      if strcmp(parms.datatype,'dmri_rsi')
        if isempty(parms.measlist)
          parms.measlist = {'RNI','RND','RNT','HNT','FNI'};
        end
      elseif isempty(parms.measlist)
        parms.measlist = {'FA','MD','LD','TD'};
      end
    case 'tfmri_stats'
      parms.modality = 'func';
      parms.outsubdir = 'vertexwise/tfmri';
      parms.data_patterns = {'*_hemi-*_atlas-ico*_desc-*_beta.mat'};
      if isempty(parms.measlist), parms.measlist = {'beta','sem'}; end
    case 'fmri_corr'
      parms.modality = 'func';
      parms.outsubdir = 'corrmat';
      parms.data_patterns = {'*_task-*_atlas-*_corrmat.mat'};
    case 'fmri_tseries'
      parms.modality = 'func';
      parms.outsubdir = 'corrmat';
      parms.data_patterns = {'*_task-*_atlas-*_tseries.mat'};
    otherwise
      errcode = 1;
      errmsg = sprintf('unrecognized datatype %s',parms.datatype);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
  end
  parms.folders = get_folders(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function folders = get_folders(parms)
  % the output folders, each with its tasks and its data patterns: outsubdir
  %   alone, except for fmri_corr and fmri_tseries, whose rest outputs go to
  %   corrmat/rsfmri and those of the other tasks to corrmat/tfmri, as in the
  %   release 7.0 layout; a folder's data patterns name its tasks
  folders = struct('subdir',{},'tasklist',{},'data_patterns',{});
  if ~ismember(parms.datatype,{'fmri_corr','fmri_tseries'})
    folders(1).subdir = parms.outsubdir;
    folders(1).tasklist = parms.tasklist;
    folders(1).data_patterns = parms.data_patterns;
    return;
  end
  isrest = strcmp(parms.tasklist,'rest');
  sets = {{'rsfmri',parms.tasklist(isrest)},{'tfmri',parms.tasklist(~isrest)}};
  for s = 1:length(sets)
    tasks = sets{s}{2};
    if isempty(tasks), continue; end
    pats = {};
    for t = 1:length(tasks)
      pats = [pats,regexprep(parms.data_patterns,'task-\*',...
        sprintf('task-%s',tasks{t}))]; %#ok<AGROW>
    end
    folders(end+1).subdir = sprintf('%s/%s',parms.outsubdir,sets{s}{1}); %#ok<AGROW>
    folders(end).tasklist = tasks;
    folders(end).data_patterns = pats;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sessions,errcode,errmsg] = get_session_list(parms)
  % the session list comes from the derivatives trees rather than an
  %   inventory. Sessions of several roots are sorted by participant_id and
  %   session_id, the order dir gives within one root, so the rows do not
  %   depend on the order of the roots
  sessions = []; errcode = 0; errmsg = [];
  visits = [];
  if ~isempty(parms.id_events_file)
    [visits,errcode,errmsg] = read_id_events(parms.id_events_file);
    if errcode, return; end
  end
  pids = {}; sids = {}; indirs = {};
  for r = 1:length(parms.analysis_rootdir)
    rootdir = parms.analysis_rootdir{r};
    if ~exist(rootdir,'dir')
      errcode = 1;
      errmsg = sprintf('analysis root %s not found',rootdir);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    fprintf('%s: scanning %s...\n',mfilename,rootdir);
    slist = dir(sprintf('%s/sub-*',rootdir));
    slist = slist([slist.isdir]);
    for i = 1:length(slist)
      pid = slist(i).name;
      vlist = dir(sprintf('%s/%s/ses-*',rootdir,pid));
      vlist = vlist([vlist.isdir]);
      for j = 1:length(vlist)
        sid = vlist(j).name;
        if ~isempty(parms.seslist) && ~ismember(sid,parms.seslist), continue; end
        indir = sprintf('%s/%s/%s/%s',rootdir,pid,sid,parms.modality);
        if ~exist(indir,'dir'), continue; end
        pids{end+1} = pid; sids{end+1} = sid; indirs{end+1} = indir; %#ok<AGROW>
      end
    end
  end
  if isempty(pids), return; end
  ids = strcat(pids,'_',sids);
  if ~isempty(parms.id_events_file)
    keep = ismember(ids,visits);
    left_out = indirs(~keep);
    pids = pids(keep); sids = sids(keep); indirs = indirs(keep); ids = ids(keep);
    fprintf('%s: %d sessions not in %s, left out\n',...
      mfilename,length(left_out),parms.id_events_file);
    for i = 1:length(left_out)
      fprintf('%s:   left out %s\n',mfilename,left_out{i});
    end
  end
  if isempty(pids), return; end
  % the same session under two roots has no single input; which to keep is
  %   for the caller to resolve
  [ids_u,~,k] = unique(ids);
  ndup = accumarray(k(:),1);
  if any(ndup > 1)
    for d = find(ndup(:)' > 1)
      fprintf('%s: ERROR: %s found in: %s\n',mfilename,ids_u{d},...
        strjoin(indirs(k==d),', '));
    end
    errcode = 1;
    errmsg = sprintf('%d sessions found under more than one analysis root',...
      sum(ndup > 1));
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  % sort is stable, so sorting by session_id and then by participant_id
  %   orders by participant_id, then session_id
  [~,o1] = sort(sids);
  [~,o2] = sort(pids(o1));
  order = o1(o2);
  sessions = struct('participant_id',pids(order),'session_id',sids(order),...
                    'indir',indirs(order));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [visits,errcode,errmsg] = read_id_events(fname)
  % the participant_id and session_id columns of an id_events release csv,
  %   joined as {pid}_{sid}; it lists every kind of visit, imaging or not
  visits = {}; errcode = 0; errmsg = [];
  if ~exist(fname,'file')
    errcode = 1;
    errmsg = sprintf('id_events_file %s not found',fname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  fid = fopen(fname,'rt');
  hdr = strtrim(strsplit(fgetl(fid),',','CollapseDelimiters',false));
  i_pid = find(strcmp(hdr,'participant_id'));
  i_sid = find(strcmp(hdr,'session_id'));
  if isempty(i_pid) || isempty(i_sid)
    fclose(fid);
    errcode = 1;
    errmsg = sprintf('no participant_id and session_id columns in %s',fname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  C = textscan(fid,repmat('%s',1,length(hdr)),'Delimiter',',');
  fclose(fid);
  visits = unique(strcat(strtrim(C{i_pid}),'_',strtrim(C{i_sid})));
  fprintf('%s: %d visits in %s\n',mfilename,length(visits),fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = get_targets(parms)
  % one entry per output file: infix is the per-session file name with the
  %   participant and session prefix removed, outstem is the output file stem
  targets = [];
  switch parms.datatype
    case 'smri_measures'
      for m = 1:length(parms.measlist)
        for s = 1:length(parms.smoothlist)
          fwhm = round(cmig_smoothsteps_to_fwhm(parms.smoothlist(s)));
          for h = 1:length(parms.hemilist)
            hemi = parms.hemilist{h};
            bids_meas = strrep(parms.measlist{m},'_','');
            targets = add_target(targets,...
              sprintf('hemi-%s_atlas-ico%d_desc-fwhm%d_%s',...
                cmig_get_bids_hemi(hemi),parms.ico_order,fwhm,bids_meas),...
              sprintf('%s_fwhm%d_%s',parms.measlist{m},fwhm,hemi),'measmat');
          end
        end
      end
    case 'smri_surfaces'
      for s = 1:length(parms.surflist)
        for h = 1:length(parms.hemilist)
          hemi = parms.hemilist{h};
          targets = add_target(targets,...
            sprintf('hemi-%s_atlas-ico%d_%s',...
              cmig_get_bids_hemi(hemi),parms.ico_order,parms.surflist{s}),...
            sprintf('%s_%s',parms.surflist{s},hemi),'coordmat');
        end
      end
    case {'dmri_dti','dmri_isdti','dmri_rsi'}
      % the model prefix is part of the per-session file name, so it selects
      %   dtiFAgm against isdtiFAgm against rsiRNIgm
      switch parms.datatype
        case 'dmri_dti',   model = 'dti';
        case 'dmri_isdti', model = 'isdti';
        case 'dmri_rsi',   model = 'rsi';
      end
      if strcmp(parms.datatype,'dmri_isdti')
        outfix = 'isdti_';
      else
        outfix = '';
      end
      for m = 1:length(parms.measlist)
        for l = 1:length(parms.layerlist)
          layer = parms.layerlist{l};
          for s = 1:length(parms.smoothlist)
            fwhm = round(cmig_smoothsteps_to_fwhm(parms.smoothlist(s)));
            for h = 1:length(parms.hemilist)
              hemi = parms.hemilist{h};
              targets = add_target(targets,...
                sprintf('hemi-%s_atlas-ico%d_desc-fwhm%d_%s%s%s',...
                  cmig_get_bids_hemi(hemi),parms.ico_order,fwhm,...
                  model,parms.measlist{m},layer),...
                sprintf('%s%s_%s_fwhm%d_%s',outfix,lower(parms.measlist{m}),...
                  layer,fwhm,hemi),...
                'measmat');
            end
          end
        end
      end
    case 'tfmri_stats'
      % conditions are discovered per session rather than listed, since the
      %   set differs by task; the first session with the task supplies them.
      %   Per-run outputs are by slot, r01 the first selected run, as in the
      %   tabulated columns; the per-session file names carry the BIDS run,
      %   which the infix leaves as %s for concat_target to fill per session
      for t = 1:length(parms.tasklist)
        task = parms.tasklist{t};
        conds = find_tfmri_conds(parms,task);
        nslots = size(parms.run_slots.(task),2);
        for c = 1:length(conds)
          for r = 0:nslots
            if r==0, run_in = ''; run_out = ''; else
              run_in = '_run-%s'; run_out = sprintf('_r%02d',r);
            end
            for m = 1:length(parms.measlist)
              for h = 1:length(parms.hemilist)
                hemi = parms.hemilist{h};
                targets = add_target(targets,...
                  sprintf('task-%s%s_hemi-%s_atlas-ico%d_desc-%s_%s',...
                    task,run_in,cmig_get_bids_hemi(hemi),...
                    parms.ico_order,conds{c},parms.measlist{m}),...
                  sprintf('%s_%s_%s%s_%s',task,conds{c},parms.measlist{m},run_out,hemi),...
                  'measmat',task,r);
              end
            end
          end
        end
      end
    case 'fmri_corr'
      % release stems, per the map worked out against the release inventory:
      %   rest  {atlas_}corr{_trimN}
      %   task  {task}_{atlas_}corr_{data|residual|task}{_trim5}
      for t = 1:length(parms.tasklist)
        task = parms.tasklist{t};
        if strcmp(task,'rest')
          % rest carries no task token and no _data infix, and has trim10
          tstem = '';
          variants = {...
            {'',                  'corr'},...
            {'desc-trim5_',       'corr_trim5'},...
            {'desc-trim10_',      'corr_trim10'}};
        else
          tstem = sprintf('%s_',task);
          variants = {...
            {'',                     'corr_data'},...
            {'desc-trim5_',          'corr_data_trim5'},...
            {'desc-residual_',       'corr_residual'},...
            {'desc-residualtrim5_',  'corr_residual_trim5'},...
            {'desc-taskpred_',       'corr_task'},...
            {'desc-taskpredtrim5_',  'corr_task_trim5'}};
        end
        for a = 1:length(parms.atlaslist)
          atlas = parms.atlaslist{a};
          pre = release_atlas_prefix(atlas);
          for v = 1:length(variants)
            % a trimmed corrmat comes from one concatenated, censored and
            %   subsampled series per session, so its run axis is one wide
            if isempty(regexp(variants{v}{1},'trim','once'))
              nruns = [];
            else
              nruns = 1;
            end
            targets = add_target(targets,...
              sprintf('task-%s_atlas-%s_%scorrmat',task,atlas,variants{v}{1}),...
              sprintf('%s%s%s',tstem,pre,variants{v}{2}),[],task,0,nruns);
          end
        end
      end
    case 'fmri_tseries'
      % release stems: {atlas_}tseries for rest, {task}_{atlas_}tseries for the
      %   task families. Kept apart from fmri_corr so the atlas list can differ,
      %   since the time series are far larger than the correlations
      for t = 1:length(parms.tasklist)
        task = parms.tasklist{t};
        if strcmp(task,'rest')
          tstem = '';
        else
          tstem = sprintf('%s_',task);
        end
        for a = 1:length(parms.atlaslist)
          atlas = parms.atlaslist{a};
          % the network atlases have no time series: their correlations are
          %   derived from the ROI-ROI corr files rather than from tseries
          if ismember(atlas,{'gpnet','gpnetaseg'}), continue; end
          targets = add_target(targets,...
            sprintf('task-%s_atlas-%s_tseries',task,atlas),...
            sprintf('%s%stseries',tstem,release_atlas_prefix(atlas)),[],task);
        end
      end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pre = release_atlas_prefix(atlas)
  % the release names every atlas except the gp+dst+dsk+aseg union, whose
  %   files carry no prefix at all. ico3 and ico5 have no release stem, so the
  %   same scheme is extended to them
  switch atlas
    case 'gpdstdskaseg', pre = '';
    case 'dskaseg',      pre = 'dsk_aseg_';
    case 'dstaseg',      pre = 'dst_aseg_';
    case 'gpaseg',       pre = 'gp_aseg_';
    case 'gpnetaseg',    pre = 'gpnet_aseg_';
    case 'gpnet',        pre = 'gpnet_';
    otherwise,           pre = sprintf('%s_',atlas);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = add_target(targets,infix,outstem,varname,task,slot,nruns)
  % task is empty except for the fmri datatypes, where it selects the nruns
  %   and ntpoints parameters. slot is the run slot of a per-run tfmri target,
  %   0 for every other target. nruns, if not empty, replaces the task's
  %   nruns parameter for this target
  if ~exist('task','var'), task = []; end
  if ~exist('slot','var'), slot = 0; end
  if ~exist('nruns','var'), nruns = []; end
  targets(end+1).infix = infix;
  targets(end).outstem = outstem;
  targets(end).varname = varname;
  targets(end).task = task;
  targets(end).slot = slot;
  targets(end).nruns = nruns;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conds = find_tfmri_conds(parms,task)
  % read the condition labels off one session's session-level file names
  %   (the per-run names have run- after task-, so the hemi- glob skips them)
  conds = {};
  for i = 1:length(parms.sessions)
    flist = dir(sprintf('%s/%s_%s_task-%s_hemi-*_desc-*_%s.mat',...
      parms.sessions(i).indir,parms.sessions(i).participant_id,...
      parms.sessions(i).session_id,task,parms.measlist{1}));
    if isempty(flist), continue; end
    for j = 1:length(flist)
      c = regexp(flist(j).name,'_desc-([^_]+)_','tokens','once');
      if ~isempty(c) && ~ismember(c{1},conds), conds{end+1} = c{1}; end %#ok<AGROW>
    end
    break;
  end
  conds = sort(conds);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = get_run_slots(parms)
  % for each task, the BIDS run label in each session's run slots (nsess x
  %   nruns_{task}, empty where the slot is unfilled), read once so the
  %   per-run targets need not reread it
  parms.run_slots = struct();
  for t = 1:length(parms.tasklist)
    task = parms.tasklist{t};
    nruns = task_dims(parms,task);
    slots = cell(length(parms.sessions),0);
    for i = 1:length(parms.sessions)
      runs = session_run_slots(parms,i,task);
      slots(i,1:length(runs)) = runs;
    end
    if ~isempty(nruns), slots(:,end+1:nruns) = {''}; slots = slots(:,1:nruns); end
    slots(cellfun(@isempty,slots)) = {''};
    parms.run_slots.(task) = slots;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runs = session_run_slots(parms,i,task)
  % the run in slot k is the run with slot_num k in the session's
  %   _desc-select_runs.tsv, written by cmig_select_bids_func_runs.py. Without
  %   one, the session's per-run files fill the slots in run order
  runs = {};
  pid = parms.sessions(i).participant_id;
  sid = parms.sessions(i).session_id;
  indir = parms.sessions(i).indir;
  fname = sprintf('%s/%s_%s_task-%s_desc-select_runs.tsv',indir,pid,sid,task);
  if exist(fname,'file')
    fid = fopen(fname,'rt');
    hdr = strsplit(fgetl(fid),'\t','CollapseDelimiters',false);
    i_run = find(strcmp(hdr,'run'));
    i_slot = find(strcmp(hdr,'mr_y_qc__post__series__tfmri__slot_num'));
    if ~isempty(i_run) && ~isempty(i_slot)
      while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        vals = strsplit(tline,'\t','CollapseDelimiters',false);
        if length(vals) < max(i_run,i_slot), continue; end
        k = str2double(vals{i_slot});
        if isfinite(k) && k >= 1, runs{k} = vals{i_run}; end %#ok<AGROW>
      end
      fclose(fid);
      return;
    end
    fclose(fid);
    fprintf('%s: WARNING: no run or slot_num column in %s\n',mfilename,fname);
  end
  flist = dir(sprintf('%s/%s_%s_task-%s_run-*_desc-*_%s.mat',...
    indir,pid,sid,task,parms.measlist{1}));
  for j = 1:length(flist)
    r = regexp(flist(j).name,'_run-([^_]+)_','tokens','once');
    if ~isempty(r) && ~ismember(r{1},runs), runs{end+1} = r{1}; end %#ok<AGROW>
  end
  runs = sort(runs);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname = session_fname(parms,i,infix)
  fname = sprintf('%s/%s_%s_%s.mat',parms.sessions(i).indir,...
    parms.sessions(i).participant_id,parms.sessions(i).session_id,infix);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,result] = concat_target(target,parms,outdir)
  % result is 'written', 'existing' (kept, forceflag 0) or 'nodata' (no
  %   session has an input file). A {outstem}_status.tsv is written beside
  %   the .mat when it is written or the target fails; one kept from an
  %   earlier run is left as it is
  startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');
  [errcode,errmsg,result,tcounts] = concat_target_data(target,parms,outdir);
  if errcode || strcmp(result,'written')
    write_concat_status(sprintf('%s/%s_status.tsv',outdir,target.outstem),...
      errcode,errmsg,startdate,tcounts);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg,result,tcounts] = concat_target_data(target,parms,outdir)
  % one output file per target, holding every variable the per-session files
  %   hold. The session axis is already present in each per-session file as a
  %   singleton, so concatenating grows an existing axis and never reshapes.
  %   The run and time extents come from the parameters rather than from the
  %   sessions, so the output is filled in a single pass, allocating on the
  %   first session that has a value for a variable. tcounts are the target
  %   status columns: sessions with and without an input file, and sessions
  %   padded or truncated on the run and time axes (0 where there is none)
  errcode = 0; errmsg = []; result = 'nodata';
  nsess = length(parms.sessions);
  tcounts = struct('nsess',nsess,'nfound',[],'nmissing',[],...
    'npad_run',[],'ntrunc_run',[],'npad_time',[],'ntrunc_time',[]);
  fname_out = sprintf('%s/%s.mat',outdir,target.outstem);
  if exist(fname_out,'file') && ~parms.forceflag, result = 'existing'; return; end
  [nruns,ntpoints] = task_dims(parms,target.task);
  if ~isempty(target.nruns), nruns = target.nruns; end
  S = struct(); varnames = {}; nvar = 0; nfound = 0;
  rules = {}; dims = []; truncs = {}; rundims = []; timedims = []; present = [];
  npad = []; ntrunc = []; nwarn = [];
  % run-axis padding is warned by session, since all of a session's per-run
  %   variables pad together; padding on other axes is warned by variable.
  %   The status counts sessions on the run and time axes
  padrun = false(1,nsess); truncrun = false(1,nsess);
  padtime = false(1,nsess); trunctime = false(1,nsess);
  for i = 1:nsess
    if target.slot
      % an unfilled slot leaves the session's row as allocated
      run = parms.run_slots.(target.task){i,target.slot};
      if isempty(run), continue; end
      fname = session_fname(parms,i,sprintf(target.infix,run));
    else
      fname = session_fname(parms,i,target.infix);
    end
    if ~exist(fname,'file')
      fprintf('%s: WARNING: missing file %s\n',mfilename,fname);
      continue;
    end
    try
      tmp = load(fname);
    catch
      errcode = 1;
      errmsg = sprintf('corrupted file %s',fname);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    if ~nfound
      % the first session found fixes the variable set and the rules
      varnames = fieldnames(tmp);
      nvar = length(varnames);
      rules = cell(1,nvar); dims = zeros(1,nvar); truncs = cell(1,nvar);
      rundims = zeros(1,nvar); timedims = zeros(1,nvar); present = false(1,nvar);
      npad = zeros(1,nvar); ntrunc = zeros(1,nvar); nwarn = zeros(1,nvar);
      for v = 1:nvar
        [rules{v},dims(v),truncs{v},rundims(v),timedims(v)] = ...
          placement_rule(varnames{v},tmp.(varnames{v}));
      end
      if ~isempty(target.varname) && ~isfield(tmp,target.varname)
        errcode = 1;
        errmsg = sprintf('%s not found in %s',target.varname,fname);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
    elseif ~isequal(sort(fieldnames(tmp)),sort(varnames))
      errcode = 1;
      errmsg = sprintf('variables in %s differ from the first session',fname);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      return;
    end
    nfound = nfound + 1;
    for v = 1:nvar
      varname = varnames{v};
      val = tmp.(varname);
      switch rules{v}
        case 'count'
          S.(varname) = nsess;
        case 'param'
          % these record the extent the output was allocated to, not what any
          %   one session had
          if strcmp(varname,'nruns')
            pval = nruns;
          else
            pval = ntpoints;
          end
          if isempty(pval)
            if ~isfield(S,varname), S.(varname) = val; end
          else
            S.(varname) = pval;
          end
        case {'once','first'}
          if ~isfield(S,varname)
            % first seen, which may be an absent value
            S.(varname) = val;
            present(v) = is_present(val);
          elseif is_present(val)
            if ~present(v)
              S.(varname) = val;
              present(v) = true;
            elseif ~isequal(val,S.(varname))
              if strcmp(rules{v},'once')
                errcode = 1;
                errmsg = sprintf('%s in %s differs from the first session',...
                  varname,fname);
                fprintf('%s: ERROR: %s\n',mfilename,errmsg);
                return;
              end
              nwarn(v) = nwarn(v) + 1;
            end
          end
        case 'stack'
          % a variable the sessions so far left empty is not allocated yet:
          %   an empty value has nothing to place, so nothing is lost
          if isempty(val), continue; end
          if ~isfield(S,varname)
            sz = stack_size(val,dims(v),varname,nsess,nruns,ntpoints);
            S.(varname) = alloc_stack(val,sz);
          end
          [S.(varname),padk,trunck,errcode,errmsg] = place_session(...
            S.(varname),dims(v),truncs{v},i,val,varname,fname);
          if errcode, return; end
          rd = rundims(v); td = timedims(v);
          % padding off the run axis, and truncation, are warned by variable
          other = true(size(padk));
          if rd, other(rd) = false; end
          npad(v) = npad(v) + any(padk(other));
          ntrunc(v) = ntrunc(v) + any(trunck);
          if rd
            padrun(i) = padrun(i) || padk(rd);
            truncrun(i) = truncrun(i) || trunck(rd);
          end
          if td && td <= length(padk)
            padtime(i) = padtime(i) || padk(td);
            trunctime(i) = trunctime(i) || trunck(td);
          end
      end
    end
  end
  if ~nfound
    fprintf('%s: WARNING: no input files for %s, skipping\n',...
      mfilename,target.outstem);
    return;
  end
  % a variable every session left empty is written empty, as the sessions have it
  for v = 1:nvar
    if ~isfield(S,varnames{v}), S.(varnames{v}) = []; end
    if nwarn(v)
      fprintf('%s: WARNING: %s differs from the first session for %d sessions, keeping the first\n',...
        mfilename,varnames{v},nwarn(v));
    end
    if npad(v)
      fprintf('%s: WARNING: %s is smaller than the output for %d sessions, padded\n',...
        mfilename,varnames{v},npad(v));
    end
    if ntrunc(v)
      fprintf('%s: WARNING: %s is larger than the output for %d sessions, truncated\n',...
        mfilename,varnames{v},ntrunc(v));
    end
  end
  if any(padrun)
    fprintf('%s: WARNING: %d sessions have fewer runs than the output for %s, padded\n',...
      mfilename,sum(padrun),target.outstem);
  end
  save(fname_out,'-struct','S','-v7.3','-nocompression');
  fprintf('%s: wrote file %s\n',mfilename,fname_out);
  result = 'written';
  tcounts.nfound = nfound;
  tcounts.nmissing = nsess - nfound;
  tcounts.npad_run = sum(padrun);
  tcounts.ntrunc_run = sum(truncrun);
  tcounts.npad_time = sum(padtime);
  tcounts.ntrunc_time = sum(trunctime);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nruns,ntpoints] = task_dims(parms,task)
  % the run count and the run length differ between tasks, so each has its own
  %   pair of parameters
  nruns = []; ntpoints = [];
  if isempty(task), return; end
  fname = sprintf('nruns_%s',task);
  if isfield(parms,fname), nruns = parms.(fname); end
  fname = sprintf('ntpoints_%s',task);
  if isfield(parms,fname), ntpoints = parms.(fname); end
  if isempty(nruns) || isempty(ntpoints)
    fprintf('%s: WARNING: no nruns or ntpoints parameter for task %s, using the first session\n',...
      mfilename,task);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rule,dim,truncdims,rundim,timedim] = placement_rule(varname,val)
  % which axis of the output the variable grows along, if any, and its run
  %   and time axes (0 if none). nruns and ntpoints come from the parameters the output
  %   was allocated to. nroi and ntpointvec are both 1x1 per session, so the
  %   fixed variables have to be named rather than detected from shape
  fixed_vars = {'roinames','roi1vec','roi2vec','nroi','nnodes',...
                'nbase','nmotion','nreg','nstim','stimregs','stim_labels'};
  % the tseries arrays carry runs on dim 1 and time on dim 3, so their session
  %   axis is dim 2 and both of those extents may be truncated
  dim2_vars = {'datamat_tsdata','datamat_motion','censvec','censvec_deconv',...
               'datamat_design','datamat_duration','datamat_onset'};
  dim = 0; truncdims = []; rundim = 0; timedim = 0;
  if strcmp(varname,'ndirs')
    rule = 'count';
  elseif ismember(varname,{'nruns','ntpoints'})
    rule = 'param';
  elseif ismember(varname,fixed_vars)
    rule = 'once';
  elseif ~(isnumeric(val) || islogical(val))
    rule = 'first';
  elseif ismember(varname,dim2_vars)
    rule = 'stack';
    dim = 2;
    truncdims = [1 3];
    rundim = 1; timedim = 3;
  else
    rule = 'stack';
    dim = 1;
    if ~isempty(regexp(varname,'_run$','once')), truncdims = 2; rundim = 2; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tf = is_present(val)
  % a zero or an empty value means the session has no design, which is not a
  %   disagreement with a session that has one
  if isempty(val)
    tf = false;
  elseif isnumeric(val) || islogical(val)
    tf = any(val(:)~=0);
  else
    tf = true;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sz = stack_size(val,dim,varname,nsess,nruns,ntpoints)
  % the run and time extents come from the parameters, so a session with more
  %   or fewer of either does not set the size of the output
  sz = size(val);
  if ~isempty(nruns) && ~isempty(ntpoints)
    if dim==2
      sz(1) = nruns;
      if length(sz)>=3, sz(3) = ntpoints; end
    elseif length(sz)>=2 && ~isempty(regexp(varname,'_run$','once'))
      sz(2) = nruns;
    end
  end
  sz(dim) = nsess;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = alloc_stack(val,sz)
  % NaN marks a session with no input file, and an extent a session did not
  %   fill. logical and integer have no NaN, so they are zero filled: censvec
  %   is logical in the release files too, so it is kept logical
  if islogical(val)
    out = false(sz);
  elseif isinteger(val)
    out = zeros(sz,class(val));
  else
    out = NaN(sz,'single');
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [datamat,padk,trunck,errcode,errmsg] = place_session(datamat,dim,truncdims,i,val,varname,fname)
  % the extent of every dimension but the session axis comes from the session,
  %   so a shorter session fills its own extent and leaves the rest allocated.
  %   A longer session is truncated on the run and time axes and is an error
  %   anywhere else, since truncating a design would drop regressors whose
  %   identity comes from the contrast set. padk and trunck flag, by
  %   dimension, where the session was padded or truncated
  errcode = 0; errmsg = [];
  padk = false(1,ndims(datamat)); trunck = false(1,ndims(datamat));
  n = ndims(datamat);
  sz = size(val);
  sz(end+1:n) = 1;
  idx = cell(1,n); vidx = cell(1,n);
  for k = 1:n
    if k==dim
      idx{k} = i; vidx{k} = 1;
      continue;
    end
    m = size(datamat,k);
    if sz(k) > m
      if ~ismember(k,truncdims)
        errcode = 1;
        errmsg = sprintf('%s in %s is %d on dim %d against %d in the output',...
          varname,fname,sz(k),k,m);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
      sz(k) = m;
      trunck(k) = true;
    elseif sz(k) < m
      padk(k) = true;
    end
    idx{k} = 1:sz(k); vidx{k} = 1:sz(k);
  end
  datamat(idx{:}) = val(vidx{:});
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = write_vol_info(parms,outdir)
  % written before any measure file, and read back by every later call, so all
  %   the files in outdir share one session axis however many jobs write them.
  %   parms.scan is the session list of a tree scan, if one was made
  errcode = 0; errmsg = [];
  parms.vol_info_written = false;
  parms.scan = [];
  fname_out = sprintf('%s/vol_info.mat',outdir);
  if exist(fname_out,'file') && ~parms.forceflag
    fprintf('%s: reading sessions from %s...\n',mfilename,fname_out);
    vol_info = load(fname_out);
    parms.sessions = struct('participant_id',vol_info.participant_id,...
                            'session_id',vol_info.session_id,...
                            'indir',vol_info.dirlist);
    % no data now is a no-data skip, but a vol_info this run did not write
    %   is left as it is
    if ~isempty(parms.sessions) && ~has_data(parms), parms.sessions = []; end
    return;
  end
  if isempty(parms.sessions)
    [parms.sessions,errcode,errmsg] = get_session_list(parms);
    if errcode, return; end
    parms.scan = parms.sessions;
  end
  % sessions none of which has one of the datatype's files are no data
  if ~isempty(parms.sessions) && ~has_data(parms), parms.sessions = []; end
  nsess = length(parms.sessions);
  if ~nsess
    % a vol_info left by an earlier run is stale (forceflag): removed,
    %   unless other files of the folder use it
    if exist(fname_out,'file'), remove_vol_info(outdir); end
    return;
  end
  mmil_mkdir(outdir);
  participant_id = cell(nsess,1);
  session_id     = cell(nsess,1);
  dirlist        = cell(nsess,1);
  for i = 1:nsess
    participant_id{i} = parms.sessions(i).participant_id;
    session_id{i}     = parms.sessions(i).session_id;
    dirlist{i}        = parms.sessions(i).indir;
  end
  % previously, visitidvec was the MMPS VisitID, site code and year label,
  %   which the BIDS pipeline no longer produces. Kept as the participant and
  %   session joined, for a single unique identifier per row
  visitidvec = strcat(participant_id,'_',session_id);
  % the cortex masks come from fsaverage and do not vary by session, so they
  %   are stored once, as in the release files
  vol_mask_lh = []; vol_mask_rh = [];
  if parms.ico_order > 0
    ico_nverts = [42 162 642 2562 10242 40962 163842];
    nverts = ico_nverts(parms.ico_order);
    labelpath = sprintf('%s/subjects/fsaverage/label',getenv('FREESURFER_HOME'));
    if exist(labelpath,'dir')
      verts_lh = fs_read_label(sprintf('%s/lh.cortex.label',labelpath));
      verts_rh = fs_read_label(sprintf('%s/rh.cortex.label',labelpath));
    else
      fprintf('%s: WARNING: label path %s not found\n',mfilename,labelpath);
      verts_lh = 1:nverts; verts_rh = 1:nverts;
    end
    vol_mask_lh = zeros(nverts,1); vol_mask_lh(intersect(verts_lh,1:nverts)) = 1;
    vol_mask_rh = zeros(nverts,1); vol_mask_rh(intersect(verts_rh,1:nverts)) = 1;
  end
  % with forceflag, a vol_info with the same rows is kept rather than
  %   rewritten, so that the jobs sharing it, and fmri_corr and fmri_tseries,
  %   do not write it at once
  if exist(fname_out,'file')
    try
      old = load(fname_out);
      same = isequal(old.participant_id,participant_id) &&...
             isequal(old.session_id,session_id) &&...
             isequal(old.dirlist,dirlist) &&...
             isequal(old.vol_mask_lh,vol_mask_lh) &&...
             isequal(old.vol_mask_rh,vol_mask_rh);
    catch
      same = false;
    end
    if same
      fprintf('%s: kept %s, which has the same sessions\n',mfilename,fname_out);
      return;
    end
  end
  save(fname_out,'participant_id','session_id','visitidvec','dirlist',...
    'vol_mask_lh','vol_mask_rh','-v7.3','-nocompression');
  fprintf('%s: wrote file %s\n',mfilename,fname_out);
  parms.vol_info_written = true;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tf = has_data(parms)
  % true at the first session with a file matching one of data_patterns
  tf = false;
  for i = 1:length(parms.sessions)
    for p = 1:length(parms.data_patterns)
      flist = dir(sprintf('%s/%s_%s%s',parms.sessions(i).indir,...
        parms.sessions(i).participant_id,parms.sessions(i).session_id,...
        parms.data_patterns{p}));
      if ~isempty(flist), tf = true; return; end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function remove_vol_info(outdir)
  % remove vol_info.mat, unless another .mat in outdir uses its rows (another
  %   datatype's, as fmri_corr and fmri_tseries share a corrmat folder, or
  %   one kept from an earlier run), then outdir and its parent if that
  %   leaves them empty (rmdir with an output argument does not error on a
  %   folder that is not empty)
  fname = sprintf('%s/vol_info.mat',outdir);
  flist = dir(sprintf('%s/*.mat',outdir));
  if any(~strcmp({flist.name},'vol_info.mat'))
    fprintf('%s: kept file %s, which other files in %s use\n',...
      mfilename,fname,outdir);
    return;
  end
  if exist(fname,'file')
    delete(fname);
    fprintf('%s: removed file %s\n',mfilename,fname);
  end
  [~] = rmdir(outdir);
  [~] = rmdir(fileparts(outdir));
return;
