function [outdir,errcode,errmsg] = cmig_postqc_bids_anat(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_postqc_bids_anat(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Automated post-processing QC for BIDS anat derivatives
%   Processes T1w and T2w scans; checks registration to T1w (T2w only)
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%     processed data at {proc_rootdir}/{participant_id}/{session_id}/anat
%     output from cmig_proc_bids_anat
%   fsurf_rootdir: BIDS FreeSurfer derivatives root directory
%     FreeSurfer recons at {fsurf_rootdir}/{participant_id}_{session_id}
%     may be empty to skip FreeSurfer-dependent checks
%   output_rootdir: BIDS output root directory for QC results
%
% Optional Input ('key',value):
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'verbose': [0|1] display status messages and warnings
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     form: {output_rootdir}/{participant_id}/{session_id}/anat
%     one session TSV and one series TSV per scan type:
%       {participant_id}_{session_id}_desc-session_postqc.tsv
%       {participant_id}_{session_id}_desc-{T1w,T2w}_postqc.tsv
%     the series TSVs have one row per scan
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/11/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 08/07/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/29/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Last Mod: 11/15/2026 by Don Hagler
%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,5), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin);
parms.status_extra = [];

[parms,errcode,errmsg] = check_input(parms);
if errcode
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-postqc',errcode,errmsg,parms.status_extra,startdate);
  return;
end

[parms,procflag] = check_output(parms);

if procflag
  [series_qcstruct,t1w_session,t2w_session,errcode,errmsg] = qc_anat(parms);
  if ~errcode
    combined_row = make_anat_combined_row(t1w_session,t2w_session,...
      parms.participant_id,parms.session_id);
    fname_combined = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);
    cmig_write_tabdat_postqc({combined_row},fname_combined,parms.forceflag);
    % the desc-{suffix} files are series-level, one row per scan, as
    %   desc-bold is for func and desc-dwi for dwi; the per-modality
    %   aggregates feed the session row only
    stypes = {'T1w','T2w'};
    for s=1:length(stypes)
      stype = stypes{s};
      stype_struct = series_qcstruct(cellfun(@(r) strcmp(r.modality,stype), series_qcstruct));
      if isempty(stype_struct), continue; end
      fname_series = sprintf('%s/%s_desc-%s_postqc.tsv',parms.outdir,parms.outstem,stype);
      cmig_write_tabdat_postqc(stype_struct,fname_series,parms.forceflag);
    end
  end
end

if parms.cleanup_flag
  cleanup_errmsg = cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir);
  if isempty(errmsg), errmsg = cleanup_errmsg; end
end
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-postqc',errcode,errmsg,parms.status_extra,startdate,procflag);

outdir = parms.outdir;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,options)
  options = cmig_parse_deployed_args(options);
  parms = mmil_args2parms(options,{...
    'participant_id',participant_id,[],...
    'session_id',session_id,[],...
    'proc_rootdir',proc_rootdir,[],...
    'fsurf_rootdir',fsurf_rootdir,[],...
    'output_rootdir',output_rootdir,[],...
  ...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });

  parms.outdir  = sprintf('%s/%s/%s/anat',...
    parms.output_rootdir,parms.participant_id,parms.session_id);
  parms.tmpdir  = sprintf('%s/autoQC',parms.outdir);
  parms.outstem = sprintf('%s_%s',parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];
  pid = parms.participant_id;
  sid = parms.session_id;

  fprintf('%s: checking input for %s_%s (modality=anat)\n',mfilename,pid,sid);

  parms.fnames_t1w = []; parms.fnames_info_t1w = [];
  [parms.fnames_t1w,parms.fnames_info_t1w] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'anat','T1w');

  parms.fnames_t2w = []; parms.fnames_info_t2w = [];
  [parms.fnames_t2w,parms.fnames_info_t2w] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'anat','T2w');

  if isempty(parms.fsurf_rootdir)
    parms.fspath = [];
  else
    [parms.fspath,parms.fs_version,parms.fs_build,parms.fs_complete_flag,errcode,errmsg] = ...
      cmig_check_bids_fsurf(parms.participant_id,parms.session_id,parms.fsurf_rootdir);
    if errcode
      fprintf('%s: WARNING: %s\n',mfilename,errmsg);
      parms.fspath = [];
      errcode = 0;
    else
      parms.status_extra = struct('fs_version',parms.fs_version,...
        'fs_build',parms.fs_build,'fs_complete_flag',parms.fs_complete_flag);
    end
  end

  for f=1:length(parms.fnames_t1w)
    fprintf('%s: T1w processed: %s\n',mfilename,parms.fnames_t1w{f});
  end
  for f=1:length(parms.fnames_t2w)
    fprintf('%s: T2w processed: %s\n',mfilename,parms.fnames_t2w{f});
  end

  if isempty(parms.fnames_t1w) && isempty(parms.fnames_t2w)
    errcode = 1;
    errmsg = 'no processed anat data found';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;
  output_files = {};

  output_files{end+1} = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);
  if ~isempty(parms.fnames_t1w)
    output_files{end+1} = sprintf('%s/%s_desc-T1w_postqc.tsv',parms.outdir,parms.outstem);
  end
  if ~isempty(parms.fnames_t2w)
    output_files{end+1} = sprintf('%s/%s_desc-T2w_postqc.tsv',parms.outdir,parms.outstem);
  end

  if ~exist(parms.outdir,'dir')
    procflag = 1;
  elseif parms.forceflag
    procflag = 1;
  else
    for f=1:length(output_files)
      if ~exist(output_files{f},'file')
        [~,name,ext] = fileparts(output_files{f});
        fprintf('%s: missing output file %s\n',mfilename,[name ext]);
        procflag = 1;
        break;
      end
    end
  end

  if procflag
    mmil_mkdir(parms.tmpdir);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [series_qcstruct,t1w_session,t2w_session,errcode,errmsg] = qc_anat(parms)
  series_qcstruct = {}; t1w_session = []; t2w_session = [];
  errcode = 0; errmsg = [];

  n = length(parms.fnames_t1w) + length(parms.fnames_t2w);
  errcodes = zeros(n,1);
  errmsgs  = cell(n,1);
  k = 0;

  t1w_rows = {};
  for f=1:length(parms.fnames_t1w)
    k = k + 1;
    fname = parms.fnames_t1w{f};
    fname_info = parms.fnames_info_t1w{f};
    [qcinfo,errcodes(k),errmsgs{k}] = qc_smri_scan(fname,'T1w',fname_info,parms);
    if ~isempty(qcinfo)
      row = convert_smri_row(qcinfo,parms.participant_id,parms.session_id);
      series_qcstruct{end+1} = rename_series_cols(row);
      t1w_rows{end+1} = row;
    end
  end

  t2w_rows = {};
  for f=1:length(parms.fnames_t2w)
    k = k + 1;
    fname = parms.fnames_t2w{f};
    fname_info = parms.fnames_info_t2w{f};
    [qcinfo,errcodes(k),errmsgs{k}] = qc_smri_scan(fname,'T2w',fname_info,parms);
    if ~isempty(qcinfo)
      row = convert_smri_row(qcinfo,parms.participant_id,parms.session_id);
      series_qcstruct{end+1} = rename_series_cols(row);
      t2w_rows{end+1} = row;
    end
  end

  if ~isempty(t1w_rows)
    t1w_session = make_smri_session_row(t1w_rows,parms.participant_id,parms.session_id);
  end
  if ~isempty(t2w_rows)
    t2w_session = make_smri_session_row(t2w_rows,parms.participant_id,parms.session_id);
  end

  nonempty_msgs = errmsgs(~cellfun(@isempty,errmsgs));
  if ~isempty(nonempty_msgs)
    errmsg = strtrim(sprintf('%s ',nonempty_msgs{:}));
  end
  if all(errcodes), errcode = 1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auto_qcinfo,errcode,errmsg] = qc_smri_scan(fname,scantype,fname_info,parms)
  auto_qcinfo = [];
  errcode = 0; errmsg = [];

  if parms.verbose
    fprintf('%s: running QC for %s scan: %s\n',mfilename,scantype,fname);
  end

  auto_qcinfo = init_smri_qcinfo(fname,scantype,parms);
  [auto_qcinfo,errcode,errmsg] = read_proc_json(auto_qcinfo,fname,fname_info,parms);
  if errcode, return; end

  if strcmp(scantype,'T2w') && ~isempty(parms.fnames_t1w)
    [auto_qcinfo,errcode,errmsg] = check_regT1_t2w(auto_qcinfo,fname,parms);
    if errcode, return; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function auto_qcinfo = init_smri_qcinfo(fname,scantype,parms)
  auto_qcinfo = [];
  [~,fstem] = cmig_get_fileparts(fname);
  auto_qcinfo.VisitID                = sprintf('%s_%s',parms.participant_id,parms.session_id);
  auto_qcinfo.stype                  = scantype;
  auto_qcinfo.fstem                  = fstem;
  auto_qcinfo.StudyDate              = '';
  auto_qcinfo.Manufacturer           = '';
  auto_qcinfo.ManufacturersModelName = '';
  auto_qcinfo.DeviceSerialNumber     = '';
  if strcmp(scantype,'T1w')
    auto_qcinfo.atlas_reg_failed     = NaN;
    auto_qcinfo.atlas_reg_corr       = NaN;
    auto_qcinfo.atlas_reg_fvox_miss  = NaN;
    auto_qcinfo.atlas_reg_corr_fail      = NaN;
    auto_qcinfo.atlas_reg_fvox_miss_fail = NaN;
    auto_qcinfo.atlas_reg_rotat      = NaN;
    auto_qcinfo.atlas_reg_trans      = NaN;
  end
  if strcmp(scantype,'T2w')
    auto_qcinfo.regT1_corr     = NaN;
    auto_qcinfo.regT1_norm     = NaN;
    auto_qcinfo.regT1_warp     = NaN;
    auto_qcinfo.regT1_imgcorr1 = NaN;
    auto_qcinfo.regT1_imgcorr2 = NaN;
    auto_qcinfo.regT1_imgdist3 = NaN;
    auto_qcinfo.regT1_imgdist4 = NaN;
    auto_qcinfo.regT1_cost0    = NaN;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auto_qcinfo,errcode,errmsg] = read_proc_json(auto_qcinfo,fname,fname_info,parms)
  errcode = 0; errmsg = [];
  if isempty(fname_info) || ~exist(fname_info,'file')
    fname_info = regexprep(fname,'\.nii.*','.json');
  end
  if ~exist(fname_info,'file')
    errmsg = sprintf('proc JSON sidecar not found: %s',fname_info);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end
  try
    json = jsondecode(fileread(fname_info));
  catch me
    errmsg = sprintf('failed to read JSON %s: %s',fname_info,me.message);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end
  auto_qcinfo.StudyDate              = mmil_getfield(json,'AcquisitionDateTime','');
  auto_qcinfo.Manufacturer           = mmil_getfield(json,'Manufacturer','');
  auto_qcinfo.ManufacturersModelName = mmil_getfield(json,'ManufacturersModelName','');
  auto_qcinfo.DeviceSerialNumber     = mmil_getfield(json,'DeviceSerialNumber','');
  if strcmp(auto_qcinfo.stype,'T1w')
    auto_qcinfo.atlas_reg_failed     = mmil_getfield(json,'AtlasRegistrationFailed',NaN);
    auto_qcinfo.atlas_reg_corr       = mmil_getfield(json,'AtlasRegistrationCorrelation',NaN);
    auto_qcinfo.atlas_reg_fvox_miss  = mmil_getfield(json,'AtlasRegistrationFracVoxMissing',NaN);
    auto_qcinfo.atlas_reg_corr_fail      = get_num(json,'AtlasRegistrationFailedCorrelation');
    auto_qcinfo.atlas_reg_fvox_miss_fail = get_num(json,'AtlasRegistrationFailedFracVoxMissing');
    % total rotation and translation of the atlas registration applied
    %   (for the fallback, the offset that centers the volume)
    M = cmig_get_reg_T1w(json);
    if isnumeric(M) && isequal(size(M),[4,4])
      [U,~,V] = svd(M(1:3,1:3));
      R = U*V';
      auto_qcinfo.atlas_reg_rotat = acosd(min(1,max(-1,(trace(R)-1)/2)));
      auto_qcinfo.atlas_reg_trans = norm(M(1:3,4));
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = get_num(json,fld)
  % a numeric sidecar field; NaN if missing or null (jsonencode writes NaN as null)
  val = mmil_getfield(json,fld,NaN);
  if isempty(val) || ~isnumeric(val), val = NaN; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [auto_qcinfo,errcode,errmsg] = check_regT1_t2w(auto_qcinfo,fname_t2w,parms)
  errcode = 0; errmsg = [];
  checkname = 'T2w registration to T1w';
  [~,fstem_t2w] = cmig_get_fileparts(fname_t2w);

  fname_t1w = [];
  for f=1:length(parms.fnames_t1w)
    if exist(parms.fnames_t1w{f},'file')
      fname_t1w = parms.fnames_t1w{f};
      break;
    end
  end
  if isempty(fname_t1w)
    errmsg = sprintf('no T1w proc file found, skipping %s',checkname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end

  fname_out = sprintf('%s/%s_regT1.mat',parms.tmpdir,fstem_t2w);
  if ~exist(fname_out,'file') || parms.forceflag
    fname_t2w_tmp = cmig_import_to_mgz(fname_t2w,parms.tmpdir,[],parms.forceflag);
    if isempty(fname_t2w_tmp)
      errmsg = sprintf('failed to import T2w, skipping %s',checkname);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
    fname_t1w_tmp = cmig_import_to_mgz(fname_t1w,parms.tmpdir,[],parms.forceflag);
    if isempty(fname_t1w_tmp)
      errmsg = sprintf('T1w tmp file missing, skipping %s',checkname);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end

    % mri_segreg boundary cost of the registration, given a FreeSurfer recon
    regT1_cost0 = NaN;
    if ~isempty(parms.fspath)
      [regT1_cost0,ec,em] = cmig_check_regT1_segreg(fname_t2w_tmp,parms.fspath,...
        'contrast','t2','tmpdir',parms.tmpdir,...
        'verbose',parms.verbose,'forceflag',parms.forceflag);
      if ec
        fprintf('%s: WARNING: no boundary cost for %s: %s\n',mfilename,checkname,em);
      end
    end

    fname_mask = [];
    if ~isempty(parms.fspath)
      fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
      if exist(fname_aseg,'file')
        fname_cc = cmig_create_aseg_ccmask(parms.fspath,parms.tmpdir,parms.verbose,parms.forceflag);
        if ~isempty(fname_cc) && exist(fname_cc,'file')
          fname_mask = cmig_resample_mask_to_proc(fname_cc,fname_t1w_tmp,parms.tmpdir,parms.verbose,parms.forceflag);
        end
      end
    end

    % without a CC mask, whole-volume correlations are not comparable to
    %   the CC-masked values, so record NaN metrics instead
    if isempty(fname_mask) || ~exist(fname_mask,'file')
      fprintf('%s: WARNING: no CC mask for %s, setting regT1 metrics to NaN\n',...
        mfilename,checkname);
      regT1_corr = NaN; regT1_norm = NaN; regT1_warp = NaN; results = [];
      save(fname_out,'regT1_corr','regT1_norm','regT1_warp','regT1_cost0','results');
    else
      if parms.verbose
        fprintf('%s: checking T2w registration to T1w...\n',mfilename);
        tic;
      end
      p = [];
      entities_t2w = regexprep(rm_desc_preproc(fstem_t2w),'_T2w$','');
      p.fname_plot = sprintf('%s/%s_desc-T2w_regT1.jpg',parms.outdir,entities_t2w);
      p.forceflag = parms.forceflag;
      p.fname_mask = fname_mask;
      a = mmil_parms2args(p);
      try
        [regT1_corr,regT1_norm,regT1_warp,results] = ...
          mmil_check_regT1(fname_t2w_tmp,fname_t1w_tmp,a{:});
        if parms.verbose, toc; end
        save(fname_out,'regT1_corr','regT1_norm','regT1_warp','regT1_cost0','results');
      catch me
        errmsg = sprintf('mmil_check_regT1 failed for T2w: %s',me.message);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        errcode = 1;
        return;
      end
    end
  else
    tmp = load(fname_out);
    regT1_corr = tmp.regT1_corr;
    regT1_norm = tmp.regT1_norm;
    regT1_warp = tmp.regT1_warp;
    regT1_cost0 = mmil_getfield(tmp,'regT1_cost0',NaN);
    results    = tmp.results;
  end

  auto_qcinfo.regT1_corr = regT1_corr;
  auto_qcinfo.regT1_norm = regT1_norm;
  auto_qcinfo.regT1_warp = regT1_warp;
  auto_qcinfo.regT1_cost0 = regT1_cost0;
  if isstruct(results)
    auto_qcinfo.regT1_imgcorr1 = mmil_getfield(results,'ImgCorr1',NaN);
    auto_qcinfo.regT1_imgcorr2 = mmil_getfield(results,'ImgCorr2',NaN);
    auto_qcinfo.regT1_imgdist3 = mmil_getfield(results,'ImgDist3',NaN);
    auto_qcinfo.regT1_imgdist4 = mmil_getfield(results,'ImgDist4',NaN);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = convert_smri_row(auto_qcinfo, participant_id, session_id)
  row.participant_id = participant_id;
  row.session_id     = session_id;
  % series-identifying entity columns
  row.modality = auto_qcinfo.stype;
  fstem = rm_desc_preproc(auto_qcinfo.fstem);
  row.acq_time = auto_qcinfo.StudyDate;
  row.filename = sprintf('anat/%s.nii.gz', fstem);
  tok = regexp(fstem,'(?:^|_)run-([^_]+)','tokens','once');
  if ~isempty(tok), row.run = tok{1}; else row.run = 'n/a'; end
  tok = regexp(fstem,'(?:^|_)acq-([^_]+)','tokens','once');
  if ~isempty(tok), row.acq = tok{1}; else row.acq = 'n/a'; end
  tok = regexp(fstem,'(?:^|_)rec-([^_]+)','tokens','once');
  if ~isempty(tok), row.rec = tok{1}; else row.rec = 'n/a'; end
  % atlas registration (T1w only)
  if strcmp(auto_qcinfo.stype,'T1w')
    row.mr_y_qc__post__aut__smri__atlas__reg__fail_indicator  = auto_qcinfo.atlas_reg_failed;
    row.mr_y_qc__post__aut__smri__atlas__reg__corr_score      = auto_qcinfo.atlas_reg_corr;
    row.mr_y_qc__post__aut__smri__atlas__reg__vox__miss_prop = auto_qcinfo.atlas_reg_fvox_miss;
    row.mr_y_qc__post__aut__smri__atlas__reg__fail__corr_score = auto_qcinfo.atlas_reg_corr_fail;
    row.mr_y_qc__post__aut__smri__atlas__reg__fail__vox__miss_prop = auto_qcinfo.atlas_reg_fvox_miss_fail;
    row.mr_y_qc__post__aut__smri__atlas__reg_rotat = auto_qcinfo.atlas_reg_rotat;
    row.mr_y_qc__post__aut__smri__atlas__reg_trans = auto_qcinfo.atlas_reg_trans;
  end
  % T2w-only registration metrics
  if isfield(auto_qcinfo,'regT1_corr')
    row.mr_y_qc__post__aut__smri__t2__imgcorr_score  = auto_qcinfo.regT1_corr;
    row.mr_y_qc__post__aut__smri__t2__rigid_score     = auto_qcinfo.regT1_norm;
    row.mr_y_qc__post__aut__smri__t2__warp_score      = auto_qcinfo.regT1_warp;
    row.mr_y_qc__post__aut__smri__t2__imgcorr1_score  = auto_qcinfo.regT1_imgcorr1;
    row.mr_y_qc__post__aut__smri__t2__imgcorr2_score  = auto_qcinfo.regT1_imgcorr2;
    row.mr_y_qc__post__aut__smri__t2__imgdist1_score  = auto_qcinfo.regT1_imgdist3;
    row.mr_y_qc__post__aut__smri__t2__imgdist2_score  = auto_qcinfo.regT1_imgdist4;
    row.mr_y_qc__post__aut__smri__t2__segreg__cost_score = auto_qcinfo.regT1_cost0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function series_row = rename_series_cols(row)
  % series columns are mr_y_qc__post__series__{measure}, with no modality;
  %   the session row is aggregated from the unrenamed rows
  series_row = struct();
  flds = fieldnames(row);
  for k=1:length(flds)
    fld = regexprep(flds{k},'^mr_y_qc__post__aut__smri__t2__(img|segreg__)',...
      'mr_y_qc__post__series__regt1__$1');
    fld = regexprep(fld,'^mr_y_qc__post__aut__smri__(t2__)?','mr_y_qc__post__series__');
    series_row.(fld) = row.(flds{k});
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = make_smri_session_row(rows,participant_id,session_id)
  session_row.participant_id = participant_id;
  session_row.session_id     = session_id;
  skip = {'participant_id','session_id','modality','run','acq','rec','acq_time','filename'};
  flds = fieldnames(rows{1});
  for k=1:length(flds)
    if ~ismember(flds{k},skip)
      session_row.(flds{k}) = rows{1}.(flds{k});
    end
  end
  if length(rows) > 1
    for k=1:length(flds)
      fld = flds{k};
      if ismember(fld,skip), continue; end
      vals = cellfun(@(r) mmil_getfield(r,fld,NaN), rows);
      if isnumeric(vals) && ~all(isnan(vals))
        good = vals(~isnan(vals));
        if ~isempty(regexp(fld,'_indicator$','once'))
          session_row.(fld) = max(good);
        else
          session_row.(fld) = mean(good);
        end
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_row = make_anat_combined_row(t1w_session,t2w_session,participant_id,session_id)
  combined_row.participant_id = participant_id;
  combined_row.session_id     = session_id;
  skip = {'participant_id','session_id'};
  t1w_fields = {'mr_y_qc__post__aut__smri__atlas__reg__fail_indicator',...
                'mr_y_qc__post__aut__smri__atlas__reg__corr_score',...
                'mr_y_qc__post__aut__smri__atlas__reg__vox__miss_prop',...
                'mr_y_qc__post__aut__smri__atlas__reg__fail__corr_score',...
                'mr_y_qc__post__aut__smri__atlas__reg__fail__vox__miss_prop',...
                'mr_y_qc__post__aut__smri__atlas__reg_rotat',...
                'mr_y_qc__post__aut__smri__atlas__reg_trans'};
  t2w_fields = {'mr_y_qc__post__aut__smri__t2__imgcorr_score',...
                'mr_y_qc__post__aut__smri__t2__rigid_score',...
                'mr_y_qc__post__aut__smri__t2__warp_score',...
                'mr_y_qc__post__aut__smri__t2__imgcorr1_score',...
                'mr_y_qc__post__aut__smri__t2__imgcorr2_score',...
                'mr_y_qc__post__aut__smri__t2__imgdist1_score',...
                'mr_y_qc__post__aut__smri__t2__imgdist2_score',...
                'mr_y_qc__post__aut__smri__t2__segreg__cost_score'};
  if ~isempty(t1w_session)
    flds = fieldnames(t1w_session);
    for k=1:length(flds)
      if ~ismember(flds{k},skip), combined_row.(flds{k}) = t1w_session.(flds{k}); end
    end
  else
    for k=1:length(t1w_fields), combined_row.(t1w_fields{k}) = NaN; end
  end
  if ~isempty(t2w_session)
    flds = fieldnames(t2w_session);
    for k=1:length(flds)
      if ~ismember(flds{k},skip), combined_row.(flds{k}) = t2w_session.(flds{k}); end
    end
  else
    for k=1:length(t2w_fields), combined_row.(t2w_fields{k}) = NaN; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = rm_desc_preproc(fstem)
  % the rawdata name of a proc file, which adds desc-preproc: the filename
  %   column joins with raw QC on it, and the jpg names carry the run entities
  fstem = regexprep(fstem,'_desc-preproc(?=_[^_]+$)','');
return;
