function [outdir,errcode,errmsg] = cmig_postqc_bids_fmap(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_postqc_bids_fmap(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Automated post-processing QC for BIDS fmap derivatives
%   Computes stdev/tSNR separately for func and dwi fmaps, AP and PA;
%   calculates AP-PA mean difference and per-scan RMS from inter-scan mean;
%   calculates per-scan correlation and centre of mass distance to the
%   other scans of the same direction;
%   calculates RMS of the B0 unwarp displacement fields (in mm);
%   computes brain cutoff for the field map scans;
%   checks the registration of each field map to the T1w
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%     processed fmap EPI data at {proc_rootdir}/{participant_id}/{session_id}/fmap
%     output from cmig_proc_bids_func
%   fsurf_rootdir: BIDS FreeSurfer derivatives root directory
%     used for the brain cutoff check and the CC mask and boundary cost of
%     the registration check
%     may be empty to skip FreeSurfer-dependent checks
%   output_rootdir: BIDS output root directory for QC results
%
% Optional Input ('key',value):
%   'acq_types': cell array of acquisition types to process
%     {default = {'func','dwi'}}
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'verbose': [0|1] display status messages and warnings
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     form: {output_rootdir}/{participant_id}/{session_id}/fmap
%     one session TSV and one series TSV per acq type:
%       {participant_id}_{session_id}_acq-{acq}_desc-session_postqc.tsv
%       {participant_id}_{session_id}_acq-{acq}_desc-epi_postqc.tsv
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  06/17/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 07/22/2026 by Don Hagler
% Prev Mod: 08/05/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/29/2026 by Don Hagler
% Prev Mod: 10/01/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 10/08/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 11/02/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Last Mod: 11/15/2026 by Don Hagler
%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,5), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin);

[parms,errcode,errmsg] = check_input(parms);
if errcode
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,parms.status_module,errcode,errmsg,parms.status_extra,startdate);
  return;
end

[parms,procflag] = check_output(parms);

if procflag
  [series_qcstruct,session_rows,errcode,errmsg] = qc_fmap(parms);
  if ~errcode
    for k=1:length(parms.acq_types)
      acq = parms.acq_types{k};
      acq_struct = series_qcstruct(cellfun(@(r) strcmp(r.acq,acq), series_qcstruct));
      if ~isempty(acq_struct)
        fname_series = sprintf('%s/%s_acq-%s_desc-epi_postqc.tsv',parms.outdir,parms.outstem,acq);
        cmig_write_tabdat_postqc(acq_struct,fname_series,parms.forceflag);
      end
      if isfield(session_rows,acq) && ~isempty(session_rows.(acq))
        fname_session = sprintf('%s/%s_acq-%s_desc-session_postqc.tsv',parms.outdir,parms.outstem,acq);
        cmig_write_tabdat_postqc({session_rows.(acq)},fname_session,parms.forceflag);
      end
    end
  end
end

if parms.cleanup_flag
  cleanup_errmsg = cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir);
  if isempty(errmsg), errmsg = cleanup_errmsg; end
end
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,parms.status_module,errcode,errmsg,parms.status_extra,startdate,procflag);

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
    'acq_types',{'func','dwi'},[],...
    'cleanup_flag',true,[false true],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });

  if ischar(parms.acq_types)
    parms.acq_types = strsplit(strtrim(parms.acq_types));
  end

  parms.outdir  = sprintf('%s/%s/%s/fmap',...
    parms.output_rootdir,parms.participant_id,parms.session_id);
  parms.tmpdir  = sprintf('%s/autoQC',parms.outdir);
  parms.outstem = sprintf('%s_%s',parms.participant_id,parms.session_id);

  if length(parms.acq_types) == 1
    parms.status_module = sprintf('acq-%s_desc-postqc',parms.acq_types{1});
  else
    parms.status_module = 'desc-postqc';
  end

  parms.status_extra = [];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];
  pid = parms.participant_id;
  sid = parms.session_id;

  fprintf('%s: checking input for %s_%s (modality=fmap)\n',mfilename,pid,sid);

  [fnames_fmap,fnames_info_fmap] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'fmap','epi');

  parms.fnames_func_ap = {}; parms.fnames_info_func_ap = {};
  parms.fnames_func_pa = {}; parms.fnames_info_func_pa = {};
  parms.fnames_dwi_ap  = {}; parms.fnames_info_dwi_ap  = {};
  parms.fnames_dwi_pa  = {}; parms.fnames_info_dwi_pa  = {};

  for f=1:length(fnames_fmap)
    fname      = fnames_fmap{f};
    fname_info = fnames_info_fmap{f};

    [~,fstem] = cmig_get_fileparts(fname);
    if isempty(regexp(fstem,'_epi$','once')), continue; end
    tok = regexp(fstem,'(?:^|_)acq-([^_]+)','tokens','once');
    if isempty(tok), continue; end
    acq_type = lower(tok{1});
    if ~ismember(acq_type,parms.acq_types), continue; end

    direction = get_pe_direction(fname,fname_info);
    if isempty(direction), continue; end

    fprintf('%s: fmap %s %s: %s\n',mfilename,acq_type,direction,fname);

    switch acq_type
      case 'func'
        if strcmp(direction,'AP')
          parms.fnames_func_ap{end+1}      = fname;
          parms.fnames_info_func_ap{end+1} = fname_info;
        else
          parms.fnames_func_pa{end+1}      = fname;
          parms.fnames_info_func_pa{end+1} = fname_info;
        end
      case 'dwi'
        if strcmp(direction,'AP')
          parms.fnames_dwi_ap{end+1}      = fname;
          parms.fnames_info_dwi_ap{end+1} = fname_info;
        else
          parms.fnames_dwi_pa{end+1}      = fname;
          parms.fnames_info_dwi_pa{end+1} = fname_info;
        end
    end
  end

  parms.acq_types_present = {};
  for k=1:length(parms.acq_types)
    acq    = parms.acq_types{k};
    fld_ap = sprintf('fnames_%s_ap',acq);
    fld_pa = sprintf('fnames_%s_pa',acq);
    has_ap = isfield(parms,fld_ap) && ~isempty(parms.(fld_ap));
    has_pa = isfield(parms,fld_pa) && ~isempty(parms.(fld_pa));
    if has_ap || has_pa
      parms.acq_types_present{end+1} = acq;
    end
  end

  all_empty = isempty(parms.fnames_func_ap) && isempty(parms.fnames_func_pa) && ...
              isempty(parms.fnames_dwi_ap)  && isempty(parms.fnames_dwi_pa);
  if all_empty
    errcode = 1;
    errmsg  = 'no processed fmap EPI data found for requested acq_types';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  % processed T1w, for the registration check
  [parms.fnames_t1w,parms.fnames_info_t1w] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'anat','T1w');

  % FreeSurfer recon is optional: used for the brain cutoff check and the
  %   CC mask and boundary cost of the registration check
  parms.fspath = [];
  if ~isempty(parms.fsurf_rootdir)
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
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;
  output_files = {};
  acq_types_required = mmil_getfield(parms,'acq_types_present',{});
  for k=1:length(acq_types_required)
    output_files{end+1} = sprintf('%s/%s_acq-%s_desc-session_postqc.tsv',parms.outdir,parms.outstem,acq_types_required{k});
    output_files{end+1} = sprintf('%s/%s_acq-%s_desc-epi_postqc.tsv',parms.outdir,parms.outstem,acq_types_required{k});
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

function [series_qcstruct,session_rows,errcode,errmsg] = qc_fmap(parms)
  series_qcstruct = {}; session_rows = struct(); errcode = 0; errmsg = [];
  errcodes = []; errmsgs = {};
  fname_means = struct();
  fname_masks = struct();
  group_rows  = struct();
  regT1_qcinfos = struct();

  group_defs = {'func','AP'; 'func','PA'; 'dwi','AP'; 'dwi','PA'};

  for g=1:size(group_defs,1)
    acq  = group_defs{g,1};
    dire = group_defs{g,2};
    if ~ismember(acq,parms.acq_types), continue; end
    label       = sprintf('%s_%s',acq,dire);
    fnames_fld  = sprintf('fnames_%s_%s',acq,lower(dire));
    info_fld    = sprintf('fnames_info_%s_%s',acq,lower(dire));
    fnames      = parms.(fnames_fld);
    fnames_info = parms.(info_fld);
    if isempty(fnames), continue; end

    [stdev_row,scan_rms_vals,scan_corr_vals,scan_com_dists,fname_mean,fname_mask,ec,em] = ...
      calc_fmap_stdev(fnames,fnames_info,label,parms);
    errcodes(end+1) = ec;
    errmsgs{end+1}  = em;
    if ec, continue; end

    fname_means.(label) = fname_mean;
    fname_masks.(label) = fname_mask;

    % interscan aggregates are session-level, one set per direction
    if ~isempty(stdev_row)
      rms_vals = scan_rms_vals(~isnan(scan_rms_vals));
      if ~isempty(rms_vals)
        stdev_row.mean_diff_rms = mean(rms_vals);
      else
        stdev_row.mean_diff_rms = NaN;
      end
      % the least correlated and the farthest field map of the direction
      stdev_row.min_corr = NaN; stdev_row.max_com_dist = NaN;
      if any(~isnan(scan_corr_vals)), stdev_row.min_corr = min(scan_corr_vals); end
      if any(~isnan(scan_com_dists)), stdev_row.max_com_dist = max(scan_com_dists); end
      group_rows.(label) = convert_fmap_group_row(stdev_row,acq,dire);
    end

    % spatial smoothness of each processed field map, measured with no
    % detrending as the raw field maps are: these are one or a few frames, so
    % there is nothing to detrend
    fwhm_vals = calc_fmap_fwhm(fnames,label,fname_mask,parms);

    % one series row per field map scan, with its registration to T1w
    if ~isfield(regT1_qcinfos,acq), regT1_qcinfos.(acq) = {}; end
    for f=1:length(fnames)
      row = convert_fmap_series_row(fnames{f},fnames_info{f},...
        acq,dire,scan_rms_vals(f),scan_corr_vals(f),scan_com_dists(f),fwhm_vals(f,:),...
        parms.participant_id,parms.session_id);
      regT1_qcinfo = check_fmap_regT1(fnames{f},fnames_info{f},parms);
      series_qcstruct{end+1} = add_fmap_regT1_series_cols(row,regT1_qcinfo);
      regT1_qcinfos.(acq){end+1} = regT1_qcinfo;
    end

    % session-level value is the mean across scans of the per-scan values
    if isfield(group_rows,label)
      modality = get_fmap_modality(acq);
      pfx = sprintf('mr_y_qc__post__aut__%s__fm__%s',modality,lower(dire));
      axes_names = {'x','y','z'};
      for k=1:3
        group_rows.(label).(sprintf('%s__%sdir_fwhm',pfx,axes_names{k})) = ...
          mean(fwhm_vals(:,k),'omitnan');
      end
    end
  end

  for k=1:length(parms.acq_types)
    at       = parms.acq_types{k};
    label_ap = sprintf('%s_AP',at);
    label_pa = sprintf('%s_PA',at);
    if ~isfield(group_rows,label_ap) && ~isfield(group_rows,label_pa), continue; end
    modality = get_fmap_modality(at);

    % difference between the AP and PA interscan means
    rms_diff = NaN;
    if isfield(fname_means,label_ap) && isfield(fname_means,label_pa)
      rms_diff = calc_fmap_AP_PA_diff(...
        fname_means.(label_ap),fname_means.(label_pa),at,parms);
    end

    % RMS of the B0 unwarp displacement fields, within the field map brain
    %   mask; prefer the PA mask, since PA is the forward direction and the
    %   displacement fields are exported in the forward reference frame
    fname_mask = [];
    if isfield(fname_masks,label_pa) && ~isempty(fname_masks.(label_pa))
      fname_mask = fname_masks.(label_pa);
    elseif isfield(fname_masks,label_ap) && ~isempty(fname_masks.(label_ap))
      fname_mask = fname_masks.(label_ap);
    end
    [b0_rms,b0_rms_vals,fnames_B0dx] = cmig_calc_b0dx_rms(parms.proc_rootdir,...
      parms.participant_id,parms.session_id,at,parms.tmpdir,...
      'fname_mask',fname_mask,'verbose',parms.verbose,'forceflag',parms.forceflag);

    % brain cutoff for the field map scans, both
    %   phase-encode directions together
    [cutoff_qcinfo,errcode,errmsg] = calc_fmap_cutoff(at,parms);
    if errcode, return; end

    session_rows.(at) = build_fmap_session_row(group_rows,at,rms_diff,b0_rms,...
      parms.participant_id,parms.session_id);
    session_rows.(at) = add_fmap_cutoff_cols(session_rows.(at),cutoff_qcinfo,modality);
    if isfield(regT1_qcinfos,at)
      session_rows.(at) = add_fmap_regT1_cols(session_rows.(at),regT1_qcinfos.(at),modality);
    end

    % b0_rms is a property of the pair, so both scans of a pair carry it
    pfx = 'mr_y_qc__post__series';
    for j=1:length(series_qcstruct)
      if ~strcmp(series_qcstruct{j}.acq,at), continue; end
      series_qcstruct{j}.(sprintf('%s__b0_rms',pfx)) = ...
        select_b0dx_val(series_qcstruct{j}.dir,series_qcstruct{j}.run,...
                        b0_rms_vals,fnames_B0dx);
    end
  end

  nonempty_msgs = errmsgs(~cellfun(@isempty,errmsgs));
  if ~isempty(nonempty_msgs)
    errmsg = strtrim(sprintf('%s ',nonempty_msgs{:}));
  end
  if ~isempty(errcodes) && all(errcodes), errcode = 1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stdev_row,scan_rms_vals,scan_corr_vals,scan_com_dists,fname_mean,fname_mask,errcode,errmsg] = calc_fmap_stdev(fnames,fnames_info,label,parms)
  stdev_row = []; scan_rms_vals = nan(length(fnames),1);
  scan_corr_vals = nan(length(fnames),1); scan_com_dists = nan(length(fnames),1);
  fname_mean = []; fname_mask = []; errcode = 0; errmsg = [];

  stdev_tmpdir = sprintf('%s/fmap_stdev_%s',parms.tmpdir,label);
  stdev_outdir = sprintf('%s/results',stdev_tmpdir);
  mmil_mkdir(stdev_tmpdir);

  fnamelist_mgz = {};
  for f=1:length(fnames)
    fname_mgz = cmig_import_to_mgz(fnames{f},stdev_tmpdir,[],parms.forceflag);
    if isempty(fname_mgz)
      errmsg = sprintf('failed to import fmap file: %s',fnames{f});
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
    fnamelist_mgz{end+1} = fname_mgz;
  end

  if parms.verbose
    fprintf('%s: calculating fmap stdev for %s...\n',mfilename,label);
    tic;
  end

  p = [];
  p.outdir        = stdev_outdir;
  p.outstem       = 'fmap';
  p.plotflag      = 1;
  p.verbose       = parms.verbose;
  p.forceflag     = parms.forceflag;
  a = mmil_parms2args(p);
  results = mmil_interscan_stdev(fnamelist_mgz,a{:});
  if parms.verbose, toc; end

  fname_mean = sprintf('%s/fmap_interscan_first_mean.mgz',stdev_outdir);
  if ~exist(fname_mean,'file'), fname_mean = []; end

  [mask,fname_mask] = get_fmap_brainmask(fname_mean,parms);
  scan_rms_vals = calc_scan_diff_rms_fmap(fnamelist_mgz,stdev_outdir,mask);
  [scan_corr_vals,scan_com_dists] = calc_scan_alignment_fmap(fnamelist_mgz,stdev_outdir,mask);
  copy_fmap_jpgs(fnamelist_mgz,stdev_outdir,label,parms);

  if isempty(results), return; end

  stdev_row.VisitID = sprintf('%s_%s',parms.participant_id,parms.session_id);
  stdev_row.stype   = 'fmap_stdev';
  stdev_row.label   = label;
  stdev_row.nscans  = length(fnamelist_mgz);
  stdev_row.fstem   = sprintf('%s_fmap_stdev_%s',parms.outstem,label);
  inter_mean  = mmil_getfield(results,'mean_mean',NaN);
  inter_stdev = mmil_getfield(results,'stdev_mean',NaN);
  stdev_row.inter_mean  = inter_mean;
  stdev_row.inter_stdev = inter_stdev;
  if ~isnan(inter_mean) && ~isnan(inter_stdev) && inter_stdev > 0
    stdev_row.inter_tSNR = inter_mean / inter_stdev;
    stdev_row.inter_cov  = inter_stdev / inter_mean;
  else
    stdev_row.inter_tSNR = NaN;
    stdev_row.inter_cov  = NaN;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fwhm_vals = calc_fmap_fwhm(fnames,label,fname_mask,parms)
  % FWHM per processed field map scan, no detrending, with the SynthStrip mask
  % built from this direction's interscan mean. outstem keeps the two scans of
  % a pair from colliding when they share the output directory
  fwhm_vals = nan(length(fnames),3);
  if isempty(fname_mask) || ~exist(fname_mask,'file')
    fprintf('%s: WARNING: no brain mask for %s, skipping FWHM\n',mfilename,label);
    return;
  end
  % a folder per direction, since cmig_3dFWHM keeps an existing converted mask
  tmpdir = sprintf('%s/fmap_fwhm_%s',parms.tmpdir,label);
  mmil_mkdir(tmpdir);
  for f=1:length(fnames)
    [~,fstem] = cmig_get_fileparts(fnames{f});
    try
      fwhm_vals(f,:) = cmig_3dFWHM(fnames{f},...
        'fname_mask',fname_mask,'outdir',tmpdir,'outstem',fstem,...
        'detrend_flag',0,'verbose',0,'forceflag',parms.forceflag);
    catch me
      fprintf('%s: WARNING: FWHM failed for %s: %s\n',mfilename,fstem,me.message);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mask,fname_mask] = get_fmap_brainmask(fname_mean,parms)
  mask = []; fname_mask = [];
  if isempty(fname_mean) || ~exist(fname_mean,'file'), return; end
  [fdir,fstem,fext] = fileparts(fname_mean);

  % SynthStrip brain mask, created directly on the fmap reference volume
  fname_tmp = sprintf('%s/%s_mask%s',fdir,fstem,fext);
  try
    cmig_create_brainmask(fname_mean,...
      'fname_mask',fname_tmp,'outdir',fdir,'binarize_flag',1,...
      'warning_flag',1,'verbose',parms.verbose,'forceflag',parms.forceflag);
  catch me
    fprintf('%s: WARNING: brain mask creation failed for %s: %s\n',mfilename,fname_mean,me.message);
    return;
  end

  if exist(fname_tmp,'file')
    fname_mask = fname_tmp;
    vol_mask = fs_load_mgh(fname_mask);
    mask = (vol_mask > 0);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scan_rms_vals = calc_scan_diff_rms_fmap(fnamelist_mgz,stdev_outdir,mask)
  scan_rms_vals = nan(length(fnamelist_mgz),1);
  fname_inter_mean = sprintf('%s/fmap_interscan_first_mean.mgz',stdev_outdir);
  if ~exist(fname_inter_mean,'file'), return; end
  vol_inter_mean = fs_load_mgh(fname_inter_mean);
  if isempty(mask), mask = true(size(vol_inter_mean)); end
  for f=1:length(fnamelist_mgz)
    [~,fstem] = fileparts(fnamelist_mgz{f});
    fname_scan_mean = sprintf('%s/%s_first.mgz',stdev_outdir,fstem);
    if ~exist(fname_scan_mean,'file'), continue; end
    vol_scan_mean = fs_load_mgh(fname_scan_mean);
    diff_vol = vol_scan_mean - vol_inter_mean;
    scan_rms_vals(f) = sqrt(mean(diff_vol(mask).^2));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [scan_corr_vals,scan_com_dists] = calc_scan_alignment_fmap(fnamelist_mgz,stdev_outdir,mask)
  % between-scan alignment of each field map's first frame, within the
  %   direction's brain mask, against the other field maps of the direction
  %   (leave one out): corr is the Pearson r with the voxelwise median of the
  %   others; com_dist is the distance (mm) from the intensity-weighted centre
  %   of mass (scanner coordinates, via the vox2ras) to the median of the
  %   others' centres of mass; NaN with fewer than 2 field maps; with 2, both
  %   get the same
  nscans = length(fnamelist_mgz);
  scan_corr_vals = nan(nscans,1); scan_com_dists = nan(nscans,1);
  if nscans < 2, return; end
  fnames_first = cell(nscans,1);
  for f=1:nscans
    [~,fstem] = fileparts(fnamelist_mgz{f});
    fnames_first{f} = sprintf('%s/%s_first.mgz',stdev_outdir,fstem);
  end
  [scan_corr_vals,scan_com_dists] = calc_image_alignment(fnames_first,mask);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [corr_vals,com_dists] = calc_image_alignment(fnames,mask)
  % leave-one-out correlation and centre of mass distance of a set of
  %   images; a missing image is NaN and is left out of the others' references
  n = length(fnames);
  corr_vals = nan(n,1); com_dists = nan(n,1);
  vals = []; coms = nan(n,3); found = false(n,1); xyz = [];
  for f=1:n
    if ~exist(fnames{f},'file'), continue; end
    [vol,M] = fs_load_mgh(fnames{f});
    if isempty(mask), mask = true(size(vol)); end
    if isempty(xyz)
      ind = find(mask);
      [i,j,k] = ind2sub(size(mask),ind);
      xyz = (M(1:3,:) * [i j k ones(length(ind),1)]')';
      vals = nan(length(ind),n);
    end
    vals(:,f) = vol(ind);
    % negative intensities would pull the weighted centre outside the image
    w = max(vals(:,f),0);
    if sum(w) > 0
      coms(f,:) = sum(bsxfun(@times,xyz,w),1) / sum(w);
    end
    found(f) = true;
  end
  for f=1:n
    others = find(found & (1:n)' ~= f);
    if ~found(f) || isempty(others), continue; end
    ref = median(vals(:,others),2);
    x = vals(:,f) - mean(vals(:,f));
    y = ref - mean(ref);
    corr_vals(f) = sum(x.*y) / sqrt(sum(x.^2) * sum(y.^2));
    com_dists(f) = norm(coms(f,:) - median(coms(others,:),1));
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rms_diff = calc_fmap_AP_PA_diff(fname_mean_ap,fname_mean_pa,acq_type,parms)
  rms_diff = NaN;
  if isempty(fname_mean_ap) || ~exist(fname_mean_ap,'file'), return; end
  if isempty(fname_mean_pa) || ~exist(fname_mean_pa,'file'), return; end

  if parms.verbose
    fprintf('%s: computing AP-PA difference for %s fmaps...\n',mfilename,acq_type);
  end

  try
    [vol_ap,M] = fs_load_mgh(fname_mean_ap);
    vol_pa     = fs_load_mgh(fname_mean_pa);
    diff_vol   = vol_ap - vol_pa;

    fname_diff_mgz = sprintf('%s/%s_desc-APvsPA_diff.mgz',parms.tmpdir,parms.outstem);
    fs_save_mgh(diff_vol,fname_diff_mgz,M);
    fname_jpg = sprintf('%s/%s_acq-%s_desc-APvsPA_diff.jpg',parms.outdir,parms.outstem,acq_type);
    if ~exist(fname_jpg,'file') || parms.forceflag
      outstem = sprintf('%s_acq-%s_desc-APvsPA_diff',parms.outstem,acq_type);
      mmil_subplots(fname_diff_mgz,'outdir',parms.outdir,'outstem',outstem,...
        'forceflag',parms.forceflag);
    end

    mask_ap = get_fmap_brainmask(fname_mean_ap,parms);
    mask_pa = get_fmap_brainmask(fname_mean_pa,parms);
    if ~isempty(mask_ap) && ~isempty(mask_pa)
      mask = mask_ap & mask_pa;
    elseif ~isempty(mask_ap)
      mask = mask_ap;
    elseif ~isempty(mask_pa)
      mask = mask_pa;
    else
      mask = true(size(diff_vol));
    end

    rms_diff = sqrt(mean(diff_vol(mask).^2));
  catch me
    fprintf('%s: WARNING: AP-PA diff failed for %s: %s\n',mfilename,acq_type,me.message);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = convert_fmap_group_row(qcinfo,acq,dire)
  % interscan aggregates for one acq type and phase-encode direction
  row = struct();
  pfx = sprintf('mr_y_qc__post__aut__%s__fm__%s',get_fmap_modality(acq),lower(dire));
  row.(sprintf('%s__nscans_count',pfx))          = qcinfo.nscans;
  % no interscan columns for dwi, which has one field map per direction
  if strcmpi(acq,'dwi'), return; end
  row.(sprintf('%s__inter__signal_mean',pfx))    = qcinfo.inter_mean;
  row.(sprintf('%s__inter__signal_std',pfx))     = qcinfo.inter_stdev;
  row.(sprintf('%s__inter__temp__snr_mean',pfx)) = qcinfo.inter_tSNR;
  row.(sprintf('%s__inter__signal_cov',pfx))     = qcinfo.inter_cov;
  row.(sprintf('%s__inter__diff_rms',pfx))       = mmil_getfield(qcinfo,'mean_diff_rms',NaN);
  row.(sprintf('%s__inter__corr_min',pfx))       = mmil_getfield(qcinfo,'min_corr',NaN);
  row.(sprintf('%s__inter__com__dist_max',pfx))  = mmil_getfield(qcinfo,'max_com_dist',NaN);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = convert_fmap_series_row(fname,fname_info,acq,dire,diff_rms,corr_val,com_dist,fwhm,pid,sid)
  % one row per field map scan; the dir entity column identifies the
  %   phase-encode direction, so it is not repeated in the column names;
  %   series columns are mr_y_qc__post__series__{measure}, with no modality
  row.participant_id = pid;
  row.session_id     = sid;
  [~,fstem] = cmig_get_fileparts(fname);
  row.modality = 'epi';
  row.acq_time = get_fmap_acq_time(fname,fname_info);
  row.filename = sprintf('fmap/%s.nii.gz',rm_desc_preproc(fstem));
  tok = regexp(fstem,'(?:^|_)run-([^_]+)','tokens','once');
  if ~isempty(tok), row.run = tok{1}; else row.run = 'n/a'; end
  row.acq = lower(acq);
  row.dir = lower(dire);
  if ~strcmpi(acq,'dwi')
    row.mr_y_qc__post__series__inter__diff_rms = diff_rms;
    row.mr_y_qc__post__series__inter__corr_score = corr_val;
    row.mr_y_qc__post__series__inter__com_dist = com_dist;
  end
  pfx = 'mr_y_qc__post__series';
  axes_names = {'x','y','z'};
  for k=1:3
    if numel(fwhm)>=k, v = fwhm(k); else v = NaN; end
    row.(sprintf('%s__%sdir_fwhm',pfx,axes_names{k})) = v;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qcinfo = check_fmap_regT1(fname,fname_info,parms)
  % registration of one field map to T1w, as MMIL_auto_postQC_Exam checked
  %   the registration reference field map; the field maps are resampled to
  %   the reference's grid and share its RegistrationMatrixT1w
  %   failures are warnings, leaving the regT1 values NaN
  qcinfo = cmig_init_epi_qcinfo(fname,'fmap',parms.participant_id,parms.session_id);
  if isempty(parms.fnames_t1w)
    fprintf('%s: WARNING: no processed T1w, skipping the registration check of %s\n',...
      mfilename,fname);
    return;
  end
  [qcinfo,ec,em] = cmig_read_proc_json_epi(qcinfo,fname,fname_info);
  if ec || isempty(qcinfo.registration_matrix_T1)
    fprintf('%s: WARNING: no registration matrix for %s, skipping its registration check %s\n',...
      mfilename,fname,em);
    return;
  end
  [~,fstem] = cmig_get_fileparts(fname);
  entities = regexprep(rm_desc_preproc(fstem),'_[^_]+$','');
  fname_plot = sprintf('%s/%s_regT1.jpg',parms.outdir,entities);
  [qcinfo,ec,em] = cmig_check_regT1_epi(qcinfo,fname,...
    'fnames_t1w',parms.fnames_t1w,'tmpdir',parms.tmpdir,'outdir',parms.outdir,...
    'fname_plot',fname_plot,...
    'fspath',parms.fspath,'verbose',parms.verbose,'forceflag',parms.forceflag);
  if ec
    fprintf('%s: WARNING: registration check of %s failed: %s\n',mfilename,fname,em);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = add_fmap_regT1_series_cols(row,qcinfo)
  % per field map regT1 values, named as postqc func and dwi name theirs
  pfx = 'mr_y_qc__post__series';
  row.(sprintf('%s__regt1__imgcorr_score',pfx))  = qcinfo.regT1_corr;
  row.(sprintf('%s__rigid_score',pfx))           = qcinfo.regT1_norm;
  row.(sprintf('%s__warp_score',pfx))            = qcinfo.regT1_warp;
  row.(sprintf('%s__regt1__imgcorr1_score',pfx)) = qcinfo.regT1_imgcorr1;
  row.(sprintf('%s__regt1__imgcorr2_score',pfx)) = qcinfo.regT1_imgcorr2;
  row.(sprintf('%s__regt1__imgdist1_score',pfx)) = qcinfo.regT1_imgdist3;
  row.(sprintf('%s__regt1__imgdist2_score',pfx)) = qcinfo.regT1_imgdist4;
  row.(sprintf('%s__regt1__segreg__cost_score',pfx)) = qcinfo.regT1_cost0;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = add_fmap_regT1_cols(session_row,qcinfos,modality)
  % session regT1 values: the median over the field maps of the acq type,
  %   both phase-encode directions, as postqc func takes over its runs
  src  = {'regT1_corr','regT1_norm','regT1_warp',...
          'regT1_imgcorr1','regT1_imgcorr2','regT1_imgdist3','regT1_imgdist4',...
          'regT1_cost0'};
  dest = {'regt1__imgcorr_score','rigid_score','warp_score',...
          'regt1__imgcorr1_score','regt1__imgcorr2_score',...
          'regt1__imgdist1_score','regt1__imgdist2_score',...
          'regt1__segreg__cost_score'};
  % the fMRI field maps' check is the fMRI registration check, so its columns
  %   take the fmri names (rigid_score and warp_score as released); the dMRI
  %   scores are of the dMRI scans, so the dMRI field maps' keep __fm
  if strcmp(modality,'fmri')
    pfx = sprintf('mr_y_qc__post__aut__%s',modality);
  else
    pfx = sprintf('mr_y_qc__post__aut__%s__fm',modality);
  end
  for k=1:length(src)
    vals = cellfun(@(q) mmil_getfield(q,src{k},NaN), qcinfos);
    vals = vals(~isnan(vals));
    if isempty(vals)
      session_row.(sprintf('%s__%s',pfx,dest{k})) = NaN;
    else
      session_row.(sprintf('%s__%s',pfx,dest{k})) = median(vals);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modality = get_fmap_modality(acq)
  % the data dictionary has no fmap token; field map variables are named
  %   for the modality whose distortion they describe
  switch lower(acq)
    case 'func'
      modality = 'fmri';
    case 'dwi'
      modality = 'dmri';
    otherwise
      modality = lower(acq);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function acq_time = get_fmap_acq_time(fname,fname_info)
  acq_time = '';
  if isempty(fname_info) || ~exist(fname_info,'file')
    fname_info = regexprep(fname,'\.nii.*','.json');
  end
  if ~exist(fname_info,'file'), return; end
  try
    info = jsondecode(fileread(fname_info));
    acq_time = mmil_getfield(info,'AcquisitionDateTime','');
  catch
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = select_b0dx_val(dir_label,run_label,b0_rms_vals,fnames_B0dx)
  % b0_rms is a property of the forward/reverse pair; the B0dx file name
  %   carries both scans' run labels in its desc entity, in fixed AP then PA
  %   slots, so match the scan's own direction against the matching slot
  val = NaN;
  if isempty(b0_rms_vals), return; end
  if length(b0_rms_vals) == 1
    val = b0_rms_vals(1);
    return;
  end
  switch upper(dir_label)
    case 'AP', j = 1;
    case 'PA', j = 2;
    otherwise, return;
  end
  % a scan without a run entity has run label 'n/a' and an empty desc slot
  if strcmp(run_label,'n/a'), run_label = ''; end
  for k=1:min(length(fnames_B0dx),length(b0_rms_vals))
    [~,fstem] = cmig_get_fileparts(fnames_B0dx{k});
    tok = regexp(fstem,'_desc-AP([^_]*)PA([^_]*)_B0dx$','tokens','once');
    if isempty(tok), continue; end
    if strcmp(run_label,tok{j})
      val = b0_rms_vals(k);
      return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function copy_fmap_jpgs(fnamelist_mgz,stdev_outdir,label,parms)
  parts = strsplit(label,'_');
  acq = parts{1};
  dir = parts{2};
  session_map = {...
    'fmap_interscan_first_mean',      'mean';...
    'fmap_interscan_first_mean_mask', 'mask';...
    'fmap_interscan_first_stdev',     'std';...
  };
  for j=1:size(session_map,1)
    src = sprintf('%s/%s.jpg',stdev_outdir,session_map{j,1});
    dst = sprintf('%s/%s_acq-%s_dir-%s_desc-inter_%s.jpg',...
      parms.outdir,parms.outstem,acq,dir,session_map{j,2});
    if exist(src,'file') && (~exist(dst,'file') || parms.forceflag)
      copyfile(src,dst);
    end
  end
  for f=1:length(fnamelist_mgz)
    [~,fstem] = fileparts(fnamelist_mgz{f});
    toks     = regexp(rm_desc_preproc(fstem),'(?:sub|ses|task|acq|rec|dir|run|hemi|space|res|den|label|desc)-[^_]+','match');
    entities = strjoin(toks,'_');
    src = sprintf('%s/%s_first.jpg',stdev_outdir,fstem);
    dst = sprintf('%s/%s_first.jpg',parms.outdir,entities);
    if exist(src,'file') && (~exist(dst,'file') || parms.forceflag)
      copyfile(src,dst);
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = build_fmap_session_row(group_rows,acq,diff_rms,b0_rms,pid,sid)
  % one session row per acq type, flattening both phase-encode directions
  session_row = struct();
  session_row.participant_id = pid;
  session_row.session_id     = sid;
  modality = get_fmap_modality(acq);
  labels = {sprintf('%s_AP',acq),sprintf('%s_PA',acq)};
  for k=1:length(labels)
    if ~isfield(group_rows,labels{k}), continue; end
    row  = group_rows.(labels{k});
    flds = fieldnames(row);
    for j=1:length(flds)
      session_row.(flds{j}) = row.(flds{j});
    end
  end
  session_row.(sprintf('mr_y_qc__post__aut__%s__fm__ap__pa__diff_rms',modality)) = diff_rms;
  session_row.(sprintf('mr_y_qc__post__aut__%s__b0_rms',modality)) = b0_rms;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qcinfo,errcode,errmsg] = calc_fmap_cutoff(acq,parms)
  % brain cutoff across the field map scans of one acq type,
  %   both phase-encode directions together, as postqc func and dwi do for
  %   their scans
  qcinfo = []; errcode = 0; errmsg = [];
  fnames      = cat(2,parms.(sprintf('fnames_%s_ap',acq)),...
                      parms.(sprintf('fnames_%s_pa',acq)));
  fnames_info = cat(2,parms.(sprintf('fnames_info_%s_ap',acq)),...
                      parms.(sprintf('fnames_info_%s_pa',acq)));
  if isempty(fnames), return; end

  qcinfo.brain_cutoff           = NaN;
  qcinfo.brain_cutoff_dorsal    = NaN;
  qcinfo.brain_cutoff_ventral   = NaN;

  fnamestem = sprintf('fmap_%s',acq);
  [qcinfo,errcode,errmsg] = cmig_detect_brain_cutoff(qcinfo,...
    fnames,fnames_info,fnamestem,...
    'tmpdir',parms.tmpdir,'outdir',parms.outdir,'outstem',parms.outstem,...
    'fname_plot',sprintf('%s/%s_acq-%s_braincutoff.jpg',parms.outdir,parms.outstem,acq),...
    'fspath',parms.fspath,'verbose',parms.verbose,'forceflag',parms.forceflag);
  if errcode, return; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function session_row = add_fmap_cutoff_cols(session_row,qcinfo,modality)
  % append session-level cutoff columns for one acq type
  if isempty(qcinfo), return; end
  pfx = sprintf('mr_y_qc__post__aut__%s__fm',modality);
  session_row.(sprintf('%s__brain__cutoff_score',pfx))         = qcinfo.brain_cutoff;
  session_row.(sprintf('%s__brain__cutoff__dorsal_score',pfx)) = qcinfo.brain_cutoff_dorsal;
  session_row.(sprintf('%s__brain__cutoff__ventral_score',pfx)) = qcinfo.brain_cutoff_ventral;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post = pe_toward_posterior(info)
  % true when the phase encoding runs toward posterior (AP): the
  %   PhaseEncodingDirection axis and sign through the image's vox2ras
  %   (Affine), so right in any output orientation (std_orient of
  %   cmig_proc_bids_func); without an Affine, or for an axis with no A/P
  %   component, the sign as in dcm2niix's LAS, where j points anterior
  post = [];
  pe = mmil_getfield(info,'PhaseEncodingDirection','');
  if isempty(pe) || ~ischar(pe), return; end
  a = strfind('ijk',pe(1));
  if isempty(a), return; end
  s = 1;
  if any(pe=='-'), s = -1; end
  M = mmil_getfield(info,'Affine',[]);
  if isnumeric(M) && isequal(size(M),[4 4]) && ...
      abs(M(2,a)) > 0.5*norm(M(1:3,a))
    post = s*M(2,a) < 0;
  else
    post = s < 0;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function direction = get_pe_direction(fname,fname_info)
  direction = [];
  [~,fstem] = cmig_get_fileparts(fname);
  tok = regexp(fstem,'(?:^|_)dir-([^_]+)','tokens','once');
  if ~isempty(tok)
    val = upper(tok{1});
    if ismember(val,{'AP','PA'})
      direction = val;
      return;
    end
  end
  if isempty(fname_info) || ~exist(fname_info,'file')
    fname_info = regexprep(fname,'\.nii.*','.json');
  end
  if ~exist(fname_info,'file'), return; end
  try
    info = jsondecode(fileread(fname_info));
    pe   = mmil_getfield(info,'PhaseEncodingDirection','');
    if isempty(pe) || ~strcmp(pe(1),'j'), return; end
    if pe_toward_posterior(info)
      direction = 'AP';
    else
      direction = 'PA';
    end
  catch
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fstem = rm_desc_preproc(fstem)
  % the rawdata name of a proc file, which adds desc-preproc: the filename
  %   column joins with raw QC on it, and the jpg names carry the run entities
  fstem = regexprep(fstem,'_desc-preproc(?=_[^_]+$)','');
return;
