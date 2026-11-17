function [outdir,errcode,errmsg] = cmig_postqc_bids_fsurf(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_postqc_bids_fsurf(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Automated post-processing QC for BIDS FreeSurfer derivatives
%   Checks topology defects and optionally compares proc T1w to FreeSurfer input
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%     processed data at {proc_rootdir}/{participant_id}/{session_id}/anat
%     output from cmig_proc_bids_anat (used to locate proc T1w for input match check)
%   fsurf_rootdir: BIDS FreeSurfer derivatives root directory
%     FreeSurfer recons at {fsurf_rootdir}/{participant_id}_{session_id}
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
%     form: {output_rootdir}/{participant_id}/{session_id}/fsurf
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 06/10/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/18/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 07/09/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/12/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 10/12/2026 by Don Hagler
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
  [all_qcstruct,errcode,errmsg] = qc_fsurf(parms);
  if ~errcode && ~isempty(all_qcstruct)
    fname_tsv = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);
    cmig_write_tabdat_postqc(all_qcstruct,fname_tsv,parms.forceflag);
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
  ... % hidden
    'vox2ras_maxdiff',1e-3,[1e-16,1],...
  });

  parms.outdir  = sprintf('%s/%s/%s/fsurf',...
    parms.output_rootdir,parms.participant_id,parms.session_id);
  parms.tmpdir  = sprintf('%s/autoQC',parms.outdir);
  parms.outstem = sprintf('%s_%s',parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];
  pid = parms.participant_id;
  sid = parms.session_id;

  fprintf('%s: checking input for %s_%s (modality=fsurf)\n',mfilename,pid,sid);

  [parms.fspath,parms.fs_version,parms.fs_build,parms.fs_complete_flag,errcode,errmsg] = ...
    cmig_check_bids_fsurf(parms.participant_id,parms.session_id,parms.fsurf_rootdir);
  parms.status_extra = struct('fs_version',parms.fs_version,...
    'fs_build',parms.fs_build,'fs_complete_flag',parms.fs_complete_flag);
  if errcode
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end

  parms.fnames_t1w = []; parms.fnames_info_t1w = [];
  [parms.fnames_t1w,parms.fnames_info_t1w] = ...
    cmig_get_bids_data(pid,sid,parms.proc_rootdir,'anat','T1w');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;
  output_files = {};

  output_files{end+1} = sprintf('%s/%s_desc-session_postqc.tsv',parms.outdir,parms.outstem);

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

function [all_qcstruct,errcode,errmsg] = qc_fsurf(parms)
  all_qcstruct = {};
  errcode = 0; errmsg = [];

  if parms.verbose
    fprintf('%s: running FreeSurfer QC for: %s\n',mfilename,parms.fspath);
  end

  fsurf_qcinfo = init_fsurf_qcinfo(parms);

  if ~isempty(parms.fnames_t1w)
    [fsurf_qcinfo,tmp_errcode,tmp_errmsg] = check_fsurf_input(fsurf_qcinfo,parms);
    if tmp_errcode
      fprintf('%s: WARNING: input match check failed for %s_%s: %s\n',...
        mfilename,parms.participant_id,parms.session_id,tmp_errmsg);
    end
  end

  [fsurf_qcinfo,tmp_errcode,tmp_errmsg] = check_fsurf_topofix(fsurf_qcinfo,parms);
  if tmp_errcode
    fprintf('%s: WARNING: topofix check failed for %s_%s: %s\n',...
      mfilename,parms.participant_id,parms.session_id,tmp_errmsg);
  end

  all_qcstruct{1} = convert_fsurf_row(fsurf_qcinfo,parms.participant_id,parms.session_id);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fsurf_qcinfo = init_fsurf_qcinfo(parms)
  fsurf_qcinfo = [];
  fsurf_qcinfo.VisitID                   = sprintf('%s_%s',parms.participant_id,parms.session_id);
  fsurf_qcinfo.stype                     = 'FSURF';
  fsurf_qcinfo.fstem                     = sprintf('%s_fsurf',parms.outstem);
  fsurf_qcinfo.fspath                    = parms.fspath;
  fsurf_qcinfo.fs_version                = parms.fs_version;
  fsurf_qcinfo.fs_build                  = parms.fs_build;
  fsurf_qcinfo.fs_complete_flag          = parms.fs_complete_flag;
  fsurf_qcinfo.input_match               = NaN;
  fsurf_qcinfo.topofix_num_defects       = NaN;
  fsurf_qcinfo.topofix_num_patches       = NaN;
  fsurf_qcinfo.topofix_min_patch_size    = NaN;
  fsurf_qcinfo.topofix_max_patch_size    = NaN;
  fsurf_qcinfo.topofix_mean_patch_size   = NaN;
  fsurf_qcinfo.topofix_median_patch_size = NaN;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fsurf_qcinfo,errcode,errmsg] = check_fsurf_input(fsurf_qcinfo,parms)
  errcode = 0; errmsg = [];
  % NOTE: cmig_fsurf_bids runs recon-all on the proc T1w, and
  %   mri/orig/001.mgz is the unconformed copy of the first input, so with
  %   a single input T1w the two volumes should be identical and
  %   input_match=1 is expected; a mismatch means the recon was built from
  %   different input than the current proc T1w
  fname_fsurf = sprintf('%s/mri/orig/001.mgz',parms.fspath);

  fname_out = sprintf('%s/fsurf_input_match.mat',parms.tmpdir);
  if ~exist(fname_out,'file') || parms.forceflag
    if ~exist(fname_fsurf,'file')
      errmsg = sprintf('%s not found, skipping input match check',fname_fsurf);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
    input_match = 1;
    try
      [vol_proc,M_proc,volsz_proc] = cmig_read_nifti(parms.fnames_t1w{1});
      [vol_fsurf,M_fsurf,~,volsz_fsurf] = fs_load_mgh(fname_fsurf);
      if any(volsz_proc(1:3) ~= volsz_fsurf(1:3))
        if parms.verbose
          fprintf('%s: NOTE: volume dims differ: proc %s vs 001.mgz %s\n',...
            mfilename,num2str(volsz_proc(1:3)),num2str(volsz_fsurf(1:3)));
        end
        input_match = 0;
      end
      if input_match && any(abs(M_proc(:)-M_fsurf(:)) > parms.vox2ras_maxdiff)
        if parms.verbose
          fprintf('%s: NOTE: vox2ras matrices differ for proc T1w vs 001.mgz\n',mfilename);
        end
        input_match = 0;
      end
      if input_match && any(vol_proc(:) ~= vol_fsurf(:))
        if parms.verbose
          fprintf('%s: NOTE: voxel values differ for proc T1w vs 001.mgz\n',mfilename);
        end
        input_match = 0;
      end
    catch me
      errmsg = sprintf('input match check failed: %s',me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
    save(fname_out,'input_match');
  else
    tmp = load(fname_out,'input_match');
    input_match = tmp.input_match;
  end
  fsurf_qcinfo.input_match = input_match;
  if input_match==1 && parms.verbose
    fprintf('%s: proc T1w matches FreeSurfer 001.mgz\n',mfilename);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fsurf_qcinfo,errcode,errmsg] = check_fsurf_topofix(fsurf_qcinfo,parms)
  % defects are the surface holes (fs_euler); patches are the defects
  %   mris_fix_topology corrected, sized in vertices, from ?h.defect_labels
  errcode = 0; errmsg = [];
  [subjdir,subj,ext] = fileparts(parms.fspath);
  subj = [subj ext];

  fname_out = sprintf('%s/fsurf_topofix_defects.mat',parms.tmpdir);
  if ~exist(fname_out,'file') || parms.forceflag
    num_defects = NaN;
    try
      num_defects = fs_euler(subj,subjdir,parms.tmpdir,parms.verbose);
    catch me
      errmsg = sprintf('fs_euler failed: %s',me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
    [patch_sizes,errcode,errmsg] = cmig_topofix_defects(parms.fspath);
    if errcode, return; end
    save(fname_out,'num_defects','patch_sizes');
  else
    tmp = load(fname_out);
    num_defects = tmp.num_defects;
    patch_sizes = tmp.patch_sizes;
  end

  fsurf_qcinfo.topofix_num_defects = num_defects;
  fsurf_qcinfo.topofix_num_patches = length(patch_sizes);
  if ~isempty(patch_sizes)
    fsurf_qcinfo.topofix_min_patch_size    = min(patch_sizes);
    fsurf_qcinfo.topofix_max_patch_size    = max(patch_sizes);
    fsurf_qcinfo.topofix_mean_patch_size   = mean(patch_sizes);
    fsurf_qcinfo.topofix_median_patch_size = median(patch_sizes);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function row = convert_fsurf_row(fsurf_qcinfo, participant_id, session_id)
  row.participant_id = participant_id;
  row.session_id     = session_id;
  row.mr_y_qc__post__aut__smri__fsurf_ver    = fsurf_qcinfo.fs_version;
  row.mr_y_qc__post__aut__smri__fsurf_build  = fsurf_qcinfo.fs_build;
  row.mr_y_qc__post__aut__smri__fsurf__complete_indicator = fsurf_qcinfo.fs_complete_flag;
  row.mr_y_qc__post__aut__smri__input__match_indicator = fsurf_qcinfo.input_match;
  row.mr_y_qc__post__aut__smri__topodfct_count         = fsurf_qcinfo.topofix_num_defects;
  row.mr_y_qc__post__aut__smri__topopatch_count         = fsurf_qcinfo.topofix_num_patches;
  row.mr_y_qc__post__aut__smri__topopatch_min  = fsurf_qcinfo.topofix_min_patch_size;
  row.mr_y_qc__post__aut__smri__topopatch_max  = fsurf_qcinfo.topofix_max_patch_size;
  row.mr_y_qc__post__aut__smri__topopatch_mean = fsurf_qcinfo.topofix_mean_patch_size;
  row.mr_y_qc__post__aut__smri__topopatch_med  = fsurf_qcinfo.topofix_median_patch_size;
return;
