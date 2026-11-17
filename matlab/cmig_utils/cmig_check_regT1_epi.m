function [auto_qcinfo,errcode,errmsg] = cmig_check_regT1_epi(auto_qcinfo, fname_epi, varargin)
%function [auto_qcinfo,errcode,errmsg] = cmig_check_regT1_epi(auto_qcinfo, fname_epi, [options])
%
% Purpose: check registration of an EPI scan to T1w
%   resamples EPI first frame to T1 space and computes registration metrics
%
% Required Input:
%   auto_qcinfo: EPI QC info struct (must contain registration_matrix_T1 field)
%   fname_epi: EPI file name (NIfTI)
%
% Optional Input ('key',value):
%   'fnames_t1w': cell array of processed T1w file names
%     {default = {}}
%   'tmpdir': temporary/working directory for intermediate files
%     {default = ''}
%   'outdir': output directory for plots
%     {default = ''}
%   'fspath': FreeSurfer subject directory (used to create CC mask and
%     for the mri_segreg boundary cost, regT1_cost0)
%     {default = ''}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   auto_qcinfo: updated QC info struct with regT1 metrics
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 06/08/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Last Mod: 10/05/2026 by Don Hagler
%

errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,2), errcode = 1; return; end
parms = mmil_args2parms(varargin,{...
  'fnames_t1w',{},[],...
  'tmpdir','',[],...
  'outdir','',[],...
  'fname_plot','',[],...
  'fspath','',[],...
  'verbose',false,[],...
  'forceflag',false,[],...
});

checkname = sprintf('%s registration to T1',upper(auto_qcinfo.stype));
[~,fstem] = cmig_get_fileparts(fname_epi);

fname_t1w = [];
for f=1:length(parms.fnames_t1w)
  if exist(parms.fnames_t1w{f},'file')
    fname_t1w = parms.fnames_t1w{f};
    break;
  end
end
if isempty(fname_t1w)
  errmsg = sprintf('no proc T1w found, skipping %s',checkname);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  errcode = 1;
  return;
end

% registration_matrix_T1 = M_T1_to_EPI stored in proc JSON
M_T1_to_EPI = auto_qcinfo.registration_matrix_T1;
if isempty(M_T1_to_EPI)
  errmsg = sprintf('registration_matrix_T1 not found in proc JSON, skipping %s',checkname);
  fprintf('%s: ERROR: %s\n',mfilename,errmsg);
  errcode = 1;
  return;
end

fname_out = sprintf('%s/%s_regT1.mat',parms.tmpdir,fstem);
if ~exist(fname_out,'file') || parms.forceflag
  % convert T1w to mgz
  fname_t1w_mgz = cmig_import_to_mgz(fname_t1w,parms.tmpdir);
  if isempty(fname_t1w_mgz)
    errmsg = sprintf('T1w mgz missing, skipping %s',checkname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end

  % convert EPI to mgz
  fname_epi_mgz = cmig_import_to_mgz(fname_epi,parms.tmpdir);
  if isempty(fname_epi_mgz)
    errmsg = sprintf('failed to import EPI, skipping %s',checkname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end

  % extract first frame from EPI
  fname_epi_f0 = sprintf('%s/%s_f0.mgz',parms.tmpdir,fstem);
  if ~exist(fname_epi_f0,'file') || parms.forceflag
    try
      [vol_epi,M_epi] = fs_load_mgh(fname_epi_mgz,[],1);
      fs_save_mgh(vol_epi,fname_epi_f0,M_epi);
    catch me
      errmsg = sprintf('failed to extract EPI first frame: %s',me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
  end

  % resample EPI f0 to T1 space using M_T1_to_EPI
  % M_reg for mmil_resample_vol maps FROM ref (T1) TO input (EPI) = M_T1_to_EPI
  fname_epi_resT1 = sprintf('%s/%s_f0_resT1.mgz',parms.tmpdir,fstem);
  if ~exist(fname_epi_resT1,'file') || parms.forceflag
    try
      [vol_epi_f0,M_epi_f0] = fs_load_mgh(fname_epi_f0);
      [M_T1,nvox_T1] = mmil_load_mgh_info(fname_t1w_mgz,parms.forceflag,parms.tmpdir);
      [vol_epi_resT1,M_epi_resT1] = mmil_resample_vol(vol_epi_f0,M_epi_f0,...
        'M_ref',M_T1,'nvox_ref',nvox_T1(1:3),...
        'M_reg',M_T1_to_EPI,'interpm',1,'bclamp',0);
      fs_save_mgh(vol_epi_resT1,fname_epi_resT1,M_epi_resT1);
    catch me
      errmsg = sprintf('failed to resample EPI to T1 space: %s',me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
  end
  if ~exist(fname_epi_resT1,'file')
    errmsg = sprintf('EPI resT1 file missing, skipping %s',checkname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end

  % mri_segreg boundary cost of the registration, given a FreeSurfer recon
  regT1_cost0 = NaN;
  if ~isempty(parms.fspath)
    [regT1_cost0,ec,em] = cmig_check_regT1_segreg(fname_epi_resT1,parms.fspath,...
      'contrast','bold','tmpdir',parms.tmpdir,...
      'verbose',parms.verbose,'forceflag',parms.forceflag);
    if ec
      fprintf('%s: WARNING: no boundary cost for %s: %s\n',mfilename,checkname,em);
    end
  end

  % optionally create CC mask
  fname_mask = [];
  if ~isempty(parms.fspath)
    fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
    if exist(fname_aseg,'file')
      fname_cc = cmig_create_aseg_ccmask(parms.fspath,parms.tmpdir,parms.verbose,parms.forceflag);
      if ~isempty(fname_cc) && exist(fname_cc,'file')
        fname_mask = cmig_resample_mask_to_proc(fname_cc,fname_t1w_mgz,parms.tmpdir,parms.verbose,parms.forceflag);
      end
    end
  end

  % without a CC mask, whole-volume correlations are not comparable to the
  %   CC-masked values, so record NaN metrics instead
  if isempty(fname_mask) || ~exist(fname_mask,'file')
    fprintf('%s: WARNING: no CC mask for %s, setting regT1 metrics to NaN\n',...
      mfilename,checkname);
    regT1_corr = NaN; regT1_norm = NaN; regT1_warp = NaN; results = [];
    save(fname_out,'regT1_corr','regT1_norm','regT1_warp','regT1_cost0','results');
  else
    if parms.verbose
      fprintf('%s: checking %s registration to T1...\n',mfilename,upper(auto_qcinfo.stype));
      tic;
    end
    p = [];
    if ~isempty(parms.fname_plot)
      p.fname_plot = parms.fname_plot;
    else
      p.fname_plot = sprintf('%s/%s_regT1.jpg',parms.outdir,fstem);
    end
    p.forceflag = parms.forceflag;
    p.fname_mask = fname_mask;
    a = mmil_parms2args(p);
    try
      [regT1_corr,regT1_norm,regT1_warp,results] = ...
        mmil_check_regT1(fname_epi_resT1,fname_t1w_mgz,a{:});
      if parms.verbose, toc; end
      save(fname_out,'regT1_corr','regT1_norm','regT1_warp','regT1_cost0','results');
    catch me
      errmsg = sprintf('mmil_check_regT1 failed for EPI: %s',me.message);
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
