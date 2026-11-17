function epi_topup_estimateB0(fname_for,fname_rev,varargin)
%function epi_topup_estimateB0(fname_for,fname_rev,[options])
%
% Purpse: estimate distortions caused by B0 susceptibility artifacts
%   using FSL's topup
%
% Usage: epi_topup_estimateB0(fname_for,fname_rev,'key','value'...)
%
% Required Parameters:
%   fname_for: input file name with "forward" phase encode polarity
%   fname_rev: input file name with "reverse" phase encode polarity
%
% Optional Parameters:
%   'fname_dx': output file name for displacement field
%     If empty, will construct from fname_for and fname_rev
%     Will use stem of fname_dx as file stem for field coefficients and motion parameters
%     {default = []}
%   'fname_for_out': output file name for forward scan after B0 unwarping
%     If empty, will not save this output
%     {default = []}
%   'fname_rev_out': output file name for reverse scan after B0 unwarping
%     If empty, will not save this output
%     {default = []}
%   'swapxy_flag': [0|1] whether to swap x and y dimensions before estimating B0
%     {default = 0}
%   'npad': number of slices to add to each end of volume
%     before estimating field (extra slices removed after)
%     {default = 5}
%   'readout_time': total readout time (s) of fname_for and fname_rev
%     written to topup's acquisition file; the displacement field is
%     topup's field (Hz) times readout_time, in voxels
%     1 gives the field in Hz, numerically equal to voxels
%     {default = 1}
%   'forceflag': [0|1] whether to overwrite existing fname_dx
%     {default = 0}
%
% Created:  07/03/19 by Don Hagler
% Prev Mod: 07/10/19 by Josh Kuperman
% Prev Mod: 11/02/19 by Don Hagler
% Prev Mod: 03/21/20 by Don Hagler
% Prev Mod: 06/30/20 by Don Hagler
% Prev Mod: 07/01/20 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Prev Mod: 09/25/2026 by Don Hagler
% Last Mod: 10/09/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on epi_estimateB0
% Created:  02/02/07 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
if ~mmil_check_nargs(nargin,2), return; end;
parms = check_input_parms(fname_for,fname_rev,varargin);

% check output files
all_exist = check_output_files(parms);
if all_exist && ~parms.forceflag, return; end;

% check that input volumes share dimensions and orientations
parms = check_inputs_match(parms);

% prepare input files
parms = prep_inputs(parms);

% estimate field
parms = estimate_field(parms);

% convert output filese
convert_output(parms);

% remove temporary files
cleanup_tmpdir(parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input_parms(fname_for,fname_rev,options)
  parms = mmil_args2parms(options,{...
    'fname_for',fname_for,[],...
    'fname_rev',fname_rev,[],...
  ...
    'fname_dx',[],[],...
    'fname_for_out',[],[],...
    'fname_rev_out',[],[],...
    'swapxy_flag',false,[false true],...
    'npad',5,[0,20],...
    'readout_time',1,[0.001,1],...
    'forceflag',false,[false true],...
  ...
    'cleanupflag',true,[false true],...
    'outext','.mgz',{'.mgz','.mgh'},...
  });
  % check that input volumes exist
  if ~exist(parms.fname_for,'file')
    error('file %s not found',parms.fname_for);
  end;
  if ~exist(parms.fname_rev,'file')
    error('file %s not found',parms.fname_rev);
  end;
  % create fname_dx and fstem_dx from fname_for and fname_ref stems
  [tpath,tstem_for,text] = fileparts(parms.fname_for);
  [tpath,tstem_rev,text] = fileparts(parms.fname_rev);
  if isempty(parms.fname_dx)
    parms.fname_dx = sprintf('%s/%sVS%s_B0dx%s',tpath,tstem_for,tstem_rev,parms.outext);
  end;
  [tpath,tstem_dx,text] = fileparts(parms.fname_dx);
  parms.fstem_dx = sprintf('%s/%s',tpath,tstem_dx);
  parms.fname_df = sprintf('%s_fieldcoef.nii',parms.fstem_dx);
  parms.fname_dp = sprintf('%s_movpar.txt',parms.fstem_dx);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_exist = check_output_files(parms)
  all_exist = 1;
  if ~exist(parms.fname_dx,'file')
    all_exist = 0;
    return;
  end;
  if ~exist(parms.fname_df,'file')
    all_exist = 0;
    return;
  end;
  if ~exist(parms.fname_dp,'file')
    all_exist = 0;
    return;
  end;
  if ~isempty(parms.fname_for_out) ~exist(parms.fname_for_out,'file')
    all_exist = 0;
    return;
  end;
  if ~isempty(parms.fname_rev_out) ~exist(parms.fname_rev_out,'file')
    all_exist = 0;
    return;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_inputs_match(parms)
  % create tmp directory
  [tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_dx);
  if isempty(tmp_path), tmp_path = pwd; end;
  parms.tmpdir = sprintf('%s/tmp_topup_estimateB0_%s',tmp_path,tmp_fstem);
  mmil_mkdir(parms.tmpdir);
  % compare info about fname_for and fname_rev
  [dimsmatch,orientmatch,voxmatch,rotmatch,offmatch,...
    parms.M,parms.M_rev,parms.volsz,parms.volsz_rev] =...
      fs_vol_match(parms.fname_for,parms.fname_rev,parms.tmpdir,parms.forceflag);
  if ~dimsmatch
    error('mismatch in input volume dimensions for %s and %s',...
      parms.fname_for,parms.fname_rev);
  end;
  if ~orientmatch
    error('mismatch in volume orientation for %s and %s',...
      parms.fname_for,parms.fname_rev);
  end;
  if ~rotmatch
    error('mismatch in voxel to RAS rotation components for %s and %s',...
      parms.fname_for,parms.fname_rev);
  end;
  if ~offmatch
    error('mismatch in voxel to RAS offset components for %s and %s',...
      parms.fname_for,parms.fname_rev);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_inputs(parms)
  % concat forward and reverse images
  fname_concat = sprintf('%s/for_rev_concat.mgh',parms.tmpdir);
  if ~exist(fname_concat,'file') || parms.forceflag
    % load first frame of forward and reverse images
    vol_for = fs_load_mgh(parms.fname_for,[],1);
    vol_rev = fs_load_mgh(parms.fname_rev,[],1);
    % pad slices with zeros (to avoid wrap-around artifact)
    if parms.npad
      npad_total = 2*parms.npad;
      % add an extra slice if number of slices is odd
      if mod(parms.volsz(3),2)==1
        npad_total = npad_total + 1;
      end;
      tmp_vol_for = zeros(parms.volsz(1),parms.volsz(2),parms.volsz(3)+npad_total);
      tmp_vol_rev = zeros(parms.volsz(1),parms.volsz(2),parms.volsz(3)+npad_total);
      tmp_vol_for(:,:,parms.npad+1:parms.volsz(3)+parms.npad) = vol_for;
      tmp_vol_rev(:,:,parms.npad+1:parms.volsz(3)+parms.npad) = vol_rev;
      vol_for = tmp_vol_for;
      vol_rev = tmp_vol_rev;
    end;
    % concatenate 3D matrices
    vol_concat = cat(4,vol_for,vol_rev);
    % save concatenated images to fname_concat
    fs_save_mgh(vol_concat,fname_concat,parms.M);
  end;
  % check input handedness
  % NOTE: FSL flips the x axis of images stored with a positive determinant,
  %   so topup's outputs are flipped back for those. Keying on the first
  %   orientation letter treated LAS (negative determinant, not flipped by FSL)
  %   like LPS, and mirrored the field left-right for LAS input
  parms.flipflag = (det(parms.M(1:3,1:3)) > 0);
  % convert to nifti
  parms.fname_input = sprintf('%s/for_rev_concat.nii',parms.tmpdir);
  fs_mri_convert(fname_concat,parms.fname_input,'forceflag',parms.forceflag);
  % create acquisition polarity file
  % NOTE: the readout time sets the field's units (Hz) and so the weight of
  %   b02b0.cnf's regularisation, which is in field units; convert_output
  %   multiplies the field by it to give displacements in voxels
  parms.fname_acqp = sprintf('%s/acqp.txt',parms.tmpdir);
  if ~exist(parms.fname_acqp,'file') || parms.forceflag
    fid = fopen(parms.fname_acqp,'wt');
    if fid<0, error('failed ot open file %s for writing\n',parms.fname_acqp); end;
    if parms.swapxy_flag
      fprintf(fid,'1 0 0 %g\n',parms.readout_time);
      fprintf(fid,'-1 0 0 %g\n',parms.readout_time);
    else
      fprintf(fid,'0 1 0 %g\n',parms.readout_time);
      fprintf(fid,'0 -1 0 %g\n',parms.readout_time);
    end;
    fclose(fid);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = estimate_field(parms)
  % estimate displacement field and produce corrected images
  parms.tmp_fstem_dx = sprintf('%s/dx',parms.tmpdir);
  parms.tmp_fname_dx = sprintf('%s_field.nii',parms.tmp_fstem_dx);
  parms.tmp_fname_df = sprintf('%s_fieldcoef.nii',parms.tmp_fstem_dx);
  parms.tmp_fname_dp = sprintf('%s_movpar.txt',parms.tmp_fstem_dx);
  parms.tmp_fname_out = sprintf('%s/for_rev_concat_B0uw.nii',parms.tmpdir);
  if ~exist(parms.tmp_fname_out,'file') || ~exist(parms.tmp_fname_dx,'file') || parms.forceflag
    % generate command string
    cmd = sprintf('setenv FSLOUTPUTTYPE NIFTI;');
    cmd = sprintf('%s topup -v',cmd);
    cmd = sprintf('%s --imain=%s',cmd,parms.fname_input);
    cmd = sprintf('%s --datain=%s --config=b02b0.cnf',cmd,parms.fname_acqp);
    cmd = sprintf('%s --fout=%s --out=%s --iout=%s',cmd,parms.tmp_fname_dx,parms.tmp_fstem_dx,parms.tmp_fname_out);
    % run command
    fprintf('%s:\ncmd = %s\n',mfilename,cmd);
    [status,result] = mmil_unix(cmd);
    if status
      error('cmd %s failed:\n%s\n',cmd,result);
    else
      disp(result);
    end;
  end;
  if ~exist(parms.tmp_fname_dx,'file')
    error('expected output files (%s) not found',parms.tmp_fname_dx);
  end;
  if ~exist(parms.tmp_fname_df,'file')
    error('expected output files (%s) not found',parms.tmp_fname_df);
  end;
  if ~exist(parms.tmp_fname_dp,'file')
    error('expected output files (%s) not found',parms.tmp_fname_dp);
  end;
  if ~exist(parms.tmp_fname_out,'file')
    error('expected output files (%s) not found',parms.tmp_fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convert_output(parms)
  % convert output back to mgh
  fname_out_tmp = sprintf('%s/for_rev_concat_B0uw_tmp.mgh',parms.tmpdir);
  fs_mri_convert(parms.tmp_fname_out,fname_out_tmp,'forceflag',parms.forceflag);
  fname_dx_tmp = sprintf('%s/dx_tmp.mgh',parms.tmpdir);
  fs_mri_convert(parms.tmp_fname_dx,fname_dx_tmp,'forceflag',parms.forceflag);
  fname_df = parms.tmp_fname_df;
  fname_dp = parms.tmp_fname_dp;
  % split for and rev
  fname_for_out = sprintf('%s/for_B0uw.mgh',parms.tmpdir);
  fname_rev_out = sprintf('%s/rev_B0uw.mgh',parms.tmpdir);
  if ~exist(fname_for_out,'file') || ~exist(fname_rev_out,'file') || parms.forceflag
    % load concatenated output file
    vol = fs_load_mgh(fname_out_tmp);
    % left-right flip of images
    if parms.flipflag, vol = vol(end:-1:1,:,:,:); end;
    % remove zero padded slices
    if parms.npad
      vol = vol(:,:,parms.npad+1:parms.volsz(3)+parms.npad,:);
    end;
    % save for and rev to separate files
    fs_save_mgh(vol(:,:,:,1),fname_for_out,parms.M);
    fs_save_mgh(vol(:,:,:,2),fname_rev_out,parms.M);
  end;
  % adjust dx volume
  fname_dx = sprintf('%s/dx.mgh',parms.tmpdir);
  if ~exist(fname_dx,'file') || parms.forceflag
    % load concatenated output file
    vol = fs_load_mgh(fname_dx_tmp);
    % left-right flip of images
    if parms.flipflag, vol = vol(end:-1:1,:,:); end;
    % remove zero padded slices
    if parms.npad
      vol = squeeze(vol(:,:,parms.npad+1:parms.volsz(3)+parms.npad));
    end;
    % field (Hz) to displacement (voxels), with the readout time in acqp.txt
    vol = vol * parms.readout_time;
    % save to file
    fs_save_mgh(vol,fname_dx,parms.M);
  end;
  % copy tmp output files
  flist1 = {fname_dx,fname_df,fname_dp,fname_for_out,fname_rev_out};
  flist2 = {parms.fname_dx,parms.fname_df,parms.fname_dp,parms.fname_for_out,parms.fname_rev_out};
  for f=1:length(flist1)
    if isempty(flist2{f}), continue; end;
    if ~exist(flist1{f},'file'), continue; end;
    [~,~,ext1] = cmig_get_fileparts(flist1{f});
    [~,~,ext2] = cmig_get_fileparts(flist2{f});
    if strcmp(ext1,ext2)
      cmd = sprintf('cp -p %s %s',flist1{f},flist2{f});
      [status,result,errmsg] = mmil_unix(cmd);
      if status
        error('cmd %s failed:\n%s\n%s',cmd,result,errmsg);
      end;
    elseif ismember(ext2,{'.nii','.nii.gz'})
      % NOTE: fs_copy_mgh writes MGH whatever the extension
      fs_mri_convert(flist1{f},flist2{f},'forceflag',1);
    else
      fs_copy_mgh(flist1{f},flist2{f});
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanupflag
    % remove tmp dir
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,result,errmsg] = mmil_unix(cmd);
    if status
      fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n%s\n%s\n',mfilename,cmd,result,errmsg);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
