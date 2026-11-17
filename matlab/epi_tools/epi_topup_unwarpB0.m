function epi_topup_unwarpB0(fname_in,fstem_dx,varargin)
%function epi_topup_unwarpB0(fname_in,fstem_dx,[options])
%
% Purpse: apply displacement field to correct
%   distortions caused by B0 susceptibility artifacts
%   using FSL's applytopup
%
% Usage: epi_topup_unwarpeB0(fname_in,fstem_dx,'key','value'...)
%
% Required Parameters:
%   fname_in:    input file name
%   fstem_dx:    input file stem for fieldcoef, movpar
%
% Optional Parameters:
%  'fname_out': output file name (B0 unwarped)
%     If not supplied, will append 'B0uw' to file stem of fname_in
%    {default = []}
%  'revflag': [0|1] whether to reverse direction of displacements
%    {default = 0}
%  'swapxy_flag': [0|1] whether to swap x and y dimensions before unwarping B0
%    {default = 0}
%  'readout_time': total readout time (s) of fname_in, written to
%    applytopup's acquisition file; the field is in Hz when
%    epi_topup_estimateB0 was given the true readout time, and the
%    displacement applied is that field times readout_time
%    {default = 1}
%  'verbose': [0|1] display status messages
%    {default = 0}
%  'forceflag': [0|1] whether to overwrite existing fname_out
%    {default = 0}
%
% Created:  07/06/19 by Don Hagler
% Prev Mod: 07/01/20 by Don Hagler
% Prev Mod: 07/04/20 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Last Mod: 10/09/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on epi_unwarpB0
% Created:  12/19/08 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input parameters
if ~mmil_check_nargs(nargin,2), return; end;
parms = check_input_parms(fname_in,fstem_dx,varargin);

% check output file
if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% check that input volumes share dimensions and orientations
parms = check_inputs_match(parms);

% prepare input files
parms = prep_inputs(parms);

% apply field
parms = apply_field(parms);

% convert output file
convert_output(parms);

% remove temporary files
cleanup_tmpdir(parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input_parms(fname_in,fstem_dx,options)
  parms = mmil_args2parms(options,{...
    'fname_in',fname_in,[],...
    'fstem_dx',fstem_dx,[],...
  ...
    'fname_out',[],[],...
    'revflag',false,[false true],...
    'swapxy_flag',false,[false true],...
    'npad',5,[0,20],...
    'readout_time',1,[0.001,1],...
    'verbose',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'cleanupflag',true,[false true],...
    'outext','.mgz',{'.mgz','.mgh'},...
  });
  % check that input volumes exist
  if ~exist(parms.fname_in,'file')
    error('input file %s not found',parms.fname_in);
  end;
  parms.fname_dx = sprintf('%s.mgz',parms.fstem_dx);
  parms.fname_df = sprintf('%s_fieldcoef.nii',parms.fstem_dx);
  parms.fname_dp = sprintf('%s_movpar.txt',parms.fstem_dx);
  if ~exist(parms.fname_df,'file')
    error('displacement field coefficients file %s not found',parms.fname_df);
  end;
  if ~exist(parms.fname_dp,'file')
    error('movement parameters file %s not found',parms.fname_dp);
  end;
  % set fname_out
  if isempty(parms.fname_out)
    [tpath,tstem,text] = fileparts(fname_in);
    parms.fname_out = [tpath '/' tstem '_B0uw' parms.outext];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_inputs_match(parms)
  % create tmp directory
  [tmp_path,tmp_fstem,tmp_ext] = fileparts(parms.fname_out);
  if isempty(tmp_path), tmp_path = pwd; end;
  parms.tmpdir = sprintf('%s/tmp_topup_unwarpB0_%s',tmp_path,tmp_fstem);
  mmil_mkdir(parms.tmpdir);
  % compare info about fname_in and fname_dx
  [dimsmatch,orientmatch,voxmatch,rotmatch,offmatch,...
    parms.M,parms.M_dx,parms.volsz,parms.volsz_dx] =...
      fs_vol_match(parms.fname_in,parms.fname_dx,parms.tmpdir,parms.forceflag);
  if ~dimsmatch
    error('mismatch in input volume dimensions for %s and %s',...
      parms.fname_in,parms.fname_dx);
  end;
  if ~orientmatch
    error('mismatch in volume orientation for %s and %s',...
      parms.fname_in,parms.fname_dx);
  end;
  %% todo: require voxmatch if this is forward scan but allow small offset if reverse?
  %%       also require that M_rev is same as M_ref_rev --> but is that info available?
%    if ~voxmatch
%      error('mismatch in voxel to RAS transformation for %s and %s',...
%        parms.fname_in,parms.fname_dx);
%    end;
  if ~rotmatch
    error('mismatch in voxel to RAS rotation components for %s and %s',...
      parms.fname_in,parms.fname_dx);
  end;
  if ~offmatch
    error('mismatch in voxel to RAS offset components for %s and %s',...
      parms.fname_in,parms.fname_dx);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = prep_inputs(parms)
  % save as nifti
  parms.tmp_fname_in = sprintf('%s/input.nii',parms.tmpdir);
  if ~exist(parms.tmp_fname_in,'file') || parms.forceflag
    mmil_mkdir(parms.tmpdir);
    % load input image and metadata
    vol = fs_load_mgh(parms.fname_in);
    nframes = parms.volsz(4);
    % add zero padded slices
    if parms.npad
      npad_total = 2*parms.npad;
      % add an extra slice if number of slices is odd
      if mod(parms.volsz(3),2)==1
        npad_total = npad_total + 1;
      end;
      tmp_vol = zeros(parms.volsz(1),parms.volsz(2),parms.volsz(3)+npad_total,nframes);
      tmp_vol(:,:,parms.npad+1:parms.volsz(3)+parms.npad,:) = vol;
      vol = tmp_vol;
      clear tmp_vol;
    end;
    %% todo: change fs_save_nifti to write directly to nifti instead of writing to mgh and converting to nifti
    fs_save_nifti(vol,parms.tmp_fname_in,parms.M);
  end;
  % create acquisition polarity file
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

% apply field to correct distortions
function parms = apply_field(parms)
  % apply displacement field to produce corrected images
  parms.tmp_fstem_out = sprintf('%s/output_B0uw',parms.tmpdir);
  parms.tmp_fname_out = sprintf('%s.nii',parms.tmp_fstem_out);
  if ~exist(parms.tmp_fname_out,'file') || parms.forceflag
    % generate command string
    cmd = sprintf('setenv FSLOUTPUTTYPE NIFTI;');
    cmd = sprintf('%s applytopup -v',cmd);
    cmd = sprintf('%s --imain=%s',cmd,parms.tmp_fname_in);
    cmd = sprintf('%s --datain=%s',cmd,parms.fname_acqp);
    cmd = sprintf('%s --inindex=%d',cmd,parms.revflag+1);
    cmd = sprintf('%s --topup=%s',cmd,parms.fstem_dx);
    cmd = sprintf('%s --method=jac',cmd);
    cmd = sprintf('%s --out=%s',cmd,parms.tmp_fstem_out);
    % run command
    if parms.verbose
      fprintf('%s:\ncmd = %s\n',mfilename,cmd);
    end;
    [status,result,errmsg] = mmil_unix(cmd);
    if status
      error('cmd %s failed:\n%s\n%s',cmd,result,errmsg);
    elseif parms.verbose
      disp(result);
    end;
  end;
  % check that expected output file was created
  if ~exist(parms.tmp_fname_out,'file')
    error('expected output files (%s) not found',parms.tmp_fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert output files
function convert_output(parms)
  % remove zero padded slices
  if parms.npad
    % load output file
    tmp_fname = sprintf('%s.mgz',mmil_tempfname('output_B0uw',parms.tmpdir));
    fs_mri_convert(parms.tmp_fname_out,tmp_fname);
    [vol,M] = fs_load_mgh(tmp_fname);
    delete(tmp_fname);
    %% todo: change fs_load_nifti to adjust M matrix position vector
    %%       to match output of fs_load_mgh
%    [vol,M] = fs_load_nifti(parms.tmp_fname_out);
    vol = vol(:,:,parms.npad+1:parms.volsz(3)+parms.npad,:);
    % save as fname_out
    % NOTE: fs_save_mgh writes MGH whatever the extension, so a nifti name
    %   is written through a temporary mgz
    [~,~,fext] = cmig_get_fileparts(parms.fname_out);
    if ismember(fext,{'.nii','.nii.gz'})
      tmp_fname = sprintf('%s.mgz',mmil_tempfname('output_B0uw',parms.tmpdir));
      fs_save_mgh(vol,tmp_fname,M);
      fs_mri_convert(tmp_fname,parms.fname_out,'forceflag',1);
      delete(tmp_fname);
    else
      fs_save_mgh(vol,parms.fname_out,M);
    end;
  else
    fs_mri_convert(parms.tmp_fname_out,parms.fname_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cleanup_tmpdir(parms)
  if parms.cleanupflag
    % remove tmp dir
    cmd = sprintf('rm -r %s',parms.tmpdir);
    [status,result,errmsg] = mmil_unix(cmd);
    if status
      fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n%s\n',mfilename,result,errmsg);
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
