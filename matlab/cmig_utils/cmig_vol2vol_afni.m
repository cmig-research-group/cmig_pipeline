function M_synth_to_orig = cmig_vol2vol_afni(vol1,vol2,tmpdir,fname_weight)
%function M_synth_to_orig = cmig_vol2vol_afni(vol1,vol2,tmpdir,[fname_weight])
%
% Purpose: rigid-body registration between two volumes using AFNI's
%   3dvolreg, returning the transform in the same convention as
%   rbreg_vol2vol_icorr_mask
%
% Required Input:
%   vol1: reference volume (ctx struct), passed to 3dvolreg as -base
%   vol2: volume registered to vol1 (ctx struct)
%   tmpdir: directory for the temporary nifti files, which are deleted
%     before returning
%
% Optional Input:
%   fname_weight: nifti file with voxel weights for the registration cost
%     (same grid as vol1); passed to 3dvolreg -weight
%     {default = []}
%
% Output:
%   M_synth_to_orig: 4x4 registration matrix from vol1 to vol2
%
% NOTE: the vox2ras matrices are taken from the ctx structs, so vol1 and
%   vol2 must carry them; 3dvolreg reports rotations about the volume grid
%   center rather than the scanner-coordinate origin, which is compensated
%   for below
%
% NOTE: 3dvolreg is called with -twopass; see the comment at the command
%
% See also: rbreg_vol2vol_icorr_mask, mmil_load_motion_1D
%
% Adapted from the 2023-2025 prototype by Mai Bui
%
% Created:  08/12/2026 by Don Hagler
% Prev Mod: 08/19/2026 by Don Hagler
% Prev Mod: 08/24/2026 by Don Hagler
% Prev Mod: 08/25/2026 by Don Hagler
% Prev Mod: 08/26/2026 by Don Hagler
% Prev Mod: 09/21/2026 by Don Hagler
% Last Mod: 10/13/2026 by Don Hagler
%

  if ~mmil_check_nargs(nargin,3), return; end;
  if ~exist('fname_weight','var') | isempty(fname_weight)
    fname_weight = [];
  end;

  % write the reference and input volumes as nifti files for 3dvolreg
  % NOTE: cmig_write_nifti rather than fs_save_nifti, which writes an
  %   intermediate mgh and runs mri_convert. Cast to single because
  %   fs_save_nifti wrote float32
   
  [V1,M1] = ctx_ctx2mgh(vol1);
  fstem1 = mmil_tempfname('vol1',tmpdir);
  fname_ref = [fstem1 '.nii'];
  cmig_write_nifti(single(V1),fname_ref,M1);
  
  [V2,M2] = ctx_ctx2mgh(vol2);
  fstem2 = mmil_tempfname('vol2',tmpdir);
  fname_in = [fstem2 '.nii'];
  cmig_write_nifti(single(V2),fname_in,M2);
  
  fname_motion_1D = [fstem2 '_motion.1D'];   % -dfile option for output 9 ASCII formated columns
  
  % register images
  % NOTE: -twopass is required, not an optimization. Without it large
  %   rotations collapse toward zero and read as LOW motion, so a failed
  %   registration looks like a quiet scan. Unconditional: there is no safe
  %   angle to gate on
  fname_out = [fstem2 '_mc.nii'];
  % AFNI_DONT_LOGFILE keeps 3dvolreg from appending to ~/.afni.log
  cmd = sprintf('env AFNI_DONT_LOGFILE=YES 3dvolreg -twopass -prefix %s',fname_out);
  % kernel named: the default changed to heptic in newer AFNI;
  %   -zpad 4 keeps Fourier shifts from wrapping across the edge slices
  cmd = sprintf('%s -Fourier -zpad 4',cmd);
  cmd = sprintf('%s -base ''%s''',cmd,fname_ref);
  if ~isempty(fname_weight)
    % sub-brick selector is required by 3dvolreg -weight
    cmd = sprintf('%s -weight ''%s[0]''',cmd,fname_weight);
  end
  cmd = sprintf('%s -dfile %s',cmd,fname_motion_1D);
  cmd = sprintf('%s %s',cmd,fname_in);

  % mmil_unix, not jsystem: jsystem defaults to /bin/sh, which does not source
  %   the csh setup that puts AFNI on the PATH, so 3dvolreg is not found.
  %   mmil_unix passes tcsh explicitly and returns stderr the same way
  [status,result,errmsg] = mmil_unix(cmd);
  if status
    error('3dvolreg failed:\n%s\n%s',result,errmsg);
  end
  
  % load motion_1D file (output from 3dvolreg)
  motion_tseries_afni = mmil_load_motion_1D(fname_motion_1D);
  tmp_motion = motion_tseries_afni(1,:);
  trans = tmp_motion(1:3);
  rot   = tmp_motion(4:6);
  M_reg_inv = mmil_construct_M('trans',trans,'rot',rot,...
      'order',{'rot','trans'},'rot_order',{'z','x','y'});

  % 3dvolreg rotations are about the volume grid center, not the
  % scanner-coordinate origin, so the reported translations require a
  % (I-R)*c compensation term, with c the grid center in LPH coordinates
  % (verified synthetically: test_mat/test_afni_transform_convention.m
  % and test_afni_fit_rotation_center.m, 08/12/26)
  volsz = size(V2);
  ctr = M2*[(volsz(1)+1)/2; (volsz(2)+1)/2; (volsz(3)+1)/2; 1];
  c_lph = [-ctr(1); -ctr(2); ctr(3)];
  R = M_reg_inv(1:3,1:3);
  M_reg_inv(1:3,4) = M_reg_inv(1:3,4) + (eye(3)-R)*c_lph;

  M_reg = inv(M_reg_inv);

  M_synth_to_orig = M_reg;

  % delete temporary files
  delete(fname_in);
  delete(fname_ref);
  delete(fname_motion_1D);
  delete(fname_out);

return;
