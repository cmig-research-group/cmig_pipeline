function [b0_rms,b0_rms_vals,fnames_B0dx] = cmig_calc_b0dx_rms(proc_rootdir,participant_id,session_id,acq_type,tmpdir,varargin)
%function [b0_rms,b0_rms_vals,fnames_B0dx] = cmig_calc_b0dx_rms(proc_rootdir,participant_id,session_id,acq_type,tmpdir,[options])
%
% Purpose: compute RMS of B0 unwarp displacement field(s) for a session
%
%   Averages the RMS of every B0dx.nii.gz file found for the given
%   acq_type in {proc_rootdir}/{participant_id}/{session_id}/fmap;
%   there is no BIDS-native equivalent of source's DICOM-era SessInfo
%   bookkeeping to identify a single designated forward/reverse
%   reference pair, so all matching B0dx files are averaged instead
%
%   The displacement fields are written by cmig_proc_bids_func in mm, in
%   the reference field map's frame at output resolution, so fname_mask
%   must be a brain mask on that same grid
%
% Required Input:
%   proc_rootdir: BIDS processed derivatives root directory
%   participant_id: participant identifier string
%   session_id: session identifier string
%   acq_type: acquisition type string used in the B0dx filename's acq-
%     entity (e.g. 'func', 'dwi')
%   tmpdir: directory for caching the computed result
%
% Optional Input ('key',value):
%   'fname_mask': full path of a brain mask file (mgh, mgz, nii, or nii.gz),
%     in the same voxel grid as the B0dx displacement field(s), used to
%     restrict the RMS calculation to brain tissue
%     if empty or not supplied, RMS is computed over the whole volume
%     a non-empty file name that does not exist, an unsupported file
%     extension, or a mask whose grid does not match the displacement
%     field is an error, since each can only arise from a code defect
%     {default = []}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing cached output
%     {default = 0}
%
% Output:
%   b0_rms: mean RMS of the B0 displacement field(s) in mm; NaN if no B0dx
%     file was found, none could be read, or the brain mask was empty
%   b0_rms_vals: vector of per-pair RMS values, one per B0dx file
%   fnames_B0dx: cell array of B0dx file names, matching b0_rms_vals
%
% Created:  07/10/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Last Mod: 08/05/2026 by Don Hagler
%

b0_rms = NaN; b0_rms_vals = []; fnames_B0dx = {};
if ~mmil_check_nargs(nargin,5), return; end

parms = mmil_args2parms(varargin,{...
  'fname_mask',[],[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
});

fmap_dir = sprintf('%s/%s/%s/fmap',proc_rootdir,participant_id,session_id);
dlist = dir(sprintf('%s/*_acq-%s_desc-AP*PA*_B0dx.nii.gz',fmap_dir,acq_type));
if isempty(dlist)
  fprintf('%s: WARNING: no B0dx file found for acq-%s in %s; b0_rms set to NaN\n',...
    mfilename,acq_type,fmap_dir);
  return;
end
fnames_B0dx = cell(length(dlist),1);
for k=1:length(dlist)
  fnames_B0dx{k} = sprintf('%s/%s',fmap_dir,dlist(k).name);
end

fname_out = sprintf('%s/%s_B0dx_RMS.mat',tmpdir,acq_type);
if ~exist(fname_out,'file') || parms.forceflag
  if parms.verbose
    fprintf('%s: computing RMS of B0 displacement field...\n',mfilename);
  end
  % load brain mask
  mask = [];
  if ~isempty(parms.fname_mask)
    if ~exist(parms.fname_mask,'file')
      error('%s: brain mask %s not found',mfilename,parms.fname_mask);
    end
    [~,~,fext] = cmig_get_fileparts(parms.fname_mask);
    switch fext
      case {'.mgh','.mgz'}
        vol_mask = fs_load_mgh(parms.fname_mask);
      case {'.nii','.nii.gz'}
        vol_mask = cmig_read_nifti(parms.fname_mask);
      otherwise
        error('%s: unsupported file format %s for brain mask %s',...
          mfilename,fext,parms.fname_mask);
    end
    mask = (vol_mask > 0);
    % an empty mask is a plausible data condition on artifact-heavy sessions,
    %   which are exactly the ones this measure exists to flag
    if ~any(mask(:))
      fprintf('%s: WARNING: brain mask %s contains no voxels; b0_rms set to NaN\n',...
        mfilename,parms.fname_mask);
      return;
    end
  end
  b0_rms_vals = nan(length(dlist),1);
  for k=1:length(dlist)
    vol = [];
    try
      vol = cmig_read_nifti(fnames_B0dx{k});
    catch me
      fprintf('%s: WARNING: cmig_read_nifti failed for %s: %s\n',...
        mfilename,dlist(k).name,me.message);
    end
    if isempty(vol), continue; end
    if ndims(vol) > 3, vol = vol(:,:,:,1); end
    if isempty(mask)
      b0_rms_vals(k) = rms(vol(:));
    else
      if ~isequal(size(mask),size(vol))
        error('%s: brain mask %s ([%s]) does not match B0dx %s ([%s])',...
          mfilename,parms.fname_mask,strtrim(sprintf('%d ',size(mask))),...
          dlist(k).name,strtrim(sprintf('%d ',size(vol))));
      end
      b0_rms_vals(k) = rms(vol(mask));
    end
  end
  tmp_vals = b0_rms_vals(~isnan(b0_rms_vals));
  if ~isempty(tmp_vals)
    b0_rms = mean(tmp_vals);
  end
  save(fname_out,'b0_rms','b0_rms_vals','fnames_B0dx');
else
  tmp = load(fname_out);
  b0_rms = tmp.b0_rms;
  b0_rms_vals = mmil_getfield(tmp,'b0_rms_vals',nan(length(dlist),1));
  fnames_B0dx = mmil_getfield(tmp,'fnames_B0dx',fnames_B0dx);
end

return;
