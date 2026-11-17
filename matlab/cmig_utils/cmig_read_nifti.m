function [vol,M,volsz,hdr] = cmig_read_nifti(fname,headeronly,scale_flag)
%function [vol,M,volsz,hdr] = cmig_read_nifti(fname,[headeronly],[scale_flag])
%
% Purpose: read imaging data and header info from NIfTI file
%          requires matlab's built in niftiread
%
% Required Input:
%   fname: path of the nifti file
%
% Optional Input:
%   headeronly: [0|1] whether to load header info only
%     if 1, vol will be empty
%     {default = 0}
%   scale_flag: [0|1] apply the intensity scaling in the header, i.e. return
%     vol*scl_slope + scl_inter rather than the stored values
%     also casts vol to double, whether or not the header asks for scaling
%     {default = 1}
%
% Output:
%   vol: 4D volume (if surface data, size will be [n,1,1,f])
%   M: 4x4 vox2ras transform such that
%     y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the indices are 1-based
%   volsz: vector of matrix dimensions, always with four elements
%     niftiinfo reports only the dimensions present, so trailing singleton
%     dimensions are added, matching the fs_load_mgh convention
%   hdr: struct of header information that is not carried by vol and M,
%     for passing to cmig_write_nifti so that a round trip preserves it
%     TR: time step, i.e. the 4th pixel dimension, in seconds
%     pixdim: the three spatial pixel dimensions as stored, so the writer can
%       carry them instead of recovering them from the column norms of M
%     descrip: the free-text header field, e.g. the dcm2niix acquisition stamp
%     scl_slope, scl_inter: intensity scaling
%     freq_dim, phase_dim, slice_dim: frequency, phase and slice encoding
%       directions, i.e. the dim_info byte
%
% NOTE: with scale_flag = 1, the default, vol carries physical intensities, as
%   mri_convert, AFNI, FSL and nibabel's get_fdata all return; matlab's
%   niftiread does not apply the scaling, so this departs from it deliberately.
%   ABCD Philips data carries a slope of order 100 for anat and dwi and of
%   order 10 for bold, while Siemens and GE carry 1, so reading unscaled makes
%   absolute intensities incomparable across vendors
%
% NOTE: when the scaling is applied, the returned hdr reports scl_slope = 0 and
%   scl_inter = 0, so writing vol back out with cmig_write_nifti does not apply
%   it a second time. That is what mri_convert records too. With headeronly = 1
%   no data is returned, so hdr reports the values actually in the file
%
% NOTE: with scale_flag = 1 vol is cast to double whether or not the header
%   asks for scaling, matching fs_load_mgh, which reads every integer datatype
%   into double. An integer volume cannot hold the scaled values, and matlab
%   refuses to combine an integer array with a double one at all, so a binary
%   mask read from a slope-1 file would otherwise come back uint8 and break
%   arithmetic against a double volume. With scale_flag = 0 the stored class
%   survives
%
% Created:  06/15/2026 by Don Hagler
% Prev Mod: 08/22/2026 by Don Hagler
% Last Mod: 08/28/2026 by Don Hagler
%

vol = []; M = []; volsz = []; hdr = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('headeronly','var') | isempty(headeronly), headeronly=0; end;
if ~exist('scale_flag','var') | isempty(scale_flag), scale_flag=1; end;

% read header
info = niftiinfo(fname);

% get voxel to RAS matrix
M = info.Transform.T';

% Convert from 0-based to 1-based indexing
M = M*[1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];

% get volume size (number of voxels in each dimension)
volsz = info.ImageSize;

% niftiinfo reports only the dimensions present, so pad to four
if numel(volsz)<4, volsz = [volsz(:)',ones(1,4-numel(volsz))]; end;

% header information not carried by vol and M
% NOTE: read from info.raw because niftiinfo reports only as many pixel
%   dimensions as the file declares, so the TR of a 3D volume is not in
%   info.PixelDimensions even when the file carries one
hdr = [];
hdr.TR = 0;
if numel(info.raw.pixdim)>4
  hdr.TR = double(info.raw.pixdim(5));
end;
hdr.pixdim = double(info.raw.pixdim(2:4));
hdr.descrip = info.Description;
hdr.scl_slope = double(info.raw.scl_slope);
hdr.scl_inter = double(info.raw.scl_inter);
hdr.freq_dim = info.FrequencyDimension;
hdr.phase_dim = info.PhaseDimension;
hdr.slice_dim = info.SpatialDimension;

% if only need header quit now
if headeronly, return; end

% read image data
vol = niftiread(info);

% apply the intensity scaling recorded in the header
%   the cast comes first and is unconditional: an integer times a double stays
%   integer in matlab, so scaling in place would round and saturate at intmax
%   slope 0 means no scaling, per the NIfTI spec, and slope 1 with a zero
%   offset is a no-op, so the multiply itself is skipped for both, leaving
%   hdr.scl_slope and hdr.scl_inter as stored
if scale_flag
  vol = double(vol);
  if (hdr.scl_slope~=0 && hdr.scl_slope~=1) || hdr.scl_inter~=0
    vol = vol*hdr.scl_slope + hdr.scl_inter;
    hdr.scl_slope = 0;
    hdr.scl_inter = 0;
  end
end;


