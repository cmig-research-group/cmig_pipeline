function cmig_write_nifti(vol,fname,M,hdr,keeptype_flag)
%function cmig_write_nifti(vol,fname,M,[hdr],[keeptype_flag])
%
% Purpose: write imaging data and header info to NIfTI file
%          requires matlab's built in niftiwrite
%
% Required Input:
%   vol: 3D or 4D volume
%   fname: output file name (.nii or .nii.gz)
%   M: 4x4 vox2ras transform such that
%     y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the indices are 1-based
%
% Optional Input:
%   hdr: header information that is not carried by vol and M, either the
%     struct returned by cmig_read_nifti or, for backward compatibility, a
%     scalar giving the TR; recognized fields, all optional:
%     TR: time step, written as the 4th pixel dimension, in seconds
%     pixdim: the three spatial pixel dimensions; used in place of the column
%       norms of M when the two agree to within 1e-04 relative, and ignored
%       otherwise, so a resampled volume still gets dimensions matching M
%     descrip: free-text header field, truncated to the 80 bytes the header
%       holds; niftiwrite rejects anything longer
%     scl_slope, scl_inter: intensity scaling, recorded in the header and NOT
%       applied to vol; cmig_read_nifti reports both as 0 when it has already
%       applied the scaling, so passing its hdr straight back is correct
%       without altering the values in vol
%     freq_dim, phase_dim, slice_dim: frequency, phase and slice encoding
%       directions, i.e. the dim_info byte
%     {default = [] (all zero, matching mri_convert output)}
%   keeptype_flag: [0|1] write vol in the class it was passed in, rather than
%     casting to single
%     {default = 0}
%
% Notes:
%   vol is cast to single unless keeptype_flag is set, so the output data type
%     is float32 by default, matching fs_save_mgh, which writes float32
%     whatever it is handed, and fs_save_nifti, which writes float32 because
%     it stores an intermediate mgh file
%   niftiwrite sets the sform only (sform_code = 1, qform_code = 0); this
%     matches how mri_convert output is read by AFNI and cmig_read_nifti,
%     both of which use the sform, but tools that require a qform will not
%     find one
%   a 3D volume with a nonzero TR is declared 4D with a single frame, since
%     niftiwrite requires as many pixel dimensions as the image has
%     dimensions and the TR is the 4th one; the result converts to an mgz
%     identical to converting the input directly, nframes included
%
% See also: cmig_read_nifti, whose scale_flag decides whether vol arrives
%   scaled; this function only records what hdr says and never rescales vol
%
% Created:   06/15/2026 by Don Hagler
% Prev Mod:  08/18/2026 by Don Hagler
% Prev Mod:  08/22/2026 by Don Hagler
% Last Mod:  08/28/2026 by Don Hagler
%

%% Based on fs_save_nifti.m
%   Created:  11/16/17 by Don Hagler
%   Last Mod: 11/17/17 by Don Hagler

if ~mmil_check_nargs(nargin,3), return; end;
if ~exist('hdr','var') | isempty(hdr), hdr = []; end;
if ~exist('keeptype_flag','var') | isempty(keeptype_flag), keeptype_flag = 0; end;

% accept a scalar TR in place of the struct
if isnumeric(hdr)
  TR = hdr; hdr = []; hdr.TR = TR;
end;
TR = mmil_getfield(hdr,'TR',0);
scl_slope = mmil_getfield(hdr,'scl_slope',0);
% descrip is char[80] in the header and niftiwrite errors above 80 rather than
%   truncating, so an over-long string would cost the whole output file
descrip = mmil_getfield(hdr,'descrip','');
if ~ischar(descrip), descrip = ''; end;
if numel(descrip)>80, descrip = descrip(1:80); end;
scl_inter = mmil_getfield(hdr,'scl_inter',0);
freq_dim = mmil_getfield(hdr,'freq_dim',0);
phase_dim = mmil_getfield(hdr,'phase_dim',0);
slice_dim = mmil_getfield(hdr,'slice_dim',0);

if numel(size(M))~=2 | any(size(M)~=[4,4])
  error('M must be a 4x4 matrix');
end;

[fpath,fstem,fext] = cmig_get_fileparts(fname);
if isempty(fpath)
  fpath = pwd;
elseif ~exist(fpath,'dir')
  mmil_mkdir(fpath);
end;

% check if output should be compressed
zipflag = strcmp(fext,'.nii.gz');
% niftiwrite appends the extension itself
fname_out = sprintf('%s/%s',fpath,fstem);

% single, not the class of vol: an integer or double volume would otherwise
%   set the output data type, and cmig_read_nifti returns double
if ~keeptype_flag
  vol = single(vol);
end

% get image dimensions
volsz = size(vol);

% voxel sizes are the column norms of M, with TR for the 4th dimension
pixdim = sqrt(sum(M(1:3,1:3).^2,1));

% prefer the stored pixel dimensions when hdr has them and they agree with the
%   derived ones. The stored float32 direction cosines are not exactly unit
%   norm, so the derived sizes land parts per million off, which shows up as a
%   spurious header difference downstream. The agreement test is what makes
%   this safe: a caller that resampled and passed the original hdr would
%   otherwise write pixel dimensions contradicting M
hdr_pixdim = mmil_getfield(hdr,'pixdim',[]);
if numel(hdr_pixdim)==3 && all(hdr_pixdim>0)
  if all(abs(hdr_pixdim(:)'-pixdim)./pixdim < 1e-4)
    pixdim = hdr_pixdim(:)';
  end;
end;
if numel(volsz)>3
  pixdim = [pixdim,TR*ones(1,numel(volsz)-3)];
elseif TR~=0
  % niftiwrite requires numel(PixelDimensions) == numel(ImageSize), so a 3D
  %   volume carrying a TR is declared 4D with a single frame
  volsz = [volsz,1];
  pixdim = [pixdim,TR];
end;

% convert from 1-based to 0-based indexing
M = M*[1 0 0 1; 0 1 0 1; 0 0 1 1; 0 0 0 1];

% construct the header
% NOTE: niftiinfo requires an existing file, so the struct is built directly
info = struct(...
  'Filename','',...
  'Filemoddate','',...
  'Filesize',0,...
  'Version','NIfTI1',...
  'Description',descrip,...
  'ImageSize',volsz,...
  'PixelDimensions',pixdim,...
  'Datatype',class(vol),...
  'BitsPerPixel',8*numel(typecast(cast(0,class(vol)),'uint8')),...
  'SpaceUnits','Millimeter',...
  'TimeUnits','Second',...
  'AdditiveOffset',scl_inter,...
  'MultiplicativeScaling',scl_slope,...
  'TimeOffset',0,...
  'SliceCode','Unknown',...
  'FrequencyDimension',freq_dim,...
  'PhaseDimension',phase_dim,...
  'SpatialDimension',slice_dim,...
  'DisplayIntensityRange',[0,0],...
  'TransformName','Sform',...
  'Transform',affine3d(M'),...
  'Qfactor',1);

% write the image to a .nii or .nii.gz file
niftiwrite(vol,fname_out,info,'Compressed',zipflag);

return;
