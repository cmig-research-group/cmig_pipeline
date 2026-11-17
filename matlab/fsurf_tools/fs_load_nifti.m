function [vol, M, volsz] = fs_load_nifti(fname,headeronly)
%function [vol, M, volsz] = fs_load_nifti(fname,headeronly)
% hdr = fs_load_nifti(fname,headeronly)
%
% Required Input:
%   fname: path of the nifti file
%
% Optional Input:
%   headeronly: [0|1] whether to load header info only
%     if 1, vol will be empty
%     {default = 0}
%
% Output:
%   vol: 4D volume (if surface data, size will be [n,1,1,f])
%   M: 4x4 vox2ras transform such that
%     y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
%     indices are 1-based. 
%     If the input has multiple frames,
%     only the first frame is read.
%   mr_parms: [tr flipangle te ti fov]
%   volsz: size(vol). Helpful when using headeronly as vol is [].
%
% Loads nifti header and volume. Columns and rows are not swapped.
%
% Handles compressed nifti (nii.gz) by issuing a unix command to
% uncompress the file to a temporary file, which is then deleted.
%
% copied from freesurfer: 08/15/2017
% Prev Mod: 09/14/2017 by Don Hagler
% Prev Mod: 03/02/2018 by Feng Xue
% Prev Mod: 03/02/2018 by Don Hagler
% Prev Mod: 06/24/2020 by Don Hagler
% Last Mod: 06/15/2026 by Don Hagler
%

vol = []; M = []; volsz = [];

if ~mmil_check_nargs(nargin,1), return; end
if ~exist('headeronly','var') || isempty(headeronly), headeronly = 0; end
%  if ~exist('keepsingle','var') || isempty(keepsingle), keepsingle = 0; end

% unzip if it is compressed 
ext = fname((strlen(fname)-2):strlen(fname));
if(strcmpi(ext,'.gz'))
  tempfname = mmil_tempfname;
  new_fname = sprintf('%s.nii', tempfname);
  %fprintf('Uncompressing %s to %s\n',fname,new_fname);
  if(strcmp(computer,'MAC') || strcmp(computer,'MACI') || ismac)
    mmil_unix(sprintf('gunzip -c %s > %s', fname, new_fname));
  else
    mmil_unix(sprintf('zcat %s > %s', fname, new_fname)) ;
  end
  fname = new_fname ;
  gzipped = 1 ;
else
  gzipped = 0 ;
end

hdr = fs_load_nifti_hdr(fname);

if(isempty(hdr)) 
  if(gzipped), delete(fname); end
  error('Error reading nifti file');
end

M = hdr.sform;

% Check for ico7
nspatial = prod(hdr.dim(2:4));
IsIco7 = 0;
if(nspatial == 163842) IsIco7 = 1; end

if(IsIco7)
  %fprintf('load_nifti: ico7 reshaping\n');
  hdr.dim(2) = 163842;
  hdr.dim(3) = 1;
  hdr.dim(4) = 1;
end
volsz = [hdr.dim(2) hdr.dim(3) hdr.dim(4) hdr.dim(5)];

% Convert from 0-based to 1-based indexing
M = M*[1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1];

% If only header is desired, return now
if(headeronly) 
  if(gzipped), delete(fname); end
  return; 
end

% Get total number of voxels
dim = hdr.dim(2:end);
ind0 = find(dim==0);
dim(ind0) = 1;
nvoxels = prod(dim);

% Open to read the pixel data
fp = fopen(fname,'r',hdr.endian);

% Get past the header
fseek(fp,round(hdr.vox_offset),'bof');

switch(hdr.datatype)
 % Note: 'char' seems to work upto matlab 7.1, but 'uchar' needed
 % for 7.2 and higher. 
 case   2, [hdr.vol nitemsread] = fread(fp,inf,'uchar');
 case   4, [hdr.vol nitemsread] = fread(fp,inf,'short');
 case   8, [hdr.vol nitemsread] = fread(fp,inf,'int');
 case  16, [hdr.vol nitemsread] = fread(fp,inf,'float');
 case  64, [hdr.vol nitemsread] = fread(fp,inf,'double');
 case 512, [hdr.vol nitemsread] = fread(fp,inf,'ushort');
 case 768, [hdr.vol nitemsread] = fread(fp,inf,'uint');
 otherwise,
   if(gzipped), delete(fname); end
   error('ERROR: data type %d not supported',hdr.datatype);
end

fclose(fp);
if(gzipped), delete(fname); end

% Check that that many voxels were read in
if(nitemsread ~= nvoxels) 
  error('ERROR: %s, read in %d voxels, expected %d\n',...
	  fname,nitemsread,nvoxels);
end

hdr.vol = reshape(hdr.vol, dim');
vol = hdr.vol;
scl_slope = hdr.scl_slope;
scl_inter = hdr.scl_inter;
clear hdr;
if(scl_slope ~= 0)
%  fprintf('%s: rescaling NIFTI: slope = %g, intercept = %g\n',...
%	  mfilename,scl_slope,scl_inter);
  %hdr.vol = hdr.vol * hdr.scl_slope  + hdr.scl_inter;
  vol = vol * scl_slope + scl_inter;
end

return;

function len = strlen(str)
%  len = strlen(str)
% compute the # of characters in str (ignoring 0s at the end)


%
% strlen.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

len = length(str) ;
for i=length(str):-1:1
        if (str(i) ~= 0)
                break ;
        end

        len = len-1;
end
return;
