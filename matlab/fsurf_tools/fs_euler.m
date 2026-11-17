function num_defects = fs_euler(subj,subjdir,outdir,verbose)
%function num_defects = fs_euler(subj,[subjdir],[outdir],[verbose])
%
% Purpose: Calculate number of topological defects from Euler number
%
% Required input:
%  subj is a string specifying the subject name
%
% Optional parameters:
%  subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  outdir: output directory for euler.stats file
%    {default = <subjdir>/<subj>/stats}
%  verbose: [0|1] display status messages and warnings
%    {default = 0}
%
% Output:
%   num_defects: number of topological defects in ?h.orig.nofix surfaces
%     sum of left and right hemispheres
%
% Created:  01/27/20 by Don Hagler
% Prev Mod: 01/27/20 by Don Hagler
% Prev Mod: 06/24/20 by Don Hagler
% Prev Mod: 07/02/20 by Don Hagler
% Last Mod: 01/16/24 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_defects = [];
if ~mmil_check_nargs(nargin,1), return; end

if ~exist('subjdir','var') || isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end
end
if ~exist('outdir','var') || isempty(outdir)
  outdir = sprintf('%s/%s/stats',subjdir,subj);
end
if ~exist('verbose','var') || isempty(verbose)
  verbose = 0;
end

hemilist = {'lh','rh'};
nhemi = length(hemilist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(outdir);

fname_out = sprintf('%s/euler.stats',outdir);
cmd = sprintf('mri_segstats --sd %s --annot %s lh aparc --euler --sum %s',...
  subjdir,subj,fname_out);

if verbose
  fprintf('%s: running cmd:\n%s\n',mfilename,cmd);
end
[status,result,errmsg] = mmil_unix(cmd);
if status
  error('cmd %s failed:\n%s\n%s',cmd,result,errmsg);
end
if verbose
  fprintf('%s\n',result);
end

if ~exist(fname_out,'file'), error('file %s not found',fname_out); end

num_hemi_defects = nan(nhemi,1);

fid = fopen(fname_out);
for t=1:13, tmp = fgetl(fid); end
for h=1:nhemi
  hemi = hemilist{h};
  tmp = fgetl(fid);
  k = findstr(tmp,', ');
  if isempty(k)
    error('failed to get number of %s defect holes\n',hemi);
  else
    tmp = tmp(k(3)+2:k(4)-1);
    num_hemi_defects(h) = str2double(tmp);
  end
  if verbose
    fprintf('%s: %d %s defects\n',mfilename,num_hemi_defects(h),hemi);
  end
end
fclose(fid);

num_defects = sum(num_hemi_defects);
if verbose
  fprintf('%s: %d total defects\n',mfilename,num_defects);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

