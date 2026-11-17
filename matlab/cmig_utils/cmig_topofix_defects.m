function [defect_sizes,errcode,errmsg] = cmig_topofix_defects(fspath)
%function [defect_sizes,errcode,errmsg] = cmig_topofix_defects(fspath)
%
% Purpose: number of vertices in each defect corrected by FreeSurfer's
%   topology correction, from the ?h.defect_labels recon-all writes
%
% Required parameters:
%   fspath: full path of FreeSurfer output directory
%
% Output:
%   defect_sizes: vertices of ?h.orig.nofix in each corrected defect, lh
%     then rh; its length is the number of corrected defects
%   errcode: 0 if successful, 1 if error
%   errmsg: error message string if errcode=1, otherwise []
%
% NOTE: ?h.defect_labels holds, for each vertex of ?h.orig.nofix, the number
%   of the defect it belongs to (1 to N), 0 outside any defect; the counts
%   equal those mris_fix_topology reports. A corrected defect is not the same
%   as a surface hole (mri_segstats --euler), which counts only defects that
%   change the Euler number
%
% Created:  10/12/2026 by Don Hagler
% Last Mod: 10/12/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defect_sizes = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,1), errcode = 1; return; end

hemilist = {'lh','rh'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h = 1:length(hemilist)
  fname = sprintf('%s/surf/%s.defect_labels',fspath,hemilist{h});
  if ~exist(fname,'file')
    errcode = 1;
    errmsg = sprintf('file %s not found',fname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  labels = round(fs_read_curv(fname));
  ids = unique(labels(labels > 0));
  sizes = zeros(length(ids),1);
  for d = 1:length(ids)
    sizes(d) = sum(labels == ids(d));
  end
  defect_sizes = cat(1,defect_sizes,sizes);
end

return;
