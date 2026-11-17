function complete_flag = epi_topup_check_B0dx(fname_B0dx,topup_flag)
%function complete_flag = epi_topup_check_B0dx(fname_B0dx,[topup_flag])
%
% Purpose: check that a topup B0dx output set is complete, not just the
%   displacement volume
%
% Required:
%   fname_B0dx: full path name of B0 displacement field volume
%
% Optional:
%   topup_flag: [0|1] whether the field was estimated with FSL's topup
%     if 0, only the volume is required, since epi_estimateB0 writes no
%     fieldcoef or movpar file
%     {default: 1}
%
% Output:
%   complete_flag: [0|1] whether the output set is complete
%
% NOTE: epi_topup_estimateB0 writes three files per pair, the displacement
%   volume plus {fstem}_fieldcoef.nii and {fstem}_movpar.txt, and
%   epi_topup_unwarpB0 requires all three. An interrupted copy out of its tmp
%   dir can leave the volume without the other two, and callers that gate on
%   the volume alone then take a partial set for finished work: the estimate is
%   skipped and the unwarp errors on the missing movpar file, on every rerun
%
% Created:  08/20/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('topup_flag','var') || isempty(topup_flag), topup_flag=1; end;

complete_flag = 0;
if isempty(fname_B0dx) || ~exist(fname_B0dx,'file'), return; end;
complete_flag = 1;
if ~topup_flag, return; end;

[tpath,tstem] = fileparts(fname_B0dx);
fstem_dx = sprintf('%s/%s',tpath,tstem);
if ~exist(sprintf('%s_fieldcoef.nii',fstem_dx),'file') ||...
   ~exist(sprintf('%s_movpar.txt',fstem_dx),'file')
  complete_flag = 0;
end;

return;
