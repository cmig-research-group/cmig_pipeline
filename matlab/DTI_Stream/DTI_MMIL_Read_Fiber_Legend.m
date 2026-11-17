function FiberLegend = DTI_MMIL_Read_Fiber_Legend(fname)
%function FiberLegend = DTI_MMIL_Read_Fiber_Legend([fname])
%
% Created:  03/22/10 by Don Hagler
% Prev Mod: 01/14/13 by Don Hagler
% Last Mod: 10/26/2026 by Don Hagler
%

if ~exist('fname','var') || isempty(fname)
  fname = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/DTI_Fiber/DTI_Fiber_Legend.csv'];
end;

FiberLegend = MMIL_Read_StudyInfo(fname);

