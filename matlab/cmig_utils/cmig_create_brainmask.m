function fname_mask = cmig_create_brainmask(fname,varargin)
%function fname_mask = cmig_create_brainmask(fname,[options])
%
% Required Input:
%   fname: full path of input file (mgz or nii)
%
% Optional Input ('key',value):
%   'synthstrip_flag': [0|1] use SynthStrip for brain masking; otherwise use FSL BET
%     {default = 1}
%   'warning_flag': [0|1] if 1, warn and return [] on failure; if 0, throw error
%     library errors always rethrow regardless of warning_flag
%     {default = 0}
%   'outdir': output directory
%     {default = []}
%   'fname_mask': output file name
%     {default = []}
%   'binarize_flag': [0|1] output binary mask
%     if 0, output the skull-stripped image with its original intensities
%     {default = 1}
%   'verbose': [0|1] display status messages
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  05/19/2026 by Don Hagler
% Prev Mod: 05/19/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Last Mod: 08/05/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_mask = [];
if ~mmil_check_nargs(nargin,1), return; end

parms = mmil_args2parms(varargin,{...
  'synthstrip_flag',true,[],...
  'warning_flag',false,[],...
  'outdir',[],[],...
  'fname_mask',[],[],...
  'binarize_flag',true,[],...
  'verbose',false,[],...
  'forceflag',false,[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.synthstrip_flag
  try
    fname_mask = mmil_synthstrip_brainmask(fname,...
      'outdir',parms.outdir,'fname_mask',parms.fname_mask,...
      'binarize_flag',parms.binarize_flag,...
      'verbose',parms.verbose,'forceflag',parms.forceflag);
  catch me
    if ~isempty(strfind(me.message,'error while loading shared libraries'))
      rethrow(me);
    end
    if parms.warning_flag
      fprintf('%s: WARNING: %s\n',mfilename,me.message);
      return;
    else
      rethrow(me);
    end
  end
else
  try
    fname_mask = mmil_fsl_brainmask(fname,...
      'outdir',parms.outdir,'fname_mask',parms.fname_mask,...
      'binarize_flag',parms.binarize_flag,...
      'verbose',parms.verbose,'forceflag',parms.forceflag);
  catch me
    if parms.warning_flag
      fprintf('%s: WARNING: %s\n',mfilename,me.message);
      return;
    else
      rethrow(me);
    end
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
