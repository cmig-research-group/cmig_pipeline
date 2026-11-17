function fname_info = cmig_save_vol_info(participant_id,session_id,varargin)
%function fname_info = cmig_save_vol_info(participant_id,session_id,[options])
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%
% Optional Input ('key',value):
%   'outdir': output directory
%     {default = '.'}
%   'outstem': output file stem (no path, no extension)
%     {default = []}
%   'ico_order': ico order for cortical vertex masks; 0 = no mask
%     {default = 0}
%   'fname_tag': tag appended to outstem in filename: {outstem}_{fname_tag}_volinfo.mat
%     if empty: filename is {outstem}_volinfo.mat
%     caller is responsible for full tag string (e.g. 'atlas-ico5')
%     {default = ''}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  05/19/2026 by Don Hagler
% Prev Mod: 05/29/2026 by Don Hagler
% Prev Mod: 08/28/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_info = [];
if ~mmil_check_nargs(nargin,2), return; end

parms = mmil_args2parms(varargin,{...
  'outdir','.',[],...
  'outstem',[],[],...
  'ico_order',0,[],...
  'fname_tag','',[],...
  'forceflag',false,[],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(parms.fname_tag)
  fname_info = sprintf('%s/%s_volinfo.mat',parms.outdir,parms.outstem);
else
  fname_info = sprintf('%s/%s_%s_volinfo.mat',parms.outdir,parms.outstem,parms.fname_tag);
end

if exist(fname_info,'file') && ~parms.forceflag, return; end

dirlist        = {parms.outdir};
participant_id = {participant_id};
session_id     = {session_id};
visitidvec     = strcat(participant_id,'_',session_id);

if parms.ico_order > 0
  ico_nverts = [42 162 642 2562 10242 40962 163842];
  nverts = ico_nverts(parms.ico_order);
  labelpath = sprintf('%s/subjects/fsaverage/label',getenv('FREESURFER_HOME'));
  if exist(labelpath,'dir')
    verts_lh = fs_read_label(sprintf('%s/lh.cortex.label',labelpath));
    verts_rh = fs_read_label(sprintf('%s/rh.cortex.label',labelpath));
  else
    fprintf('%s: WARNING: label path %s not found\n',mfilename,labelpath);
    verts_lh = 1:nverts;
    verts_rh = 1:nverts;
  end
  vol_mask_lh = zeros(nverts,1); vol_mask_lh(intersect(verts_lh,1:nverts)) = 1;
  vol_mask_rh = zeros(nverts,1); vol_mask_rh(intersect(verts_rh,1:nverts)) = 1;
  save(fname_info,'participant_id','session_id','visitidvec','dirlist',...
    'vol_mask_lh','vol_mask_rh','-v7');
else
  vol_mask = [];
  save(fname_info,'participant_id','session_id','visitidvec','dirlist',...
    'vol_mask','-v7');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
