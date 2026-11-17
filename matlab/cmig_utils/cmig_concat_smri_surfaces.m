function [outdir,errcode,errmsg] = cmig_concat_smri_surfaces(participant_id,session_id,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_concat_smri_surfaces(participant_id,session_id,output_rootdir,[options])
%
% Purpose: Save per-subject sphere-resampled sMRI surface coordinates in BIDS
%   .mat format
%   Reads the resampled surface files produced by cmig_analyze_bids_anat.
%   Outputs one .mat per surface x hemisphere.
%
% Required Input:
%   participant_id: participant identifier string (e.g. 'sub-XXXXXXXX')
%   session_id: session identifier string (e.g. 'ses-00A')
%   output_rootdir: root directory for analysis outputs
%
% Optional Parameters ('key',value pairs):
%   'fname_info': full path of the vol_info mat file written by
%     cmig_save_vol_info, supplying the cortex vertex masks
%     {default = []}
%   'surflist': cell array of surface names
%     {default = {'pial','white'}}
%   'hemilist': cell array of hemisphere labels
%     {default = {'lh','rh'}}
%   'ico_order': icosahedral order of the input surfaces
%     {default = 5}
%   'forceflag': [0|1] overwrite existing output files {default = 0}
%
% Output:
%   outdir: output directory
%   errcode: 0 if successful, 1 if error
%   errmsg: error message string if errcode=1, otherwise []
%
% Input surface naming:
%   {tmpdir}/{hemi}.{surf}
%
% Output BIDS .mat naming:
%   {outdir}/{pid}_{sid}_hemi-{L|R}_atlas-ico{N}_{surf}.mat
%   e.g., hemi-L_atlas-ico5_white.mat
%
% NOTE: coordmat is [1 x nverts x 3], unlike the [1 x nverts] measmat written
%   by cmig_concat_smri_measures
%
% Created:  08/28/2026 by Don Hagler
% Prev Mod: 08/28/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,3), errcode = 1; return; end

parms = mmil_args2parms(varargin,{...
  'participant_id',participant_id,[],...
  'session_id',session_id,[],...
  'output_rootdir',output_rootdir,[],...
  'fname_info',[],[],...
  'surflist',{'pial','white'},{'pial','white'},...
  'hemilist',{'lh','rh'},[],...
  'ico_order',5,[0 7],...
  'forceflag',false,[false true],...
...
  'modality','anat',[],...
});

outdir = sprintf('%s/%s/%s/%s',...
  parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);
tmpdir = sprintf('%s/analyze_anat',outdir);
mmil_mkdir(outdir);

% save concatenated fsico vertex coordinates in .mat format
fprintf('%s: saving sMRI surfaces in concatenated format...\n',mfilename);
for hemii = 1:length(parms.hemilist)
  hemi = parms.hemilist{hemii};
  for surfi = 1:length(parms.surflist)
    coordmat = [];
    bids_hemi_surf = cmig_get_bids_hemi(hemi);
    fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico%d_%s.mat',...
      outdir,parms.participant_id,parms.session_id,bids_hemi_surf,parms.ico_order,parms.surflist{surfi});
    if ~exist(fname_out,'file') || parms.forceflag
      % load vol_info
      if ~exist('vol_info','var') || isempty(vol_info), vol_info = load(parms.fname_info); end
      switch hemi
        case 'lh'
          verts = find(vol_info.vol_mask_lh);
          nverts = length(vol_info.vol_mask_lh);
        case 'rh'
          verts = find(vol_info.vol_mask_rh);
          nverts = length(vol_info.vol_mask_rh);
      end
      % load resampled surface file
      fname = sprintf('%s/%s.%s',tmpdir,hemi,parms.surflist{surfi});
      if ~exist(fname,'file')
        errcode = 1;
        errmsg = sprintf('missing file %s',fname);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
      % load surface file and catch potential file corruption
      try
        vals = fs_read_surf(fname);
        vals = vals.vertices(1:nverts,:);
      catch
        errcode = 1;
        errmsg = sprintf('corrupted file %s',fname);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        return;
      end
      if isempty(coordmat)
        coordmat = NaN([1 size(vals)],'single');
      end
      coordmat(1,:,:) = vals;
      save(fname_out,'coordmat','-v7');
      fprintf('%s: wrote file %s\n',mfilename,fname_out);
    end
  end
end

return;
