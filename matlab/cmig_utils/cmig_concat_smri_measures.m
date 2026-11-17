function [outdir,errcode,errmsg] = cmig_concat_smri_measures(participant_id,session_id,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_concat_smri_measures(participant_id,session_id,output_rootdir,[options])
%
% Purpose: Save per-subject sphere-resampled sMRI surface maps in BIDS .mat format
%   Reads sphere-smoothed MGZ files produced by cmig_analyze_bids_anat.
%   Outputs one .mat per measure x smoothing x hemisphere.
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
%   'measlist': cell array of sMRI measure names
%     {default = {'thk','area','sulc','t1_gwc','t2_gwc'}}
%   'smoothlist': sphere smoothing steps (scalar or vector)
%     {default = [0,16,256,1000]}
%   'hemilist': cell array of hemisphere labels
%     {default = {'lh','rh'}}
%   'ico_order': icosahedral order of the input sphere maps
%     {default = 5}
%   'forceflag': [0|1] overwrite existing output files {default = 0}
%
% Output:
%   outdir: output directory
%   errcode: 0 if successful, 1 if error
%   errmsg: error message string if errcode=1, otherwise []
%
% Input MGZ naming:
%   {tmpdir}/{fullmeas}-sphere-{hemi}.mgz          (smoothsteps=0)
%   {tmpdir}/{fullmeas}-sphere-sm{N}-{hemi}.mgz   (smoothsteps>0)
%
% Output BIDS .mat naming:
%   {outdir}/{pid}_{sid}_hemi-{L|R}_atlas-ico{N}_desc-fwhm{N}_{meas}.mat
%   e.g., hemi-L_atlas-ico5_desc-fwhm0_thk.mat
%
% NOTE: measlist is entirely surface-derived, so callers must skip this for a
%   volume-only FreeSurfer recon; a missing input file is fatal here
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
  'measlist',{'thk','area','sulc','t1_gwc','t2_gwc'},[],...
  'smoothlist',[0,16,256,1000],[],...
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

% save concatenated vertex-wise sMRI measures in .mat format
fprintf('%s: saving sMRI measures in concatenated format...\n',mfilename);
for measi = 1:length(parms.measlist)
  % set full measure name (used in input file names)
  fullmeas = cmig_get_smri_fullmeas(parms.measlist{measi});
  for smfi = 1:length(parms.smoothlist)
    smoothsteps = parms.smoothlist(smfi);
    smoothstrin = sprintf('sm%d',smoothsteps);
    for hemii = 1:length(parms.hemilist)
      hemi = parms.hemilist{hemii};
      measmat = [];
      bids_hemi = cmig_get_bids_hemi(hemi);
      bids_meas = strrep(parms.measlist{measi},'_','');
      bids_fwhm = round(cmig_smoothsteps_to_fwhm(smoothsteps));
      fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico%d_desc-fwhm%d_%s.mat',...
        outdir,parms.participant_id,parms.session_id,...
        bids_hemi,parms.ico_order,bids_fwhm,bids_meas);
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
        % set input file name
        if smoothsteps == 0
          fname = sprintf('%s/%s-sphere-%s.mgz',tmpdir,fullmeas,hemi);
        else
          fname = sprintf('%s/%s-sphere-%s-%s.mgz',tmpdir,fullmeas,smoothstrin,hemi);
        end
        if ~exist(fname,'file')
          errcode = 1;
          errmsg = sprintf('missing file %s',fname);
          fprintf('%s: ERROR: %s\n',mfilename,errmsg);
          return;
        end
        % load mgz file and catch potential file corruption
        try
          vals = squeeze(fs_load_mgh(fname));
          tmp = zeros(size(vals));
          tmp(verts,:) = vals(verts,:);
          vals = tmp(1:nverts,:);
        catch
          errcode = 1;
          errmsg = sprintf('corrupted file %s',fname);
          fprintf('%s: ERROR: %s\n',mfilename,errmsg);
          return;
        end
        if isempty(measmat)
          measmat = NaN([1 size(vals)],'single');
        end
        measmat(1,:,:) = vals;
        save(fname_out,'measmat','-v7');
        fprintf('%s: wrote file %s\n',mfilename,fname_out);
      end
    end
  end
end

return;
