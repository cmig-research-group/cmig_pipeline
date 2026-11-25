function [outdir,errcode,errmsg] = cmig_concat_dwi_measures(participant_id,session_id,output_rootdir,model,varargin)
%function [outdir,errcode,errmsg] = cmig_concat_dwi_measures(participant_id,session_id,output_rootdir,model,[options])
%
% Purpose: Save per-subject sphere-resampled dMRI surface maps in BIDS .mat format
%   Reads sphere-smoothed MGZ files produced by cmig_analyze_bids_dwi (analyze_gwcsurf).
%   Outputs one .mat per measure x layer x smoothing x hemisphere.
%
% Required Input:
%   participant_id: participant identifier string (e.g. 'sub-XXXXXXXX')
%   session_id: session identifier string (e.g. 'ses-00A')
%   output_rootdir: root directory for analysis outputs
%   model: model type: 'dti', 'rsi', or 'dti_inner'
%
% Optional Parameters ('key',value pairs):
%   'dt_measlist_full': cell array of DTI measure names
%     {default = {'FA','MD','LD','TD'}}
%   'rsi_measlist': cell array of RSI measure names
%     {default = {'N0','ND','NT','NIf','NTs2'}}
%   'gwcsurf_outfix': gwcsurf output file stem suffix (dti and rsi)
%     {default = 'gwcsurf'}
%   'gwcsurf_outfix_inner': gwcsurf output file stem suffix for dti_inner
%     {default = 'gwcsurf_inner'}
%   'ico_order': output icosahedral order; the input sphere maps are ico7
%     and are reduced to this order (nested FreeSurfer ico ordering)
%     {default = 5}
%   'sphsmoothsteps': sphere smoothing steps (scalar or vector)
%     {default = [0,16,256,1000]}
%   'layerlist': cell array of surface layer labels
%     {default = {'gm','wm','gwc'}}
%   'hemilist': cell array of hemisphere labels
%     {default = {'lh','rh'}}
%   'fnames_coverage_label': cell array of ico7 coverage label files, one per
%     hemisphere in hemilist order; values outside the label are set to NaN
%     (if empty, non-cortical vertices are set to 0)
%     {default = []}
%   'forceflag': [0|1] overwrite existing output files {default = 0}
%
% Output:
%   outdir: output directory
%   errcode: 0 if successful, 1 if error
%   errmsg: error message string if errcode=1, or a warning, otherwise []
%
% Input MGZ naming (from mmil_gwcsurf_analysis; ico7 sphere maps):
%   {tmpdir}/{meas}_{gwcsurf_outfix}_{layer}-sphere-{hemi}.mgz          (sphsmooth=0)
%   {tmpdir}/{meas}_{gwcsurf_outfix}_{layer}-sphere-sm{N}-{hemi}.mgz   (sphsmooth>0)
%
% Output BIDS .mat naming:
%   {outdir}/{pid}_{sid}_hemi-{L|R}_atlas-ico{N}_desc-fwhm{N}_{suffix}.mat
%   suffix: dti{meas}{layer}, rsi{extmeas}{layer}, or isdti{meas}{layer}
%   e.g., hemi-L_atlas-ico5_desc-fwhm0_dtiFAgm.mat
%
% Created:  06/09/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 06/25/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/28/2026 by Don Hagler
% Prev Mod: 10/09/2026 by Don Hagler
% Prev Mod: 10/22/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,4), errcode = 1; return; end

parms = mmil_args2parms(varargin,{...
  'participant_id',participant_id,[],...
  'session_id',session_id,[],...
  'output_rootdir',output_rootdir,[],...
  'model',model,[],...
...
  'dt_measlist_full',{'FA','MD','LD','TD'},[],...
  'rsi_measlist',{'N0','ND','NT','NIf','NTs2'},[],...
  'gwcsurf_outfix','gwcsurf',[],...
  'gwcsurf_outfix_inner','gwcsurf_inner',[],...
  'ico_order',5,[0 7],...
  'sphsmoothsteps',[0,16,256,1000],[0,10000],...
  'layerlist',{'gm','wm','gwc'},[],...
  'hemilist',{'lh','rh'},[],...
  'fnames_coverage_label',[],[],...
  'forceflag',false,[false true],...
...
  'modality','dwi',[],...
});

outdir = sprintf('%s/%s/%s/%s',...
  parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);
tmpdir = sprintf('%s/analyze_dwi',outdir);
mmil_mkdir(outdir);

% RSI internal → ABCD 7.0 external measure name mapping
% the internal to external RSI measure names, as get_rsi_extmeas in
%   cmig_analyze_bids_dwi names its other outputs
rsi_int = {'N0','ND','NT','NIf','NTs2','N0s2','NDs2',...
           'F0','FD','FT','If','FTs2','F0s2','FDs2'};
rsi_ext = {'RNI','RND','RNT','FNI','HNT','HNI','HND',...
           'RI','RD','RT','FI','HT','HI','HD'};

% number of vertices per ico order; input sphere maps are ico7 and are
% reduced to ico_order using nested FreeSurfer ico vertex ordering
ico_nverts = [42 162 642 2562 10242 40962 163842];
if parms.ico_order>=1 && parms.ico_order<=7
  nverts = ico_nverts(parms.ico_order);
else
  nverts = [];
end

% load fsaverage cortex labels to zero non-cortical (medial wall) vertices
verts_lh = []; verts_rh = [];
labelpath = sprintf('%s/subjects/fsaverage/label',getenv('FREESURFER_HOME'));
if exist(labelpath,'dir')
  verts_lh = fs_read_label(sprintf('%s/lh.cortex.label',labelpath));
  verts_rh = fs_read_label(sprintf('%s/rh.cortex.label',labelpath));
else
  errmsg = sprintf('label path %s not found; skipping cortex masking',labelpath);
  fprintf('%s: WARNING: %s\n',mfilename,errmsg);
end

% dispatch on model
rsi = strcmp(parms.model,'rsi');
dti_inner = strcmp(parms.model,'dti_inner');
if rsi
  meas_list = parms.rsi_measlist;
  gwcsurf_outfix = parms.gwcsurf_outfix;
  bids_prefix = 'rsi';
elseif dti_inner
  meas_list = parms.dt_measlist_full;
  gwcsurf_outfix = parms.gwcsurf_outfix_inner;
  bids_prefix = 'isdti';
else
  meas_list = parms.dt_measlist_full;
  gwcsurf_outfix = parms.gwcsurf_outfix;
  bids_prefix = 'dti';
end

fprintf('%s: saving sphere-sampled dMRI surface maps (%s)...\n',mfilename,parms.model);

for m = 1:length(meas_list)
  meas = meas_list{m};
  if rsi
    idx_ext = find(strcmp(rsi_int,meas));
    if ~isempty(idx_ext)
      bids_meas = rsi_ext{idx_ext};
    else
      bids_meas = meas;
    end
  else
    bids_meas = meas;
  end
  gwc_outstem = sprintf('%s/%s_%s',tmpdir,meas,gwcsurf_outfix);
  for l = 1:length(parms.layerlist)
    layer = parms.layerlist{l};
    for smfi = 1:length(parms.sphsmoothsteps)
      sphsm = parms.sphsmoothsteps(smfi);
      bids_fwhm = round(1.25 * sqrt(sphsm));
      for h = 1:length(parms.hemilist)
        hemi = parms.hemilist{h};
        bids_hemi = upper(hemi(1));
        bids_suffix = sprintf('%s%s%s',bids_prefix,bids_meas,layer);
        fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico%d_desc-fwhm%d_%s.mat',...
          outdir,parms.participant_id,parms.session_id,...
          bids_hemi,parms.ico_order,bids_fwhm,bids_suffix);
        if exist(fname_out,'file') && ~parms.forceflag, continue; end
        if sphsm == 0
          fname = sprintf('%s_%s-sphere-%s.mgz',gwc_outstem,layer,hemi);
        else
          fname = sprintf('%s_%s-sphere-sm%d-%s.mgz',gwc_outstem,layer,sphsm,hemi);
        end
        if ~exist(fname,'file')
          errcode = 1;
          errmsg = cmig_append_errmsg(errmsg,sprintf('missing sphere file %s',fname));
          fprintf('%s: ERROR: missing sphere file %s\n',mfilename,fname);
          return;
        end
        try
          vals = squeeze(fs_load_mgh(fname));
        catch err
          errcode = 1;
          tmp_errmsg = sprintf('failed to load %s: %s',fname,err.message);
          errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
          fprintf('%s: ERROR: %s\n',mfilename,tmp_errmsg);
          return;
        end
        try
          vals = vals(:);
          % set non-cortical vertices to zero
          switch hemi
            case 'lh', verts = verts_lh;
            case 'rh', verts = verts_rh;
          end
          if ~isempty(verts)
            tmp = zeros(size(vals));
            keep = verts(verts <= numel(vals));
            tmp(keep) = vals(keep);
            vals = tmp;
          end
          % set vertices outside the coverage label to NaN
          if ~isempty(parms.fnames_coverage_label)
            cov = fs_read_label(parms.fnames_coverage_label{h});
            tmp = nan(size(vals));
            keep = cov(cov <= numel(vals));
            tmp(keep) = vals(keep);
            vals = tmp;
          end
          % reduce ico7 sphere map to ico_order verts (ABCD_concat style)
          if ~isempty(nverts)
            if numel(vals) >= nverts
              vals = vals(1:nverts);
            else
              fprintf('%s: WARNING: %s has %d verts, fewer than ico%d (%d); keeping all\n',...
                mfilename,fname,numel(vals),parms.ico_order,nverts);
            end
          end
          measmat = single(vals(:)');
          save(fname_out,'measmat','-v7');
          fprintf('%s: wrote %s\n',mfilename,fname_out);
        catch err
          errcode = 1;
          tmp_errmsg = sprintf('failed to save %s: %s',fname_out,err.message);
          errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
          fprintf('%s: ERROR: %s\n',mfilename,tmp_errmsg);
          return;
        end
      end
    end
  end
end

fprintf('%s: surface map concat done (%s)\n',mfilename,parms.model);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
