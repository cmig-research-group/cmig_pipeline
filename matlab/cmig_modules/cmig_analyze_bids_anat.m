function [outdir,errcode,errmsg] = cmig_analyze_bids_anat(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin)
%function [outdir,errcode,errmsg] = cmig_analyze_bids_anat(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,[options])
%
% Purpose: Analyze structural MRI data from BIDS processed derivatives
%
% Required Input:
%   participant_id: participant identifier string
%     e.g. 'sub-XXXXXXXX'
%   session_id: session identifier string
%     e.g. 'ses-00A'
%   proc_rootdir: BIDS processed derivatives root directory
%   fsurf_rootdir: BIDS freesurfer derivatives root directory
%   output_rootdir: BIDS output root directory
%
% Optional Input ('key',value):
%   'cleanup_flag': [0|1] whether to remove temporary directories
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   outdir: output directory
%     will have form of {output_rootdir}/{participant_id}/{session_id}
%   errcode: error code (0 = success, 1 = error)
%   errmsg: error message string (empty if no error)
%
% Created:  10/24/2025 by Don Hagler
% Prev Mod: 11/18/2025 by Don Hagler
% Prev Mod: 04/17/2026 by Don Hagler
% Prev Mod: 04/24/2026 by Don Hagler
% Prev Mod: 04/29/2026 by Don Hagler
% Prev Mod: 05/26/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 06/17/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 06/24/2026 by Don Hagler
% Prev Mod: 07/10/2026 by Don Hagler
% Prev Mod: 07/13/2026 by Don Hagler
% Prev Mod: 08/11/2026 by Don Hagler
% Prev Mod: 08/13/2026 by Don Hagler
% Prev Mod: 08/28/2026 by Don Hagler
% Prev Mod: 09/09/2026 by Don Hagler
% Prev Mod: 09/28/2026 by Don Hagler
% Prev Mod: 10/21/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Prev Mod: 11/11/2026 by Don Hagler
% Prev Mod: 11/12/2026 by Don Hagler
% Last Mod: 11/20/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on cmig_analyze_bids_rsfmri
%   Created:  09/04/2025 by Don Hagler
%   Last Mod: 10/22/2025 by Don Hagler
% Based on cmig_proc_bids_anat
%   Created:  08/04/2025 by Don Hagler
%   Last Mod: 08/14/2025 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outdir = []; errcode = 0; errmsg = [];
if ~mmil_check_nargs(nargin,5), return; end
% stamped at entry so the status row carries a duration of its own, rather
%   than one inferred by differencing against the module before it
startdate = datestr(now,'yyyy-mm-ddTHH:MM:SS.FFF');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parameters
parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,varargin);

% locate input data
[parms,errcode,errmsg] = check_input(parms);
if errcode
  % write status file
  cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate);
  return;
end

% check whether processing required
[parms,procflag] = check_output(parms);

% process data
if procflag

  % import T1w and T2w (convert to mgz in tmpdir)
  parms.fname_t1w = cmig_import_to_mgz(parms.fname_t1w,parms.tmpdir,'T1w',parms.forceflag);
  if ~isempty(parms.fname_t2w)
    parms.fname_t2w = cmig_import_to_mgz(parms.fname_t2w,parms.tmpdir,'T2w',parms.forceflag);
  end

  % surface smoothing and ROI averages for sMRI measures
  [errcode,tmp_errmsg] = analyze_smri(parms);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end

  % resample surfaces to ico
  errcode_surf = cmig_resamp_surfs(parms.fspath,parms.tmpdir,parms.forceflag);
  if errcode_surf
    tmp_errmsg = 'failed to resample surfaces';
    fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  end

  % save vol_info
  bids_stem = sprintf('%s_%s',parms.participant_id,parms.session_id);
  parms.fname_info = cmig_save_vol_info(parms.participant_id,parms.session_id,...
    'outdir',parms.outdir,'outstem',bids_stem,...
    'ico_order',parms.ico_order,...
    'fname_tag',sprintf('atlas-ico%d',parms.ico_order),'forceflag',parms.forceflag);

  % save sMRI measures in GIfTI format
  [errcode,tmp_errmsg] = convert_smri_measures(parms);
  errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
  if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end

  % save sMRI measures in concatenated data format
  % NOTE: measlist is entirely surface-derived, so skip for a volume-only
  %   recon; the function treats a missing input file as fatal
  if parms.fs_complete_flag
    [~,errcode,tmp_errmsg] = cmig_concat_smri_measures(parms.participant_id,parms.session_id,...
      parms.output_rootdir,'fname_info',parms.fname_info,...
      'measlist',parms.measlist,'smoothlist',parms.smoothlist,...
      'hemilist',parms.hemilist,'ico_order',parms.ico_order,...
      'forceflag',parms.forceflag);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end
  else
    fprintf('%s: WARNING: skipping sMRI measure concatenation (volume-only FreeSurfer recon)\n',mfilename);
  end

  % save sMRI surfaces in concatenated data format
  if ~errcode_surf
    [~,errcode,tmp_errmsg] = cmig_concat_smri_surfaces(parms.participant_id,parms.session_id,...
      parms.output_rootdir,'fname_info',parms.fname_info,...
      'surflist',parms.surflist,'hemilist',parms.hemilist,...
      'ico_order',parms.ico_order,'forceflag',parms.forceflag);
    errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
    if errcode, cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate); return; end
  end

  % write tabulated data
  write_tabulated_data(parms);
end

% cleanup tmpdirs
if parms.cleanup_flag
  errmsg = cmig_append_errmsg(errmsg,cmig_cleanup_tmpdirs(parms.outdir,parms.tmpdir));
end

% write status file
cmig_write_status(parms.participant_id,parms.session_id,parms.outdir,'desc-analyze',errcode,errmsg,parms.status_extra,startdate,procflag);

outdir = parms.outdir;
fprintf('%s: finished\n',mfilename);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_parms(participant_id,session_id,proc_rootdir,fsurf_rootdir,output_rootdir,options)
  options = cmig_parse_deployed_args(options);
  parms = mmil_args2parms(options,{...
    'participant_id',participant_id,[],...
    'session_id',session_id,[],...
    'proc_rootdir',proc_rootdir,[],...
    'fsurf_rootdir',fsurf_rootdir,[],...
    'output_rootdir',output_rootdir,[],...
  ...
    'cleanup_flag',true,[false true],...
    'forceflag',false,[false true],...
  ...
    'modality','anat',[],...
    'hemilist',{'lh','rh'},[],...
    'measlist',{'thk','area','sulc','t1_gwc','t2_gwc',...
                't1_gm','t2_gm','t1_wm','t2_wm'},[],...
    'surflist',{'pial','white'},{'pial','white'},...
    'tablist',{'thk','area','vol','sulc','t1','t2'},[],...
    'atlaslist',{'dsk','dst','fzy','aseg'},[],...
    'smoothlist',[0,16,256,1000],[],...
    'ico_order',5,[0 7],...
  ...
    'aparc_infix',{'aparc','aparc.a2009s'},[],...
    'aseg_roilist',[2 3 4 5 7 8 10 11 12 13 14 15 16 17 18 24 25 26 28 41 42 43 44 46 47 49 50 51 52 53 54 57 58 60 77 78 79 251 252 253 254 255],[1,Inf],...
    'aseg_roigroups_flag',true,[false true],...
    'aseg_roigroups',[],[],...
    'fname_colorlut',[],[],...
    'fuzzy_dir',[],[],...
    'fuzzy_fstem','fuzzy',[],...
    'fuzzy_order',18,[0,2,4,12,18],...
    'fuzzy_smooth',2819,[0,Inf],...
    'fuzzy_thresh',0,[0,Inf],...
    'projdist_list',[-0.2,0.2],[-5,5],...
    'status_extra',[],[],...
  });

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,errcode,errmsg] = check_input(parms)
  errcode = 0; errmsg = [];

  % set output directory
  parms.outdir = sprintf('%s/%s/%s/%s',parms.output_rootdir,parms.participant_id,parms.session_id,parms.modality);
  
  % inititalize inputs
  parms.fname_t1w = []; parms.fname_t2w = []; parms.fspath = [];

  fprintf('%s: checking input for %s_%s in %s\n',mfilename,...
    parms.participant_id,parms.session_id,parms.proc_rootdir);
   
  % get input file names
  parms.fname_t1w = get_input_data(parms.participant_id,parms.session_id,parms.proc_rootdir,'T1w');
  parms.fname_t2w = get_input_data(parms.participant_id,parms.session_id,parms.proc_rootdir,'T2w');
  % proc file names, for the sidecar Sources
  parms.fname_t1w_proc = parms.fname_t1w;
  parms.fname_t2w_proc = parms.fname_t2w;

  % check that there are input files to analyze
  if isempty(parms.fname_t1w)
    errcode = 1;
    errmsg = 'not processed: no usable T1w scan';
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
  
  % warn if no T2w scan
  if isempty(parms.fname_t2w)
    fprintf('%s: WARNING: no usable T2w scan\n',mfilename);
    parms.measlist = parms.measlist(~contains(parms.measlist,'t2'));
    parms.tablist  = parms.tablist(~contains(parms.tablist,'t2'));
  end

  % report input files found
  fprintf('%s: input T1w data:\n %s\n',mfilename,parms.fname_t1w);
  if ~isempty(parms.fname_t2w)
    fprintf('%s: input T2w data:\n %s\n',mfilename,parms.fname_t2w);
  end
    
  [parms.fspath,parms.fs_version,parms.fs_build,parms.fs_complete_flag,errcode,errmsg] = ...
    cmig_check_bids_fsurf(parms.participant_id,parms.session_id,parms.fsurf_rootdir);
  parms.status_extra = struct('fs_version',parms.fs_version,...
    'fs_build',parms.fs_build,'fs_complete_flag',parms.fs_complete_flag);
  if errcode
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_data = get_input_data(participant_id,session_id,proc_rootdir,series_type)
  fname_data = [];
  [fnames_data,fnames_info,scan_info] = cmig_get_bids_data(participant_id,session_id,proc_rootdir,'anat',series_type);
  if isempty(fnames_data), return; end
  idx = cmig_select_scan(scan_info);
  fname_data = fnames_data{idx};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,procflag] = check_output(parms)
  procflag = 0;

  % set output file stem, and temporary directory
  parms.outstem = sprintf('%s_%s_smri',parms.participant_id,parms.session_id);
  parms.tmpdir = sprintf('%s/analyze_anat',parms.outdir);

  % check whether output directory exists
  if ~exist(parms.outdir,'dir')
    procflag = 1;
  elseif parms.forceflag
    procflag = 1;
  else
    % create list of expected output files
    % NOTE: surface-derived outputs are expected only for a complete recon,
    %   matching what analyze_smri and write_tabdat produce; without this
    %   gate a volume-only session reprocesses on every invocation
    output_files = {};
    for hemii = 1:length(parms.hemilist)
      if ~parms.fs_complete_flag, break; end
      hemi = parms.hemilist{hemii};
      % sMRI measures
      for measi = 1:length(parms.measlist)
        for smfi = 1:length(parms.smoothlist)
          smoothsteps = parms.smoothlist(smfi);
          % sMRI measures in concatenated format
          bids_hemi = cmig_get_bids_hemi(hemi);
          bids_meas = strrep(parms.measlist{measi},'_','');
          bids_fwhm = round(cmig_smoothsteps_to_fwhm(smoothsteps));
          fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico%d_desc-fwhm%d_%s.mat',...
            parms.outdir,parms.participant_id,parms.session_id,...
            bids_hemi,parms.ico_order,bids_fwhm,bids_meas);
          output_files{end+1} = fname_out;
          % sMRI measures converted to GIfTI
          fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico7_desc-fwhm%d_%s.%s.gii',...
            parms.outdir,parms.participant_id,parms.session_id,bids_hemi,bids_fwhm,bids_meas,...
            get_gii_type(parms.measlist{measi}));
          output_files{end+1} = fname_out;
          output_files{end+1} = json_name(fname_out);
        end
      end
      
      % sMRI surfaces
      for surfi = 1:length(parms.surflist)
        fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico%d_%s.mat',...
          parms.outdir,parms.participant_id,parms.session_id,bids_hemi,parms.ico_order,parms.surflist{surfi});
        output_files{end+1} = fname_out;
      end
    end    

    % tabulated data files
    bids_base = sprintf('%s_%s',parms.participant_id,parms.session_id);
    for m=1:length(parms.tablist)
      meas = parms.tablist{m};
      bids_meas_tab = strrep(meas,'_','');
      for a=1:length(parms.atlaslist)
        atlas = parms.atlaslist{a};
        if ~parms.fs_complete_flag && ~strcmp(atlas,'aseg'), continue; end
        if strcmp(atlas,'aseg') && ismember(meas,{'thk','area','sulc'}), continue; end
        gwc_flag = (ismember(atlas,{'dsk','dst','fzy'}) && ismember(meas,{'t1','t2'}));
        if gwc_flag && ~strcmp(atlas,'aseg')
          suffix_list = {'wm','gm','gwc'};
          for s=1:length(suffix_list)
            fname_out = sprintf('%s/%s_atlas-%s_%s%s.tsv',...
              parms.outdir,bids_base,atlas,bids_meas_tab,suffix_list{s});
            output_files{end+1} = fname_out;
          end
        else
          fname_out = sprintf('%s/%s_atlas-%s_%s.tsv',...
            parms.outdir,bids_base,atlas,bids_meas_tab);
          output_files{end+1} = fname_out;
        end
      end
    end

    % check whether output files exist
    for i=1:length(output_files)
      if ~exist(output_files{i},'file')
        [~,name,ext] = fileparts(output_files{i});
        missing_file = [name,ext];
        fprintf('%s: missing output file %s\n',mfilename,missing_file);      
        procflag = 1;
        break;
      end
    end
  end
    
  % create temporary directory (and output directory)
  if procflag, mmil_mkdir(parms.tmpdir); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = analyze_smri(parms)
  errcode = 0; errmsg = [];
  parms = setup_smri_parms(parms);
  [errcode,errmsg] = read_stats(parms);
  if errcode, return; end
  [errcode,errmsg] = create_aseg_roigroups(parms);
  if errcode, return; end
  [errcode,errmsg] = aseg_analysis(parms);
  if errcode, return; end
  if parms.fs_complete_flag
    [errcode,errmsg] = cortsurf_analysis(parms);
    if errcode, return; end
    [errcode,errmsg] = aparc_surf_analysis(parms,'sulc');
    if errcode, return; end
    [errcode,errmsg] = gwc_analysis(parms);
    if errcode, return; end
    [errcode,errmsg] = fuzzy_analysis(parms);
  else
    fprintf('%s: WARNING: skipping surface-dependent analyses (volume-only FreeSurfer recon)\n',mfilename);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = convert_smri_measures(parms)
  errcode = 0; errmsg = [];

  % parms.measlist is entirely surface-derived; nothing to convert for
  % volume-only recon
  if ~parms.fs_complete_flag
    fprintf('%s: WARNING: skipping sMRI measure GIfTI conversion (volume-only FreeSurfer recon)\n',mfilename);
    return;
  end

  % save vertex-wise sMRI measures in .gii format
  fprintf('%s: saving sMRI measures in GIfTI format...\n',mfilename);
  for measi = 1:length(parms.measlist)
    % set full measure name (used in input file names)
    fullmeas = cmig_get_smri_fullmeas(parms.measlist{measi});
    for smfi = 1:length(parms.smoothlist)
      smoothsteps = parms.smoothlist(smfi);
      smoothstrin = sprintf('sm%d',smoothsteps);
      for hemii = 1:length(parms.hemilist)
        hemi = parms.hemilist{hemii};
        bids_hemi = cmig_get_bids_hemi(hemi);
        bids_meas_gii = strrep(parms.measlist{measi},'_','');
        bids_fwhm_gii = round(cmig_smoothsteps_to_fwhm(smoothsteps));
        fname_out = sprintf('%s/%s_%s_hemi-%s_atlas-ico7_desc-fwhm%d_%s.%s.gii',...
          parms.outdir,parms.participant_id,parms.session_id,bids_hemi,bids_fwhm_gii,bids_meas_gii,...
          get_gii_type(parms.measlist{measi}));
        if ~exist(fname_out,'file') || parms.forceflag
          % set input file name          
          if smoothsteps == 0
            fname = sprintf('%s/%s-sphere-%s.mgz',parms.tmpdir,fullmeas,hemi);
          else
            fname = sprintf('%s/%s-sphere-%s-%s.mgz',parms.tmpdir,fullmeas,smoothstrin,hemi);
          end
          if ~exist(fname,'file')
            tmp_errmsg = sprintf('missing file %s',fname);
            fprintf('%s: WARNING: %s\n',mfilename,tmp_errmsg);
            errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
            continue;
          end
          % convert mgz to gii
          fs_mri_convert(fname,fname_out);
          [errcode,tmp_errmsg] = strip_gii_metadata(fname_out);
          errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
          if errcode, return; end
          fprintf('%s: wrote file %s\n',mfilename,fname_out);
        end
        % sidecar json (needs only the parms, so restored without the tmpdir)
        fname_json = json_name(fname_out);
        if ~exist(fname_json,'file') || parms.forceflag
          info = gii_sidecar_info(parms,parms.measlist{measi},hemi,bids_fwhm_gii);
          [errcode,tmp_errmsg] = write_json(info,fname_json);
          errmsg = cmig_append_errmsg(errmsg,tmp_errmsg);
          if errcode, return; end
        end
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function info = gii_sidecar_info(parms,meas,hemi,fwhm)
  % sidecar of an ico7 sMRI measure GIfTI
  % NOTE: fs_paint and mmil_cortsurf_analysis resample to the fsaverage
  %   sphere with all its vertices, so the maps are on ico7
  [ico_ref,ico_desc] = cmig_fsaverage_ref(7);
  pdist_wm = abs(min(parms.projdist_list));
  pdist_gm = abs(max(parms.projdist_list));
  sphere_uri = fsurf_uri(parms,sprintf('surf/%s.sphere.reg',hemi));
  switch meas
    case 'thk'
      desc = 'Cortical thickness (mm)';
      sources = {fsurf_uri(parms,sprintf('surf/%s.thickness',hemi)),sphere_uri};
    case 'area'
      desc = 'Cortical surface area of the white surface (mm^2 per vertex)';
      sources = {fsurf_uri(parms,sprintf('surf/%s.area',hemi)),sphere_uri};
    case 'sulc'
      desc = 'Sulcal depth';
      sources = {fsurf_uri(parms,sprintf('surf/%s.sulc',hemi)),sphere_uri};
    otherwise
      % t1_gm, t1_wm, t1_gwc, t2_gm, t2_wm, t2_gwc
      tok = regexp(meas,'^t(\d)_(\w+)$','tokens','once');
      if strcmp(tok{1},'1')
        scantype = 'T1w'; fname_proc = parms.fname_t1w_proc;
      else
        scantype = 'T2w'; fname_proc = parms.fname_t2w_proc;
      end
      switch tok{2}
        case 'gm'
          desc = sprintf(['%s intensity sampled %g mm outside the white surface '...
            '(gray matter)'],scantype,pdist_gm);
        case 'wm'
          desc = sprintf(['%s intensity sampled %g mm inside the white surface '...
            '(white matter)'],scantype,pdist_wm);
        otherwise
          desc = sprintf(['%s gray/white contrast, 2*(wm-gm)/(wm+gm), from the '...
            'intensities sampled %g mm inside and %g mm outside the white surface'],...
            scantype,pdist_wm,pdist_gm);
      end
      sources = {proc_uri(parms,fname_proc),...
                 fsurf_uri(parms,sprintf('surf/%s.white',hemi)),sphere_uri};
  end
  if fwhm==0
    fwhm_desc = ', unsmoothed';
  else
    fwhm_desc = sprintf(', smoothed on the sphere to FWHM %d mm',fwhm);
  end
  info = [];
  info.Description = sprintf('%s, resampled to %s%s',desc,ico_desc,fwhm_desc);
  info.Sources = sources;
  info.SpatialReference = ico_ref;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = proc_uri(parms,fname)
  % BIDS URI of a file in proc_rootdir, under the dataset name cmig_proc
  [~,fstem,fext] = cmig_get_fileparts(fname);
  uri = sprintf('bids:cmig_proc:%s/%s/%s/%s%s',...
    parms.participant_id,parms.session_id,parms.modality,fstem,fext);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uri = fsurf_uri(parms,relpath)
  % BIDS URI of a file in the recon, under the dataset name freesurfer
  uri = sprintf('bids:freesurfer:%s_%s/%s',...
    parms.participant_id,parms.session_id,relpath);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fname_json = json_name(fname)
  % the sidecar of a surface file: its whole extension replaced
  fname_json = regexprep(fname,'(\.(func|shape|label))?\.gii$','.json');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = write_json(info,fname)
  errcode = 0; errmsg = [];
  % NOTE: PrettyPrint was added in R2021a (MATLAB 9.10); guarded so an
  %   older MATLAB writes a single line rather than failing
  if verLessThan('matlab','9.10')
    json = jsonencode(info);
  else
    json = jsonencode(info,'PrettyPrint',true);
  end
  fid = fopen(fname,'wt');
  if fid<0
    errmsg = sprintf('failed to open %s for writing',fname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end
  fprintf(fid,'%s',json);
  fclose(fid);
  fprintf('%s: wrote file %s\n',mfilename,fname);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GIfTI intent by content: shape for morphometry, func for intensity samples
function gii_type = get_gii_type(meas)
  if ismember(meas,{'thk','area','sulc'})
    gii_type = 'shape';
  else
    gii_type = 'func';
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove the UserName and Date metadata that mri_convert writes into a GIfTI
function [errcode,errmsg] = strip_gii_metadata(fname)
  errcode = 0; errmsg = [];
  str = fileread(fname);
  str = regexprep(str,...
    '\s*<MD>\s*<Name><!\[CDATA\[(UserName|Date)\]\]></Name>\s*<Value>.*?</Value>\s*</MD>','');
  fid = fopen(fname,'w');
  if fid<0
    errmsg = sprintf('failed to open %s for writing',fname);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
    return;
  end
  fwrite(fid,str,'char');
  fclose(fid);
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tabulated_data(parms)
  colstem = 'mr_y_smri';
  if isempty(parms.outstem)
    outstem = 'smri';
  else
    outstem = parms.outstem;
  end
  bids_base = sprintf('%s_%s',parms.participant_id,parms.session_id);
  fprintf('%s: writing tabulated data...\n',mfilename);
  tic;
  for m=1:length(parms.tablist)
    meas = parms.tablist{m};
    bids_meas = strrep(meas,'_','');
    for a=1:length(parms.atlaslist)
      atlas = parms.atlaslist{a};
      % non-aseg atlases (dsk/dst/fzy) are surface parcellations;
      % skip them entirely for volume-only recon
      if ~parms.fs_complete_flag && ~strcmp(atlas,'aseg'), continue; end
      if strcmp(atlas,'aseg') && ismember(meas,{'thk','area','sulc'}), continue; end
      bids_fname = sprintf('%s/%s_atlas-%s_%s.tsv',...
        parms.outdir,bids_base,atlas,bids_meas);
      p = [];
      p.outdir = parms.outdir;
      p.outstem = outstem;
      p.colstem = colstem;
      p.meas = meas;
      p.atlas = atlas;
      p.fname_out = bids_fname;
      if ismember(atlas,{'dsk','dst','fzy'}) && ismember(meas,{'t1','t2'})
        p.gwc_flag = 1;
      end
      p.fs_version = parms.fs_version;
      p.tabdir = []; % if empty, get tabdir from $CMIG_PIPELINE_DIR/parms/tabulated
      p.fname_colorlut = parms.fname_colorlut; % if empty, CMIG_ColorLUT.txt
      p.forceflag = parms.forceflag;
      a = mmil_parms2args(p);
      cmig_write_tabdat_smri('indir',parms.tmpdir,...
        'participant_id',parms.participant_id,...
        'session_id',parms.session_id,...
        a{:});
    end
  end
  toc;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = setup_smri_parms(parms)
  % subjdir and subj from fspath
  [parms.subjdir,subj_base,subj_ext] = fileparts(parms.fspath);
  parms.subj = [subj_base subj_ext];
  % number of hemispheres
  parms.nhemi = length(parms.hemilist);
  % aseg file
  parms.fname_aseg = sprintf('%s/mri/aseg.mgz',parms.fspath);
  % proc file lists
  parms.proc_inputlist = {parms.fname_t1w};
  parms.proc_outputlist = {'T1w'};
  parms.proc_scalefacts = 1;
  if ~isempty(parms.fname_t2w)
    parms.proc_inputlist{2} = parms.fname_t2w;
    parms.proc_outputlist{2} = 'T2w';
    parms.proc_scalefacts(2) = 1;
  end
  % ensure aparc_infix is a cell array
  if ~iscell(parms.aparc_infix)
    parms.aparc_infix = {parms.aparc_infix};
  end
  % fs_paint output type and smoothing steps
  parms.outtype = 'mgz';
  parms.sphsmoothsteps = parms.smoothlist;
  % cortvol smoothing settings
  parms.icosubj = 'fsaverage';
  parms.icosubjdir = [getenv('FREESURFER_HOME') '/subjects'];
  % fuzzy cluster setup
  if isempty(parms.fuzzy_dir)
    parms.fuzzy_dir = [cmig_getenv('CMIG_PIPELINE_ATLASDIR') '/fuzzy_clusters'];
  end
  if parms.fuzzy_order == 0
    parms.fuzzy_outfix = parms.fuzzy_fstem;
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      fname_txt = sprintf('%s/%s_%s_roinames.txt',...
        parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
      if exist(fname_txt,'file')
        parms.fuzzy_names{h} = mmil_readtext(fname_txt);
      else
        parms.fuzzy_names{h} = [];
      end
    end
  else
    parms.fuzzy_outfix = sprintf('%s%d',parms.fuzzy_fstem,parms.fuzzy_order);
    for h=1:parms.nhemi
      parms.fuzzy_names{h} = mmil_fuzzy_names(...
        'fuzzy_fstem',parms.fuzzy_fstem,'fuzzy_order',parms.fuzzy_order);
    end
  end
  parms.fnames_weights = cell(parms.nhemi,1);
  for h=1:parms.nhemi
    fname_weights = sprintf('%s/%s-%s.mgz',...
      parms.fuzzy_dir,parms.fuzzy_outfix,parms.hemilist{h});
    if ~exist(fname_weights,'file')
      error('%s: fuzzy cluster file %s not found',mfilename,fname_weights);
    end
    parms.fnames_weights{h} = fname_weights;
  end
  parms.weights_thresh = parms.fuzzy_thresh;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = read_stats(parms)
  errcode = 0; errmsg = [];
  % read aseg stats, save to mat file
  fname_mat = sprintf('%s/aseg_stats.mat',parms.tmpdir);
  if ~exist(fname_mat,'file') || parms.forceflag
    aseg_roigroups = cmig_aseg_roigroups;
    aseg_stats = fs_read_aseg_stats(parms.subj,parms.subjdir,aseg_roigroups);
    if ~isempty(aseg_stats)
      save(fname_mat,'aseg_stats');
    end
  end
  % read aparc stats, save to mat file (one per aparc_infix)
  % (aparc stats require the surface pipeline; skip for volume-only recon)
  if parms.fs_complete_flag
    for i=1:length(parms.aparc_infix)
      aparc_infix = parms.aparc_infix{i};
      fname_mat = sprintf('%s/%s_stats.mat',parms.tmpdir,aparc_infix);
      if ~exist(fname_mat,'file') || parms.forceflag
        aparc_stats = fs_read_aparc_stats(parms.subj,parms.subjdir,aparc_infix);
        if ~isempty(aparc_stats)
          save(fname_mat,'aparc_stats');
        end
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = create_aseg_roigroups(parms)
  errcode = 0; errmsg = [];
  if ~parms.aseg_roigroups_flag, return; end

  aseg_roigroups = parms.aseg_roigroups;
  if isempty(aseg_roigroups)
    aseg_roigroups = cmig_aseg_roigroups;
  end

  try
    for j=1:length(aseg_roigroups)
      fname_mask = sprintf('%s/aseg_roigroup_%d.mgz',parms.tmpdir,aseg_roigroups(j).roicode);
      if ~exist(fname_mask,'file') || parms.forceflag
        [vol,M,mr_parms] = fs_load_mgh(parms.fname_aseg);
        ind = ismember(vol(:),aseg_roigroups(j).roicodes);
        vol(ind)  = aseg_roigroups(j).roicode;
        vol(~ind) = 0;
        fs_save_mgh(vol,fname_mask,M,mr_parms);
      end
    end
  catch me
    errmsg = sprintf('aseg roigroup mask creation failed: %s',me.message);
    fprintf('%s: ERROR: %s\n',mfilename,errmsg);
    errcode = 1;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roigroups = cmig_aseg_roigroups()
  % the summary ROIs of fs_define_aseg_roigroups at the codes of
  %   CMIG_ColorLUT.txt; fs_define_aseg_roigroups keeps MMIL_ColorLUT.txt's
  %   10001-10003, which FreeSurfer 8.x's FreeSurferColorLUT reuses
  roigroups(1).roiname  = 'WholeBrain';
  roigroups(1).roicode  = 21001;
  roigroups(1).roicodes = [2 3 7 8 10 11 12 13 17 18 26 28 41 42 46 47 ...
                           49 50 51 52 53 54 58 60 78 79];
  roigroups(2).roiname  = 'LatVentricles';
  roigroups(2).roicode  = 21002;
  roigroups(2).roicodes = [4 5 43 44];
  roigroups(3).roiname  = 'AllVentricles';
  roigroups(3).roicode  = 21003;
  roigroups(3).roicodes = [4 5 14 15 43 44 72 75 76 213];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = aseg_analysis(parms)
  errcode = 0; errmsg = [];
  for i=1:length(parms.proc_inputlist)
    fname_in = parms.proc_inputlist{i};
    outstem = sprintf('%s/%s',parms.tmpdir,parms.proc_outputlist{i});
    try
      mmil_aseg_analysis(fname_in,parms.fspath,...
        'outdir',parms.tmpdir,...
        'outstem',outstem,...
        'fname_aseg',parms.fname_aseg,...
        'aseg_roilist',parms.aseg_roilist,...
        'aseg_aparc_flag',0,...
        'erode_flag',0,...
        'scalefact',parms.proc_scalefacts(i),...
        'FS_version',parms.fs_version,...
        'verbose',1,...
        'forceflag',parms.forceflag);
    catch me
      errmsg = sprintf('aseg analysis failed for %s: %s',...
        parms.proc_outputlist{i},me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = cortsurf_analysis(parms)
  errcode = 0; errmsg = [];
  % smooth surface measures to ico sphere
  measlist = {'thickness','area','sulc'};
  for mi=1:length(measlist)
    meas = measlist{mi};
    outstem = sprintf('%s/%s',parms.tmpdir,meas);
    for si=1:length(parms.sphsmoothsteps)
      sphsmoothsteps = parms.sphsmoothsteps(si);
      try
        fs_paint(parms.subj,[],...
          'outstem',outstem,...
          'meas',meas,...
          'outtype',parms.outtype,...
          'subjdir',parms.subjdir,...
          'smoothsteps',0,...
          'sphere_flag',1,...
          'sphsmoothsteps',sphsmoothsteps,...
          'mask_midbrain_flag',0,...
          'hemilist',parms.hemilist,...
          'forceflag',parms.forceflag);
      catch me
        errmsg = sprintf('cortsurf analysis failed for %s: %s',meas,me.message);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        errcode = 1;
        return;
      end
    end
  end
  % compute and smooth cortical volume (needed for tabulated vol data)
  [errcode,errmsg] = smooth_cortvol(parms);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = smooth_cortvol(parms)
  errcode = 0; errmsg = [];
  % get paths to unsmoothed sphere files (already produced by cortsurf_analysis)
  fnames_thick = cell(parms.nhemi,1);
  fnames_area  = cell(parms.nhemi,1);
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fnames_thick{h} = sprintf('%s/thickness-sphere-%s.%s',parms.tmpdir,hemi,parms.outtype);
    fnames_area{h}  = sprintf('%s/area-sphere-%s.%s',parms.tmpdir,hemi,parms.outtype);
  end
  % compute cortvol = thickness * area
  outstem = sprintf('%s/cortvol-sphere',parms.tmpdir);
  for h=1:parms.nhemi
    hemi = parms.hemilist{h};
    fname_vol = sprintf('%s-%s.%s',outstem,hemi,parms.outtype);
    if ~exist(fname_vol,'file') || parms.forceflag
      vals_thick = fs_load_mgh(fnames_thick{h});
      vals_area  = fs_load_mgh(fnames_area{h});
      fs_save_mgh(vals_thick .* vals_area,fname_vol);
    end
  end
  % smooth cortvol on sphere for each required smoothing level
  smoothvec = union(parms.sphsmoothsteps,parms.fuzzy_smooth);
  smoothvec = setdiff(smoothvec,0);
  for s=1:length(smoothvec)
    smooth = smoothvec(s);
    try
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fname_in  = sprintf('%s-%s.%s',outstem,hemi,parms.outtype);
        fname_out = sprintf('%s-sm%d-%s.%s',outstem,smooth,hemi,parms.outtype);
        fs_surf2surf(fname_in,parms.icosubj,...
          'fname_out',fname_out,...
          'hemi',hemi,...
          'outtype',parms.outtype,...
          'smooth_in',smooth,...
          'cortex_flag',true,...
          'subjdir',parms.icosubjdir,...
          'verbose',0,...
          'forceflag',parms.forceflag);
      end
    catch me
      errmsg = sprintf('failed to smooth cortvol: %s',me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1;
      return;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = aparc_surf_analysis(parms,meas)
  errcode = 0; errmsg = [];
  outstem = sprintf('%s/%s',parms.tmpdir,meas);
  for i=1:length(parms.aparc_infix)
    aparc_infix = parms.aparc_infix{i};
    fname_mat = sprintf('%s_%s_roi_data.mat',outstem,aparc_infix);
    if ~exist(fname_mat,'file') || parms.forceflag
      roi_data = [];
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fname_vals = sprintf('%s-%s.%s',outstem,hemi,parms.outtype);
        fname_aparc = sprintf('%s/label/%s.%s.annot',parms.fspath,hemi,aparc_infix);
        try
          tmp_roi_data = mmil_surf_roi(fname_vals,...
            'fname_aparc',fname_aparc,...
            'verbose',1);
        catch me
          errmsg = sprintf('aparc surf analysis failed for %s %s: %s',...
            meas,aparc_infix,me.message);
          fprintf('%s: ERROR: %s\n',mfilename,errmsg);
          errcode = 1;
          return;
        end
        roi_data = [roi_data,tmp_roi_data];
      end
      save(fname_mat,'roi_data');
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = gwc_analysis(parms)
  errcode = 0; errmsg = [];
  for i=1:length(parms.proc_inputlist)
    fname_in = parms.proc_inputlist{i};
    outstem = sprintf('%s/%s',parms.tmpdir,parms.proc_outputlist{i});
    scalefact = parms.proc_scalefacts(i);
    for j=1:length(parms.aparc_infix)
      aparc_infix = parms.aparc_infix{j};
      fnames_aparc = cell(parms.nhemi,1);
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        fnames_aparc{h} = sprintf('%s/label/%s.%s.annot',...
          parms.fspath,hemi,aparc_infix);
      end
      try
        mmil_cortsurf_analysis(fname_in,parms.fspath,...
          'outdir',parms.tmpdir,...
          'outstem',outstem,...
          'outfix',aparc_infix,...
          'fnames_aparc',fnames_aparc,...
          'fnames_weights',[],...
          'projdist_list',parms.projdist_list,...
          'smoothsteps',0,...
          'sphere_flag',1,...
          'sphsmoothsteps',parms.sphsmoothsteps,...
          'mask_midbrain_flag',0,...
          'scalefact',scalefact,...
          'verbose',1,...
          'forceflag',parms.forceflag);
      catch me
        errmsg = sprintf('gwc analysis failed for %s: %s',...
          parms.proc_outputlist{i},me.message);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        errcode = 1;
        return;
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [errcode,errmsg] = fuzzy_analysis(parms)
  errcode = 0; errmsg = [];

  % fuzzy analysis for surface measures (thickness, area, sulc)
  measlist = {'thickness','area','sulc'};
  for mi=1:length(measlist)
    meas = measlist{mi};
    outstem = sprintf('%s/%s',parms.tmpdir,meas);
    try
      fnames_out = fs_paint(parms.subj,[],...
        'meas',meas,...
        'outstem',outstem,...
        'outtype',parms.outtype,...
        'subjdir',parms.subjdir,...
        'smoothsteps',0,...
        'sphere_flag',1,...
        'sphsmoothsteps',parms.fuzzy_smooth,...
        'mask_midbrain_flag',0,...
        'hemilist',parms.hemilist,...
        'forceflag',parms.forceflag);
    catch me
      errmsg = sprintf('fuzzy resampling failed for %s: %s',meas,me.message);
      fprintf('%s: ERROR: %s\n',mfilename,errmsg);
      errcode = 1; return;
    end
    fname_mat = sprintf('%s_%s_roi_data.mat',outstem,parms.fuzzy_outfix);
    if ~exist(fname_mat,'file') || parms.forceflag
      roi_data = [];
      for h=1:parms.nhemi
        hemi = parms.hemilist{h};
        try
          tmp_roi_data = mmil_surf_roi(fnames_out{h},...
            'hemi',hemi,...
            'fname_weights',parms.fnames_weights{h},...
            'weights_thresh',parms.weights_thresh,...
            'verbose',1);
        catch me
          errmsg = sprintf('fuzzy ROI analysis failed for %s: %s',meas,me.message);
          fprintf('%s: ERROR: %s\n',mfilename,errmsg);
          errcode = 1; return;
        end
        tmp_roi_data = replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},hemi);
        roi_data = [roi_data,tmp_roi_data];
      end
      save(fname_mat,'roi_data');
    end
  end

  % fuzzy analysis for cortical volume
  instem = sprintf('%s/cortvol-sphere',parms.tmpdir);
  outstem = sprintf('%s/cortvol',parms.tmpdir);
  fname_mat = sprintf('%s_%s_roi_data.mat',outstem,parms.fuzzy_outfix);
  if ~exist(fname_mat,'file') || parms.forceflag
    roi_data = [];
    for h=1:parms.nhemi
      hemi = parms.hemilist{h};
      if parms.fuzzy_smooth > 0
        fname_in = sprintf('%s-sm%d-%s.%s',instem,parms.fuzzy_smooth,hemi,parms.outtype);
      else
        fname_in = sprintf('%s-%s.%s',instem,hemi,parms.outtype);
      end
      try
        tmp_roi_data = mmil_surf_roi(fname_in,...
          'hemi',hemi,...
          'fname_weights',parms.fnames_weights{h},...
          'weights_thresh',parms.weights_thresh,...
          'verbose',1);
      catch me
        errmsg = sprintf('fuzzy ROI analysis failed for cortvol: %s',me.message);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        errcode = 1; return;
      end
      tmp_roi_data = replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},hemi);
      roi_data = [roi_data,tmp_roi_data];
    end
    save(fname_mat,'roi_data');
  end

  % fuzzy analysis for proc files (T1w, T2w)
  for i=1:length(parms.proc_inputlist)
    fname_in = parms.proc_inputlist{i};
    for projdist=parms.projdist_list
      outstem = sprintf('%s/%s_pdist%0.1f',...
        parms.tmpdir,parms.proc_outputlist{i},projdist);
      try
        fnames_out = fs_paint(parms.subj,fname_in,...
          'outstem',outstem,...
          'outtype',parms.outtype,...
          'subjdir',parms.subjdir,...
          'projdist',projdist,...
          'smoothsteps',0,...
          'sphere_flag',1,...
          'sphsmoothsteps',parms.fuzzy_smooth,...
          'mask_midbrain_flag',0,...
          'hemilist',parms.hemilist,...
          'forceflag',parms.forceflag);
      catch me
        errmsg = sprintf('fuzzy resampling failed for %s: %s',...
          parms.proc_outputlist{i},me.message);
        fprintf('%s: ERROR: %s\n',mfilename,errmsg);
        errcode = 1; return;
      end
      fname_mat = sprintf('%s/%s_%s_pdist%0.1f_roi_data.mat',...
        parms.tmpdir,parms.proc_outputlist{i},parms.fuzzy_outfix,projdist);
      if ~exist(fname_mat,'file') || parms.forceflag
        roi_data = [];
        for h=1:parms.nhemi
          hemi = parms.hemilist{h};
          try
            tmp_roi_data = mmil_surf_roi(fnames_out{h},...
              'hemi',hemi,...
              'fname_weights',parms.fnames_weights{h},...
              'weights_thresh',parms.weights_thresh,...
              'verbose',1);
          catch me
            errmsg = sprintf('fuzzy ROI analysis failed for %s: %s',...
              parms.proc_outputlist{i},me.message);
            fprintf('%s: ERROR: %s\n',mfilename,errmsg);
            errcode = 1; return;
          end
          tmp_roi_data = replace_fuzzy_names(tmp_roi_data,parms.fuzzy_names{h},hemi);
          roi_data = [roi_data,tmp_roi_data];
        end
        save(fname_mat,'roi_data');
      end
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function roi_data = replace_fuzzy_names(roi_data,fuzzy_names,hemi)
  if isempty(fuzzy_names), return; end
  for i=1:length(roi_data)
    roi_data(i).roiname = sprintf('ctx-%s-%s',hemi,fuzzy_names{i});
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


