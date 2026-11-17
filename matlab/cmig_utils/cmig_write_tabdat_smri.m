function cmig_write_tabdat_smri(varargin)
%function cmig_write_tabdat_smri([options])
%
% Optional Input:
%   'indir': input directory containing ROI averages in mat files
%     {default = pwd}
%   'participant_id': participant ID
%     {default = []}
%   'session_id': session ID
%     {default = []}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     {default = 'mr_y_smri'}
%   'colstem': prefix attached to column names
%     {default = 'mr_y_smri'}
%   'meas': measure name; e.g., 'area', 'thk', 'vol', 'sulc', 't1', 't2'
%     {default = 'area'}
%   'atlas': atlas name / ROI type
%     'dsk' = Desikan aparc
%     'dst' = Destrieux aparc
%     'fzy' = fuzzy clusters
%     'aseg' = automated brain segmentation
%     {default = 'dsk'}
%   'gwc_flag': [0|1] whether to load averages for gray and white matter
%     and calculate gray-white contrast
%     NOTE: a value of 1 only works for measures like 't1' and 't2'
%     {default = 0}
%   'fuzzy_dir': input directory for fuzzy cluster ROIs
%     only relevant if atlas = 'fzy'
%     if empty, will use [cmig_getenv('CMIG_PIPELINE_ATLASDIR') '/fuzzy_clusters']
%     {default = []}
%   'fuzzy_fstem': input file stem for fuzzy cluster ROIs
%     only relevant if atlas = 'fzy'
%     with expected names like {fstem}{order}-{hemi}.mgz
%     {default = 'fuzzy'}
%   'fuzzy_order': [0|2|4|12|18] number of fuzzy cluster ROIs
%     only relevant if atlas = 'fzy'
%     note: set of 18 includes combined sets of 2, 4, and 12
%     if order=0, names are like {fstem}-{hemi}.mgz
%     {default = 18}
%   'fs_version': FreeSurfer version number; e.g., 450, 530, 711, 810
%     {default = 711}
%   'tabdir': directory containing ROI abbreviation lookup table CSV files
%     {default = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/tabulated']}
%   'fname_colorlut': color lookup table of the custom region codes, combined
%     with FreeSurfer's by mmil_colorlut
%     {default = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/CMIG_ColorLUT.txt']}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Created:  11/04/2025 by Don Hagler
% Prev Mod: 11/07/2025 by Don Hagler
% Prev Mod: 11/13/2025 by Don Hagler
% Prev Mod: 11/14/2025 by Don Hagler
% Prev Mod: 11/17/2025 by Don Hagler
% Prev Mod: 11/18/2025 by Don Hagler
% Prev Mod: 05/20/2026 by Don Hagler
% Prev Mod: 05/28/2026 by Don Hagler
% Prev Mod: 09/22/2026 by Don Hagler
% Prev Mod: 09/28/2026 by Don Hagler
% Prev Mod: 10/26/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',pwd,[],...
  'participant_id',[],[],...
  'session_id',[],[],...
...
  'outdir',pwd,[],...
  'outstem','mr_y_smri',[],...
  'colstem','mr_y_smri',[],...
  'fname_out',[],[],...
  'meas','area',[],...
  'atlas','dsk',{'dsk','dst','fzy','aseg'},...
  'gwc_flag',false,[false true],...
  'fs_version',711,[],...
  'tabdir',[],[],...
  'fname_colorlut',[],[],...
  'forceflag',false,[false true],...
...
  'aseg_roilist',[2,4:5,7,8,10:18,24,26,28,41,43:44,46,47,49:54,58,60,77,251:255],[1,Inf],...
  'extra_roilist',[2,3,41,42,21001:21003,20003,20004,20008,20009,20010],[1,Inf],...
  'aparc_roilist',[1001:1003,1005:1035,2001:2003,2005:2035],[],...
  'aparc2_roilist',[11101:11141,11143:11175,12101:12141,12143:12175],[],...
...
  'wm_projdist',-0.2,[-5,0],...
  'gm_projdist',0.2,[0,5],...
...
  'fuzzy_dir',[],[],...
  'fuzzy_fstem','fuzzy',[],...
  'fuzzy_order',18,[0,2,4,12,18],...
...
  'erode_flag',false,[false true],...
...
  'ext','.tsv',[],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

% the cmig color LUT, whose custom region codes clear FreeSurfer 8.x's
if isempty(parms.fname_colorlut)
  parms.fname_colorlut = [cmig_getenv('CMIG_PIPELINE_DIR') '/parms/CMIG_ColorLUT.txt'];
end

% check if output file(s) exists
runflag = check_output(parms);
if ~runflag, return; end

% get roi-related parameters
parms = get_roi_info(parms);

% write tabulated data files
if strcmp(parms.atlas,'aseg')
  write_tabdat_aseg(parms);
elseif parms.gwc_flag
  write_tabdat_gwc(parms);
else
  write_tabdat_surf(parms);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runflag = check_output(parms)
  runflag = 0;
  if parms.gwc_flag && ~strcmp(parms.atlas,'aseg')
    suffix_list = {'wm','gm','gwc'};
    for s=1:length(suffix_list)
      if ~isempty(parms.fname_out)
        [d,n,e] = fileparts(parms.fname_out);
        fname_out = fullfile(d,[n,suffix_list{s},e]);
      else
        outstem = sprintf('%s__%s__%s__%s',parms.outstem,lower(parms.meas),suffix_list{s},lower(parms.atlas));
        fname_out = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
      end
      if ~exist(fname_out,'file') || parms.forceflag, runflag = 1; end
    end
  else
    if ~isempty(parms.fname_out)
      fname_out = parms.fname_out;
    else
      outstem = sprintf('%s__%s__%s',parms.outstem,lower(parms.meas),lower(parms.atlas));
      fname_out = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
    end
    if ~exist(fname_out,'file') || parms.forceflag, runflag = 1; end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = get_roi_info(parms)
  % get ROI name abbreviations used for tabulated data
  [roi_names, roi_abbrevs] = cmig_load_abbrev('indir', parms.tabdir);

  % set ROI codes and names
  switch parms.atlas
    case {'dsk','dst'}
      parms = set_aparc_rois(parms);
    case 'fzy'
      parms = set_fuzzy_rois(parms);
    case 'aseg'
      parms = set_aseg_rois(parms);
  end

  % replace ROI names with abbreviations
  parms.shortnames = cell(size(parms.roinames));
  for i=1:length(parms.roinames)
    idx_roi = find(strcmp(parms.roinames{i},roi_names));
    if isempty(idx_roi)
      fprintf('%s: WARNING: no ROI in look-up table called %s\n',mfilename,parms.roinames{i});
      parms.shortnames{i} = parms.roinames{i};
    else
      parms.shortnames{i} = roi_abbrevs{idx_roi};
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set aparc ROI codes and names depending on atlas ('dsk' or 'dst')
function parms = set_aparc_rois(parms)
  % read roi codes and names from FreeSurfer color lookup table
  %   and combine with the cmig color LUT
  [all_roicodes,all_roinames] = mmil_colorlut('FS_version',parms.fs_version,...
    'fname_excolorlut',parms.fname_colorlut);
  % set file name infix and roilist
  switch parms.atlas
    case 'dsk'
      parms.infix = 'aparc';
      roilist = parms.aparc_roilist;
    case 'dst'
      parms.infix = 'aparc.a2009s';
      roilist = parms.aparc2_roilist;
  end
  % match with master list of roicodes to get roinames
  [tmp,ind_rois] = intersect(all_roicodes,roilist);
  parms.roicodes = mmil_rowvec(all_roicodes(ind_rois));
  parms.roinames = mmil_rowvec(all_roinames(ind_rois));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set ROI codes and names for fuzzy cluster ROIs
function parms = set_fuzzy_rois(parms)
  hemilist = {'lh','rh'};
  if isempty(parms.fuzzy_dir)
    parms.fuzzy_dir = [cmig_getenv('CMIG_PIPELINE_ATLASDIR') '/fuzzy_clusters'];
  end
  if ~exist(parms.fuzzy_dir,'dir')
    error('fuzzy cluster dir %s not found',parms.fuzzy_dir);
  end
  if parms.fuzzy_order~=0
    parms.infix = ...
      sprintf('%s%d',parms.fuzzy_fstem,parms.fuzzy_order);
    fuzzy_names = mmil_fuzzy_names('fuzzy_fstem',parms.fuzzy_fstem,...
                                   'fuzzy_order',parms.fuzzy_order);
    parms.roicodes = [];
    parms.roinames = [];
    for h=1:length(hemilist)
      hemi = hemilist{h};
      tmp_codes = 20000 + (h-1)*5000 + ...
                  parms.fuzzy_order*100 + [1:parms.fuzzy_order]';
      tmp_names = cell(parms.fuzzy_order,1);
      for f=1:parms.fuzzy_order
        tmp_names{f} = sprintf('ctx-%s-%s',hemi,fuzzy_names{f});
      end
      parms.roicodes = cat(2,parms.roicodes,mmil_rowvec(tmp_codes));
      parms.roinames = cat(2,parms.roinames,mmil_rowvec(tmp_names));
    end
  else
    parms.infix = parms.fuzzy_fstem;
    parms.roicodes = [];
    parms.roinames = [];
    for h=1:length(hemilist)
      hemi = hemilist{h};
      fname_weights = sprintf('%s/%s-%s.mgz',...
        parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
      if ~exist(fname_weights,'file')
        error('fuzzy cluster file %s not found',fname_weights);
      end
      [M,volsz] = fs_read_header(fname_weights);
      nrois = volsz(4);
      tmp_codes = 20000 + (h-1)*5000 + ...
                  nrois*100 + [1:nrois]';
      % load roinames
      fname_txt = sprintf('%s/%s_%s_roinames.txt',...
        parms.fuzzy_dir,parms.fuzzy_fstem,hemi);
      if ~exist(fname_txt,'file')
        fprintf('%s: WARNING: roinames file %s not found',...
          fname_txt);
        tmp_names = cell(nrois,1);
        for f=1:nrois
          tmp_names{f} = sprintf('ctx-%s_%d-%s',hemi,parms.fuzzy_fstem,f);
        end
      else
        tmp_names = mmil_readtext(fname_txt);
      end
      parms.roicodes = cat(2,parms.roicodes,mmil_rowvec(tmp_codes));
      parms.roinames = cat(2,parms.roinames,mmil_rowvec(tmp_names));
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set aparc ROI codes and names for aseg and extra
function parms = set_aseg_rois(parms)
  % read roi codes and names from FreeSurfer color lookup table
  %   and combine with the cmig color LUT
  [all_roicodes,all_roinames] = mmil_colorlut('FS_version',parms.fs_version,...
    'fname_excolorlut',parms.fname_colorlut);
  % set file name infix
  parms.inix = 'aseg';
  % set roilist depending measure
  if strcmp(parms.meas,'vol')
    roilist = union(parms.aseg_roilist,parms.extra_roilist);
  else
    roilist = parms.aseg_roilist;
  end
  % match with master list of roicodes to get roinames
  [tmp,ind_rois] = intersect(all_roicodes,roilist);
  parms.roicodes = mmil_rowvec(all_roicodes(ind_rois));
  parms.roinames = mmil_rowvec(all_roinames(ind_rois));
  % set string included in eroded aseg roi_data files
  if parms.erode_flag
    parms.erode_outfix = '_erode';
    if parms.erode_nvoxels>1
      parms.erode_outfix = sprintf('%s%d',...
        parms.erode_outfix,parms.erode_nvoxels);
    end
  else
    parms.erode_outfix = [];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tabdat_surf(parms)
  % load data and write tsv file
  colstem = lower(sprintf('%s__%s__%s',parms.colstem,parms.meas,parms.atlas));
  if ~isempty(parms.fname_out)
    fname_out = parms.fname_out;
  else
    outstem = sprintf('%s__%s__%s',parms.outstem,lower(parms.meas),lower(parms.atlas));
    fname_out = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
  end
  data = [];
  mmil_mkdir(parms.outdir);
  if ismember(parms.meas,{'area','thk','vol'}) && ~strcmp(parms.atlas,'fzy')
    % load aparc_stats file
    data = load_aparc_stats(parms,parms.meas);
  else
    % load roi_data file
    data = load_roi_data(parms,parms.meas);
  end
  if isempty(data), return; end
  if ~strcmp(parms.atlas,'fzy')
    % calculate global mean or total values and append to data
    if ismember(parms.meas,{'area','vol'})
      [data,shortnames] = calc_sums(data,parms.shortnames);
    else
      % load area for aparc ROIs
      area_data = load_aparc_stats(parms,'area');
      % calculate sum or mean values weighted by area
      [data,shortnames] = calc_means(data,parms.shortnames,area_data);
    end
  else
    shortnames = parms.shortnames;
  end
  % write tsv file
  write_tsv_file(fname_out,data,shortnames,colstem,parms)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tabdat_gwc(parms)
  mmil_mkdir(parms.outdir);
  % load wm roi_data file
  wm_data = load_roi_data(parms,parms.meas,parms.wm_projdist);
  % load gm roi_data file
  gm_data = load_roi_data(parms,parms.meas,parms.gm_projdist);
  % check data
  if isempty(wm_data) || isempty(gm_data), error('empty data'); end
  % calculate gray-white contrast
  gwc_data = calc_contrast(wm_data,gm_data);
  if ~strcmp(parms.atlas,'fzy')
    % load area
    area_data = load_aparc_stats(parms,'area');
    % calculate mean values weighted by area
    [wm_data,wm_shortnames] = calc_means(wm_data,parms.shortnames,area_data);
    [gm_data,gm_shortnames] = calc_means(gm_data,parms.shortnames,area_data);
    [gwc_data,gwc_shortnames] = calc_means(gwc_data,parms.shortnames,area_data);
  else
    wm_shortnames = parms.shortnames;
    gm_shortnames = parms.shortnames;
    gwc_shortnames = parms.shortnames;
  end
  % write wm tsv file
  colstem = lower(sprintf('%s__%s__wm__%s',parms.colstem,parms.meas,parms.atlas));
  if ~isempty(parms.fname_out)
    [d,n,e] = fileparts(parms.fname_out); fname_wm = fullfile(d,[n,'wm',e]);
  else
    outstem = sprintf('%s__%s__wm__%s',parms.outstem,lower(parms.meas),lower(parms.atlas));
    fname_wm = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
  end
  write_tsv_file(fname_wm,wm_data,wm_shortnames,colstem,parms);
  % write gm tsv file
  colstem = lower(sprintf('%s__%s__gm__%s',parms.colstem,parms.meas,parms.atlas));
  if ~isempty(parms.fname_out)
    [d,n,e] = fileparts(parms.fname_out); fname_gm = fullfile(d,[n,'gm',e]);
  else
    outstem = sprintf('%s__%s__gm__%s',parms.outstem,lower(parms.meas),lower(parms.atlas));
    fname_gm = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
  end
  write_tsv_file(fname_gm,gm_data,gm_shortnames,colstem,parms);
  % write gwc tsv file
  colstem = lower(sprintf('%s__%s__gwc__%s',parms.colstem,parms.meas,parms.atlas));
  if ~isempty(parms.fname_out)
    [d,n,e] = fileparts(parms.fname_out); fname_gwc = fullfile(d,[n,'gwc',e]);
  else
    outstem = sprintf('%s__%s__gwc__%s',parms.outstem,lower(parms.meas),lower(parms.atlas));
    fname_gwc = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
  end
  write_tsv_file(fname_gwc,gwc_data,gwc_shortnames,colstem,parms);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tabdat_aseg(parms)
  % load data and write tsv file
  colstem = lower(sprintf('%s__%s__%s',parms.colstem,parms.meas,parms.atlas));
  if ~isempty(parms.fname_out)
    fname_out = parms.fname_out;
  else
    outstem = sprintf('%s__%s__%s',parms.outstem,lower(parms.meas),lower(parms.atlas));
    fname_out = sprintf('%s/%s%s',parms.outdir,outstem,parms.ext);
  end
  data = [];
  mmil_mkdir(parms.outdir);
  if strcmp(parms.meas,'vol')
    % load aseg_stats file
    data = load_aseg_stats(parms);
  else
    % load roi_data file
    data = load_aseg_data(parms,parms.meas);
  end
  if isempty(data), return; end
  % write tsv file
  write_tsv_file(fname_out,data,parms.shortnames,colstem,parms)
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = load_aparc_stats(parms,meas)
  data = [];
  fname_mat = sprintf('%s/%s_stats.mat',parms.indir,parms.infix);
  if ~exist(fname_mat,'file')
    fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
  else
    aparc_stats = [];
    load(fname_mat);
    if isempty(aparc_stats)
      fprintf('%s: WARNING: empty aparc_stats in %s\n',mfilename,fname_mat);
    else
      % match ROIs by name
      tmp_roinames = {aparc_stats.roiname};
      [tmp,ind_all,ind_subj]=intersect(parms.roinames,tmp_roinames);
      % select datatype
      switch meas
        case 'area'
          tmp_data = cell2mat({aparc_stats(ind_subj).surfarea});
        case 'thk'
          tmp_data = cell2mat({aparc_stats(ind_subj).thickavg});
        case 'vol'
          tmp_data = cell2mat({aparc_stats(ind_subj).grayvol});
        otherwise
          error('unsupported meas %s',meas);
      end
      data = nan(1,length(parms.roinames));
      data(1,ind_all) = tmp_data;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = load_roi_data(parms,meas,projdist)
  data = [];
  if ~exist('projdist','var'), projdist = []; end
  switch meas
    case 'thk'
      instem = 'thickness';
    case 'vol'
      instem = 'cortvol';
    case 't1'
      instem = 'T1w';
    case 't2'
      instem = 'T2w';
    otherwise
      instem = meas;
  end
  if isempty(projdist)
    fname_mat = sprintf('%s/%s_%s_roi_data.mat',parms.indir,instem,parms.infix);
  else
    fname_mat = sprintf('%s/%s_%s_pdist%0.1f_roi_data.mat',parms.indir,instem,parms.infix,projdist);
  end
  if ~exist(fname_mat,'file')
    fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
  else
    roi_data = [];
    load(fname_mat);
    if isempty(roi_data)
      fprintf('%s: WARNING: empty roi_data in %s\n',mfilename,fname_mat);
    else
      % match ROIs by name
      tmp_roinames = {roi_data.roiname};
      [tmp,ind_all,ind_subj]=intersect(parms.roinames,tmp_roinames);
      tmp_data = cell2mat({roi_data(ind_subj).avg});
      data = nan(1,length(parms.roinames));
      data(1,ind_all) = tmp_data;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = load_aseg_stats(parms)
  % load aseg_stats
  fname_mat = sprintf('%s/aseg_stats.mat',parms.indir);
  aseg_stats = [];
  if ~exist(fname_mat,'file')
    fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
  else
    load(fname_mat);
    if isempty(aseg_stats)
      fprintf('%s: WARNING: empty aseg_stats in %s\n',...
        mfilename,fname_mat);
    else
      % match ROIs by code
      tmp_roicodes = cell2mat({aseg_stats.roicode});
      [tmp,ind_all,ind_subj]=intersect(parms.roicodes,tmp_roicodes);
      tmp_data = cell2mat({aseg_stats(ind_subj).volume});
      data = nan(1,length(parms.roinames));
      data(1,ind_all) = tmp_data;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = load_aseg_data(parms,meas)
  data = [];
  switch meas
    case 't1'
      instem = 'T1w';
    case 't2'
      instem = 'T2w';
    otherwise
      instem = meas;
  end
  % load average aseg ROI data
  fname_mat = sprintf('%s/%s_aseg%s_roi_data.mat',...
    parms.indir,instem,parms.erode_outfix);
  roi_data = [];
  if ~exist(fname_mat,'file')
    fprintf('%s: WARNING: missing file %s\n',mfilename,fname_mat);
  else
    roi_data = [];
    load(fname_mat);
    if isempty(roi_data)
      fprintf('%s: WARNING: missing %s aseg roi_data in %s\n',...
        mfilename,meas,fname_mat);
    else
      % match ROIs by name
      tmp_roinames = {roi_data.roiname};
      [tmp,ind_all,ind_subj]=intersect(parms.roinames,tmp_roinames);
      tmp_data = cell2mat({roi_data(ind_subj).avg});
      data = nan(1,length(parms.roinames));
      data(1,ind_all) = tmp_data;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gwc_data = calc_contrast(wm_data,gm_data)
  gwc_data = (wm_data - gm_data) ./ ((wm_data + gm_data)/2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,roinames] = calc_sums(data,roinames)
  % identify left and right hemisphere regions
  idx_lh = find(contains(roinames,'_lh'));
  idx_rh = find(contains(roinames,'_rh'));
  % calculate sums for lh, rh, and all
  tmp_vals = [...
    sum(data(1,idx_lh),2),...
    sum(data(1,idx_rh),2),...
    sum(data,2),...
  ];
  tmp_names = {'lh_sum','rh_sum','sum'};
  % concatenate with data and roinames
  data = cat(2,data,tmp_vals);
  roinames = cat(2,roinames,tmp_names);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,roinames] = calc_means(data,roinames,area_data)
  % identify left and right hemisphere regions
  idx_lh = find(contains(roinames,'_lh'));
  idx_rh = find(contains(roinames,'_rh'));
  % calculate mean weighted by area
  tmp_vals = [...
    mmil_wtd_mean(data(1,idx_lh),area_data(1,idx_lh),2),...
    mmil_wtd_mean(data(1,idx_rh),area_data(1,idx_rh),2),...
    mmil_wtd_mean(data,area_data,2),...
  ];
  tmp_names = {'lh_mean','rh_mean','mean'};
  % concatenate with data and roinames
  data = cat(2,data,tmp_vals);
  roinames = cat(2,roinames,tmp_names);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tsv_file(fname_out,data,roinames,prefix,parms)
  if ~isempty(data)
    if strcmp(parms.atlas,'fzy')
      vartype = 'wmean';
    elseif ismember(parms.meas,{'area','vol'})
      vartype = 'sum';
    else
      vartype = 'mean';
    end
    % prepend prefix to column labels
    colnames = strcat(prefix,'__',roinames,'_',vartype);
    % remove extra underscore and "_sum" (or "_mean") if roiname is "sum" (or "mean")
    colnames = regexprep(colnames,sprintf('_+%s_%s',vartype,vartype),sprintf('_%s',vartype));
    % add participant and session IDs as initial columns
    colnames = cat(2,{'participant_id','session_id'},colnames);
    data = cat(2,{parms.participant_id,parms.session_id},num2cell(data));
    % add column labels as first row
    data = cat(1,colnames,data);
    % write to file
    mmil_write_csv(fname_out,data,'separator','\t','quote_flag',0);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


