function cmig_write_tabdat_dmri(varargin)
%function cmig_write_tabdat_dmri([options])
%
% Purpose: write one tabulated (TSV) file of dMRI ROI measures from a
%   roi_data .mat file produced by cmig_analyze_bids_dwi
%
% Optional Input ('key',value):
%   'fname_in': full path of the roi_data .mat file to tabulate
%     {default = []}
%   'fname_out': full path of the output TSV file
%     {default = []}
%   'participant_id': participant ID string
%     {default = []}
%   'session_id': session ID string
%     {default = []}
%   'colstem': prefix attached to column names
%     {default = ''}
%   'atlas': ROI type; selects the column-name convention
%     'aseg' = FreeSurfer aseg subcortical ROIs (mean)
%     'dsk'  = Desikan gwcsurf ROIs (mean)
%     'dst'  = Destrieux gwcsurf ROIs (mean)
%     'at'   = AtlasTrack fiber tracts (weighted mean)
%     {default = 'aseg'}
%   'tabdir': directory with ROI abbreviation lookup CSV files (aseg/dsk/dst)
%     if empty, uses $MMPS_PARMS/tabulated (via cmig_load_abbrev)
%     {default = []}
%   'FiberLegend': fiber legend struct (for atlas='at'); if empty, loaded via
%     DTI_MMIL_Read_Fiber_Legend
%     {default = []}
%   'fspath': FreeSurfer recon directory (for atlas='dsk' or 'dst'); if given,
%     lh, rh and whole-cortex means weighted by each ROI's surface area are
%     appended, the area from stats/?h.{aparc_infix}.stats
%     {default = []}
%   'aparc_infix': annotation of the gwcsurf ROIs, e.g. 'aparc' or
%     'aparc.a2009s' (for atlas='dsk' or 'dst')
%     {default = []}
%
% Created:  07/22/2026 by Don Hagler
% Prev Mod: 07/22/2026 by Don Hagler
% Prev Mod: 09/28/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'fname_in',[],[],...
  'fname_out',[],[],...
  'participant_id',[],[],...
  'session_id',[],[],...
  'colstem','',[],...
  'atlas','aseg',{'aseg','dsk','dst','at'},...
  'tabdir',[],[],...
  'FiberLegend',[],[],...
  'fspath',[],[],...
  'aparc_infix',[],[],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

if isempty(parms.fname_in) || ~exist(parms.fname_in,'file')
  fprintf('%s: WARNING: missing roi_data: %s\n',mfilename,parms.fname_in);
  return;
end

roi_data = [];
load(parms.fname_in);
if isempty(roi_data), return; end

colnames = {}; avgdata = [];
switch parms.atlas
  case 'aseg'
    [roi_names,roi_abbrevs] = cmig_load_abbrev('indir',parms.tabdir);
    [colnames,avgdata] = build_aseg_colnames(roi_data,roi_names,roi_abbrevs,parms.colstem,mfilename);
  case {'dsk','dst'}
    [roi_names,roi_abbrevs] = cmig_load_abbrev('indir',parms.tabdir);
    [colnames,avgdata,roinames] = build_gwcsurf_colnames(roi_data,roi_names,roi_abbrevs,parms.colstem,mfilename);
    if ~isempty(parms.fspath) && ~isempty(colnames)
      [colnames,avgdata] = add_hemi_means(colnames,avgdata,roinames,parms);
    end
  case 'at'
    FiberLegend = parms.FiberLegend;
    if isempty(FiberLegend)
      try
        FiberLegend = DTI_MMIL_Read_Fiber_Legend();
      catch err
        fprintf('%s: WARNING: could not load fiber legend: %s\n',mfilename,err.message);
      end
    end
    [colnames,avgdata] = build_fiber_colnames(roi_data,FiberLegend,parms.colstem,mfilename);
end

if isempty(colnames), return; end

mmil_mkdir(fileparts(parms.fname_out));
write_tsv(parms.fname_out,parms.participant_id,parms.session_id,colnames,avgdata);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [colnames,avgdata] = build_aseg_colnames(roi_data,roi_names,roi_abbrevs,colstem,caller)
  colnames = {};
  avgdata = [];
  nrois = length(roi_data);
  tmp_colnames = cell(1,nrois);
  tmp_avgs = zeros(1,nrois);
  valid = false(1,nrois);
  for r = 1:nrois
    roiname = roi_data(r).roiname;
    idx = find(strcmp(roi_names, roiname));
    if isempty(idx)
      fprintf('%s: WARNING: roiname not in abbrev table: %s\n',caller,roiname);
      shortname = roiname;
    else
      shortname = roi_abbrevs{idx(1)};
    end
    tmp_colnames{r} = sprintf('%s__%s_mean',colstem,shortname);
    tmp_avgs(r) = roi_data(r).avg;
    valid(r) = true;
  end
  colnames = tmp_colnames(valid);
  avgdata = tmp_avgs(valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [colnames,avgdata,roinames] = build_gwcsurf_colnames(roi_data,roi_names,roi_abbrevs,colstem,caller)
  % same lookup as build_aseg_colnames; roi_abbrev table entries include hemi suffix
  colnames = {};
  avgdata = [];
  roinames = {};
  nrois = length(roi_data);
  tmp_colnames = cell(1,nrois);
  tmp_avgs = zeros(1,nrois);
  valid = false(1,nrois);
  % known non-cortical surface labels (medial wall) to skip silently: not
  % real ROIs and intentionally absent from the abbrev table; as production's
  % aparc_roilist and aparc2_roilist (MMIL_Compile_DTI_Analysis), which leave
  % out unknown and corpuscallosum (1000, 1004) and Unknown and Medial_wall
  % (11100, 11142)
  skip_labels = {'unknown','corpuscallosum','medial_wall'};
  for r = 1:nrois
    roiname = roi_data(r).roiname;
    % strip ctx-{lh,rh}- / ctx_{lh,rh}_ prefix to get the base label
    base_label = regexprep(roiname,'^ctx[-_](lh|rh)[-_]','');
    if any(strcmpi(base_label,skip_labels))
      continue;  % skip silently: no warning, no column
    end
    idx = find(strcmp(roi_names, roiname));
    if isempty(idx)
      fprintf('%s: WARNING: gwcsurf roiname not in abbrev table: %s\n',caller,roiname);
      shortname = roiname;
    else
      shortname = roi_abbrevs{idx(1)};
    end
    tmp_colnames{r} = sprintf('%s__%s_mean',colstem,shortname);
    tmp_avgs(r) = roi_data(r).avg;
    valid(r) = true;
  end
  colnames = tmp_colnames(valid);
  avgdata = tmp_avgs(valid);
  roinames = {roi_data(valid).roiname};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [colnames,avgdata] = add_hemi_means(colnames,avgdata,roinames,parms)
  % append lh, rh and whole-cortex means weighted by ROI surface area, as
  %   cmig_write_tabdat_smri does for the sMRI gwc measures
  [subjdir,subj] = fileparts(parms.fspath);
  aparc_infix = parms.aparc_infix;
  if isempty(aparc_infix), aparc_infix = 'aparc'; end
  aparc_stats = fs_read_aparc_stats(subj,subjdir,aparc_infix);
  if isempty(aparc_stats)
    fprintf('%s: WARNING: no %s stats in %s, skipping hemisphere means\n',...
      mfilename,aparc_infix,parms.fspath);
    return;
  end
  area = nan(1,length(roinames));
  [~,ind_roi,ind_stats] = intersect(roinames,{aparc_stats.roiname});
  area(ind_roi) = [aparc_stats(ind_stats).surfarea];
  if any(isnan(area))
    fprintf('%s: WARNING: %d of %d ROIs have no %s area, left out of the hemisphere means\n',...
      mfilename,sum(isnan(area)),length(area),aparc_infix);
  end
  ok = ~isnan(area) & ~isnan(avgdata);
  idx_lh = ok & ~cellfun(@isempty,regexp(roinames,'^ctx[-_]lh[-_]','once'));
  idx_rh = ok & ~cellfun(@isempty,regexp(roinames,'^ctx[-_]rh[-_]','once'));
  hemi_vals = [mmil_wtd_mean(avgdata(idx_lh),area(idx_lh),2),...
               mmil_wtd_mean(avgdata(idx_rh),area(idx_rh),2),...
               mmil_wtd_mean(avgdata(ok),area(ok),2)];
  hemi_names = {sprintf('%s__lh_mean',parms.colstem),...
                sprintf('%s__rh_mean',parms.colstem),...
                sprintf('%s_mean',parms.colstem)};
  colnames = cat(2,colnames,hemi_names);
  avgdata = cat(2,avgdata,hemi_vals);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [colnames,avgdata] = build_fiber_colnames(roi_data,FiberLegend,colstem,caller)
  colnames = {};
  avgdata = [];
  if isempty(roi_data), return; end
  % build roicode → fiber abbreviation map from legend
  legend_fibernums = [];
  if ~isempty(FiberLegend)
    legend_fibernums = cell2mat({FiberLegend.FiberNumber});
  end
  nfibers = length(roi_data);
  tmp_colnames = cell(1,nfibers);
  tmp_avgs = zeros(1,nfibers);
  valid = false(1,nfibers);
  for f = 1:nfibers
    if ~isfield(roi_data(f),'roicode') || ~isfield(roi_data(f),'avg')
      fprintf('%s: WARNING: roi_data(%d) missing roicode or avg field\n',caller,f);
      continue;
    end
    roicode = roi_data(f).roicode;
    % look up fiber name from legend
    if ~isempty(legend_fibernums)
      idx = find(legend_fibernums == roicode);
    else
      idx = [];
    end
    if isempty(idx)
      abbrev = sprintf('%d', roicode);
    else
      abbrev = fiber_abbrev_from_name(FiberLegend(idx(1)).FiberName);
    end
    if isempty(abbrev)
      tmp_colnames{f} = sprintf('%s_wmean', colstem);
    else
      tmp_colnames{f} = sprintf('%s__%s_wmean', colstem, abbrev);
    end
    tmp_avgs(f) = roi_data(f).avg;
    valid(f) = true;
  end
  colnames = tmp_colnames(valid);
  avgdata = tmp_avgs(valid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function abbrev = fiber_abbrev_from_name(fibername)
  % convert DTI_Fiber_Legend FiberName to ABCD 7.0 abbreviation
  % e.g. R_ATR → atr__rh, L_ATR → atr__lh, Fmaj → fmaj, AllFibers → ''
  switch fibername
    case 'AllFibers',    abbrev = '';
    case 'R_AllFib',     abbrev = 'rh';
    case 'L_AllFib',     abbrev = 'lh';
    case 'R_AllFibnoCC', abbrev = 'nocc__rh';
    case 'L_AllFibnoCC', abbrev = 'nocc__lh';
    otherwise
      if strncmp(fibername,'R_',2)
        abbrev = sprintf('%s__rh', lower(fibername(3:end)));
      elseif strncmp(fibername,'L_',2)
        abbrev = sprintf('%s__lh', lower(fibername(3:end)));
      else
        abbrev = lower(fibername);
      end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_tsv(fname_out,participant_id,session_id,colnames,avgdata)
  header = cat(2,{'participant_id','session_id'},colnames);
  datarow = cat(2,{participant_id,session_id},num2cell(avgdata));
  all_data = cat(1,header,datarow);
  mmil_write_csv(fname_out,all_data,'separator','\t','quote_flag',0);
return;
