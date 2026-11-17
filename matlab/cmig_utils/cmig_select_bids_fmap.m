function [fname_for,fname_rev,fallback_flag] = cmig_select_bids_fmap(fname_info,fnames_fmap_data,fnames_fmap_info)
%function [fname_for,fname_rev,fallback_flag] = cmig_select_bids_fmap(fname_info,fnames_fmap_data,fnames_fmap_info)
%
% Purpose: Select forward (PA, PhaseEncodingDirection j) and
%   reverse (AP, PhaseEncodingDirection j-) fmap EPI pair for a given scan
%   using B0FieldIdentifier/B0FieldSource label matching,
%   with acquisition time fallback for ambiguous cases
%
% Required Input:
%   fname_info: JSON sidecar path for the main scan (func or dwi)
%   fnames_fmap_data: cell array of fmap NIfTI file paths
%   fnames_fmap_info: cell array of fmap JSON sidecar paths (same order)
%
% Output:
%   fname_for: NIfTI path of forward (PA) fmap; [] if not found
%   fname_rev: NIfTI path of reverse (AP) fmap; [] if not found
%   fallback_flag: 1 if the scan's B0FieldSource matched no available fmap
%     label and a substitute pair was chosen, 0 otherwise
%
% Created:  06/11/2026 by Don Hagler
% Prev Mod: 08/27/2026 by Don Hagler
% Last Mod: 11/12/2026 by Don Hagler
%

fname_for = []; fname_rev = []; fallback_flag = 0;
if isempty(fnames_fmap_data) || isempty(fnames_fmap_info), return; end

% load main scan info
try
  scan_info = jsondecode(fileread(fname_info));
catch
  return;
end
scan_label = mmil_getfield(scan_info,'B0FieldSource',...
  mmil_getfield(scan_info,'B0FieldIdentifier','UNKNOWN'));

% load fmap info (B0FieldIdentifier, PhaseEncodingDirection, AcquisitionDateTime, ShimSetting)
nfmap = length(fnames_fmap_info);
fmap_info = struct('B0FieldIdentifier',cell(1,nfmap),...
                   'PhaseEncodingDirection',cell(1,nfmap),...
                   'AcquisitionDateTime',cell(1,nfmap),...
                   'ShimSetting',cell(1,nfmap));
valid = true(1,nfmap);
for i = 1:nfmap
  try
    finfo = jsondecode(fileread(fnames_fmap_info{i}));
  catch
    valid(i) = false;
    continue;
  end
  fmap_info(i).B0FieldIdentifier = mmil_getfield(finfo,'B0FieldIdentifier','UNKNOWN');
  fmap_info(i).PhaseEncodingDirection = mmil_getfield(finfo,'PhaseEncodingDirection','');
  fmap_info(i).AcquisitionDateTime = mmil_getfield(finfo,'AcquisitionDateTime','');
  fmap_info(i).ShimSetting = mmil_getfield(finfo,'ShimSetting',[]);
end
if ~all(valid)
  keep = find(valid);
  fmap_info = fmap_info(keep);
  fnames_fmap_data = fnames_fmap_data(keep);
  nfmap = length(fmap_info);
end
if nfmap == 0, return; end

% build fmap pairs: for each unique label, pair forwards with reverses by position
revflags = get_revflags(fmap_info);
fmap_labels = {fmap_info.B0FieldIdentifier};
[~,~,idx_all] = unique(fmap_labels,'stable');
n_uniq = max(idx_all);
s_fmap_for = [];
s_fmap_rev = [];
for k = 1:n_uniq
  idx = find(idx_all == k);
  f = idx(~revflags(idx));
  r = idx(revflags(idx));
  if ~isempty(f) && ~isempty(r)
    npairs = min(length(f),length(r));
    s_fmap_for(end+1:end+npairs) = f(1:npairs);
    s_fmap_rev(end+1:end+npairs) = r(1:npairs);
  end
end
if isempty(s_fmap_for), return; end

% match fmaps by B0FieldIdentifier == scan_label
% (all-UNKNOWN naturally matches all fmaps: UNKNOWN == UNKNOWN)
s_match = find(strcmp(fmap_labels,scan_label));

% select pair
if length(s_match) == 2
  % exactly one labeled AP/PA pair — use label directly
  s_for = intersect(s_match,s_fmap_for);
  s_rev = intersect(s_match,s_fmap_rev);
else
  % no pair carries the scan's label, or the set is ambiguous (multiple
  %   UNKNOWN pairs, etc.) — choose among the pairs that are available
  % NOTE: the labeled pair may have been rejected upstream, so falling back
  %   here is what keeps the scan corrected rather than skipped
  fallback_flag = isempty(s_match);
  [s_for,s_rev] = select_nearest_fmap(scan_info,fmap_info,s_fmap_for,s_fmap_rev);
end
if isempty(s_for) || isempty(s_rev), return; end

fname_for = fnames_fmap_data{s_for(1)};
fname_rev = fnames_fmap_data{s_rev(1)};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s_for,s_rev] = select_nearest_fmap(scan_info,fmap_info,s_fmap_for,s_fmap_rev)
  s_for = []; s_rev = [];
  % filter by phase-encode axis (e.g., 'j' not 'i')
  phasedir = get_phasedirs(scan_info);
  fmap_phasedirs = get_phasedirs(fmap_info);
  s_valid = find(strcmp(fmap_phasedirs,phasedir));
  if isempty(s_valid), return; end
  % prefer field maps acquired with the same shim settings as the scan
  % NOTE: a borrowed pair shimmed differently estimates a field that does not
  %   describe the scan being corrected
  s_shim = s_valid(match_shim(scan_info,fmap_info(s_valid)));
  if ~isempty(s_shim), s_valid = s_shim; end
  % compare acquisition times to find nearest fmap, in either direction
  t_diff = compare_scan_times(scan_info,fmap_info(s_valid));
  if isempty(t_diff), return; end
  [~,i_match] = min(abs(t_diff));
  s_match = s_valid(i_match);
  % find pair partner
  if ismember(s_match,s_fmap_for)
    i_pair = find(s_match == s_fmap_for);
  else
    i_pair = find(s_match == s_fmap_rev);
  end
  if isempty(i_pair), return; end
  s_for = s_fmap_for(i_pair);
  s_rev = s_fmap_rev(i_pair);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ind = match_shim(scan_info,fmap_info)
  ind = [];
  shim = mmil_getfield(scan_info,'ShimSetting',[]);
  if isempty(shim), return; end
  ind = find(arrayfun(@(x) isequal(x.ShimSetting,shim),fmap_info));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t_diff = compare_scan_times(scan_info,fmap_info)
  t_diff = [];
  try
    scan_time = scan_info.AcquisitionDateTime;
    multi_times = {fmap_info.AcquisitionDateTime};
    scan_datenum = datenum(scan_time,'yyyy-mm-ddTHH:MM:SS');
    multi_datenums = cellfun(@(x) datenum(x,'yyyy-mm-ddTHH:MM:SS'),multi_times);
    t_diff = multi_datenums - scan_datenum;
  catch
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function revflags = get_revflags(scan_info)
  dirs = {scan_info.PhaseEncodingDirection};
  revflags = contains(dirs,'-');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phasedirs = get_phasedirs(scan_info)
  dirs = {scan_info.PhaseEncodingDirection};
  phasedirs = cell(size(dirs));
  for k = 1:length(dirs)
    if ~isempty(dirs{k})
      phasedirs{k} = dirs{k}(1);
    else
      phasedirs{k} = '';
    end
  end
return;
