function results = abcd_network_corr(R,seed_roinames,varargin)
%function results = abcd_network_corr(R,seed_roinames,[options])
%
% Required Input:
%   R: correlation matrix with size = [ntargs,nseeds]
%   seed_roinames: cell array of seed ROI names
%     length must match size(R,2)
%
% Optional Parameters:
%   'targ_roinames': cell array of target ROI names
%     length must match size(R,1)
%     if not supplied, set to seed_roinames
%   'network_flag': [0|1|2|3] average correlation within networks of ROIs
%     0: no network correlations
%     1: mean correlation between each network and each ROI
%     2: mean correlation within each network
%     3: mean correlation within and between each network
%     {default = 2}
%   'networks': struct containing network definitions
%     fields of networks struct are names of networks with values
%       that are cell arrays of ROI names
%     if empty, will generate automatically based on roi_type
%     {default = []}
%   'roi_type': type of surface ROI to be used ('fparc', 'gordon', 'talrc', or 'aparc')
%     {default = 'gordon'}
%   'roinames': cell array of ROIs to be included in output
%     applies if network_flag = 1
%     if empty, use all ROIs
%     {default = []}
%   'strict_roinames_flag': [0|1] strictly match ROI names
%     otherwise use regexp to include and exclude ROIs
%     {default = 1}
%   'save_corr_flag': [0|1] save R and Z in output results
%     may be set to 0 if R is large to save memory
%     {default = 1}
%
% Output:
%   results: struct containing fields:
%                   R: correlation matrix
%                   Z: R-to-Z transformed correlation matrix
%               ntarg: number of target ROIs
%               nseed: number of seed ROIs
%               nsubs: number of subjects / visits
%       seed_roinames: cell array of seed ROI names
%       targ_roinames: cell array of target ROI names
%            networks: struct containing network membership
%       network_names: cell array of network names
%                  nn: number of networks
%            roinames: cell array of ROI names
%               nrois: number of ROIs
%                  nt: number of target networks
%        network_corr: network correlation matrix
%  network_corr_array: network correlation cell array
% network_corr_labels: cell array of network correlation labels
%
% Created:  01/01/21 by Don Hagler
% Prev Mod: 03/26/21 by Don Hagler
% Prev Mod: 04/20/21 by Don Hagler
% Prev Mod: 06/04/24 by Don Hagler
% Last Mod: 05/14/25 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = [];
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'targ_roinames',[],[],...
  'network_flag',2,[0:3],...
  'networks',[],[],...
  'roi_type','gordon',{'fparc','gordon','talrc','aparc'},...
  'roinames',[],[],...
  'strict_roinames_flag',true,[false true],...
  'save_corr_flag',true,[false true],...
});

results.ntarg = size(R,1);
results.nseed = size(R,2);
results.nsubs = size(R,3);

results.seed_roinames = seed_roinames;
if length(results.seed_roinames) ~= results.nseed
  error('first dim of R (%d) does not match length of seed_roinames (%d)',...
    results.nseed,length(results.seed_roinames));
end;
results.targ_roinames = parms.targ_roinames;
if isempty(results.targ_roinames)
  results.targ_roinames = results.seed_roinames;
end;
if length(results.targ_roinames) ~= results.ntarg
  error('second dim of R (%d) does not match length of targ_roinames (%d)',...
    results.ntarg,length(results.targ_roinames));
end;

if ~parms.network_flag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.networks = parms.networks;
if isempty(results.networks)
  results.networks = abcd_define_rsBOLD_networks(parms.roi_type);
end;
results.network_names = fieldnames(results.networks);
results.nn = length(results.network_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select values for ROIs in roinames
results.roinames = results.targ_roinames;
results.nrois = length(results.roinames);
if ~isempty(parms.roinames)
  [tmp,ind_parms,ind_results] = intersect(parms.roinames,results.roinames);
  if isempty(ind_results)
    error('input roinames are not valid');
  end;
  % maintain order in roinames
  [tmp,ind_sort] = sort(ind_parms);
  ind_results = ind_results(ind_sort);
  results.nrois = length(ind_results);
%  fprintf('%s: %d ROIs selected\n',mfilename,results.nrois);
  % keep results only for selected ROIs
  results.roinames = results.roinames(ind_results);
  R = R(ind_results,:,:);
end;

if parms.save_corr_flag
  results.R = R;
  results.Z = atanh(results.R);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate mean correlation for each network
switch parms.network_flag
  case 1
    results.nt = results.nrois;
  case 2
    results.nt = 1;
  case 3
    results.nt = results.nn;
end;
results.network_corr = zeros(results.nn,results.nt,results.nsubs);
for n=1:results.nn
  network_name = results.network_names{n};
  network_seed_roinames = results.networks.(network_name);
  ns = length(network_seed_roinames);
  switch parms.network_flag
    case 1
      targ_names = results.roinames;
    case 2
      targ_names = {network_name};
    case 3
      targ_names = results.network_names;
  end;
  for t=1:results.nt
    targ_name = targ_names{t};
    if parms.network_flag>1
      % get list of ROIs to average, excluding seed
      targ_roinames = results.networks.(targ_name);
    else
      targ_roinames = {targ_name};
    end;
    % find all ROIs with this name (e.g. lh and rh)
    ind_rois = [];
    for r=1:length(targ_roinames)
      if parms.strict_roinames_flag
        ind_rois = [ind_rois,...
          find(strcmp(results.roinames,targ_roinames{r}))];
      else
        ind_rois = [ind_rois,...
          find(~cellfun(@isempty,regexp(results.roinames,targ_roinames{r})))];
      end;
    end;
    nrois = length(ind_rois);
    for r=1:nrois
      ind_roi = ind_rois(r);
      roiname = results.roinames{ind_roi};
      % exclude seed ROI from targets
      if parms.strict_roinames_flag
        ind_tmp = find(~strcmp(roiname,network_seed_roinames));
      else
        ind_tmp = find(cellfun(@isempty,regexp(roiname,network_seed_roinames)));
      end;
      seed_roinames = network_seed_roinames(ind_tmp);
      [tmp,ind_seeds] = intersect(results.seed_roinames,seed_roinames);
      if isempty(ind_seeds)
        fprintf('%s: WARNING: ind_seeds is empty for %s and %s\n',...
          mfilename,network_name,roiname);
        continue;
      end;
      % sum Z-transformed correlations for this ROI
      tmp_Z = mean(atanh(R(ind_roi,ind_seeds,:)),2);
      results.network_corr(n,t,:) = results.network_corr(n,t,:) + tmp_Z;
    end;
    % divide by the number of ROIs for this network
    results.network_corr(n,t,:) = results.network_corr(n,t,:) / nrois;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reshape corr matrix
results.network_corr_array = num2cell(reshape(permute(results.network_corr,[2,1,3]),...
                            [results.nn*results.nt,results.nsubs])');

% create column labels
col_labels = cell(1,results.nn*results.nt);
k = 1;
for n=1:results.nn
  network_name = results.network_names{n};
  switch parms.network_flag
    case 1
      targ_names = results.roinames;
    case 2
      targ_names = {network_name};
    case 3
      targ_names = results.network_names;
  end;
  for t=1:results.nt
    targ_name = targ_names{t};
    col_labels{k} = [network_name '->' targ_name];
    k = k + 1;
  end;
end;
results.network_corr_labels = col_labels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

