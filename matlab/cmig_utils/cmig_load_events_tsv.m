function [onset_mat,duration_mat,stim_types] = cmig_load_events_tsv(fname_events,ntpoints,TR,stim_labels,onset_offset)
%function [onset_mat,duration_mat,stim_types] = cmig_load_events_tsv(fname_events,ntpoints,TR,stim_labels,[onset_offset])
%
% Purpose: load an events TSV and return per-TR onset and duration indicator
%   matrices, so that the concatenated time series carry the design
%
% Required Input:
%   fname_events: full path of the events TSV; an empty or missing file gives
%     all-zero matrices sized from stim_labels
%   ntpoints: number of time points in the run
%   TR: repetition time in seconds
%   stim_labels: cell array of condition names, fixing the column order; if
%     empty, the order is taken from the events file, first appearance first
%
% Optional Input:
%   onset_offset: seconds removed from the start of the run, since the TR grid
%     starts at the first retained frame while the events TSV is in pre-skip
%     time
%     {default = 0}
%
% Output:
%   onset_mat: [ntpoints,nstim] single, 1 at the TR an event onset falls in
%   duration_mat: [ntpoints,nstim] single, 1 at each TR the event overlaps by
%     more than half a TR at each end
%   stim_types: cell array of the condition names the columns are in
%
% Created:  09/17/2026 by Don Hagler
% Last Mod: 09/17/2026 by Don Hagler
%

onset_mat = []; duration_mat = []; stim_types = {};
if ~mmil_check_nargs(nargin,4), return; end
if ~exist('onset_offset','var') || isempty(onset_offset), onset_offset = 0; end

if isempty(fname_events) || ~exist(fname_events,'file')
  if ~isempty(fname_events)
    fprintf('%s: WARNING: events file not found: %s\n',mfilename,fname_events);
  end
  if ~isempty(stim_labels)
    stim_types   = stim_labels;
    onset_mat    = zeros(ntpoints,length(stim_labels),'single');
    duration_mat = zeros(ntpoints,length(stim_labels),'single');
  else
    onset_mat    = zeros(ntpoints,0,'single');
    duration_mat = zeros(ntpoints,0,'single');
  end
  return;
end
events = mmil_csv2struct(fname_events,'separator','\t');
events = alias_trial_types(events);
all_types = {events.trial_type};
if ~isempty(stim_labels)
  stim_types = stim_labels;
else
  stim_types = unique(all_types(~strcmp(all_types,'dummy')),'stable');
end
nstim = length(stim_types);
onset_mat    = zeros(ntpoints,nstim,'single');
duration_mat = zeros(ntpoints,nstim,'single');
if isempty(TR) || ntpoints == 0, return; end
% per-TR time grid (matching abcd_load_tfMRI_data convention)
TR_offset = linspace(TR,ntpoints*TR,ntpoints);  % end time of each TR
TR_onset  = TR_offset - TR;                     % start time of each TR
mindur = 0.5 * TR;
for c = 1:nstim
  stim_type = stim_types{c};
  idx_events = find(strcmp(all_types,stim_type));
  if isempty(idx_events), continue; end
  onset_vals = {events(idx_events).onset};
  if isnumeric(onset_vals{1})
    event_times = cell2mat(onset_vals);
  else
    event_times = str2double(onset_vals);
  end
  dur_vals = {events(idx_events).duration};
  if isnumeric(dur_vals{1})
    event_durs = cell2mat(dur_vals);
  else
    event_durs = str2double(dur_vals);
  end
  % the TR grid starts at the first retained frame while the events tsv is
  %   in pre-skip time, so shift every onset back by the frames removed
  event_times = event_times - onset_offset;
  i_keep = find(event_times >= 0);
  event_times = event_times(i_keep);
  event_durs  = event_durs(i_keep);
  stim_onset_vec    = zeros(ntpoints,1,'single');
  stim_duration_vec = zeros(ntpoints,1,'single');
  % onset: TR frame where event onset falls
  event_frames = ceil((event_times + eps) / TR);
  event_frames = event_frames(ismember(event_frames,1:ntpoints));
  stim_onset_vec(event_frames) = 1;
  % duration: all TRs that overlap event (by more than mindur at each end)
  for k = 1:length(event_times)
    ind_rep = find(event_times(k) < (TR_offset - mindur) & ...
                   (event_times(k) + event_durs(k)) > (TR_onset + mindur));
    stim_duration_vec(ind_rep) = 1;
  end
  onset_mat(:,c)    = stim_onset_vec;
  duration_mat(:,c) = stim_duration_vec;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function events = alias_trial_types(events)
  % events files written before the 08/14/2026 rename carry ssd_stop where
  %   abcd_set_contrasts_sst now says early_stop
  idx = find(strcmp({events.trial_type},'ssd_stop'));
  if ~isempty(idx)
    [events(idx).trial_type] = deal('early_stop');
  end
return;
