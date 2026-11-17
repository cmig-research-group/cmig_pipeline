function [parms,errcode,errmsg] = cmig_check_bids_func(parms,events_flag)
%function [parms,errcode,errmsg] = cmig_check_bids_func(parms,[events_flag])
%
% Purpose: locate and validate fMRI input files (bold data, motion files,
%          and optionally events files) for a given subject/session
%          from the cmig_proc_bids_func outputs: _desc-preproc_{series_type}
%          data, _desc-motion_timeseries.tsv motion estimates, and the
%          _events.tsv copies, which keep their rawdata names
%
% Required Input:
%   parms: parameter struct with fields:
%     participant_id, session_id, proc_rootdir, modality, series_type, task
%
% Optional Input:
%   events_flag: [0|1] whether to locate and check events files
%     {default = 0}
%
% Output:
%   parms: updated parameter struct with fields added:
%     fnames_data, fnames_info, scan_info, fnames_motion
%     fnames_events (if events_flag=1)
%   errcode: [0|1] error code (1 = fatal error)
%   errmsg: error message string (empty if no error)
%
% Created:  05/21/2026 by Don Hagler
% Prev Mod: 05/21/2026 by Don Hagler
% Prev Mod: 11/06/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

errcode = 0; errmsg = [];

if ~exist('events_flag','var') || isempty(events_flag), events_flag = 0; end

% get input file names
[parms.fnames_data,parms.fnames_info,parms.scan_info] = cmig_get_bids_data(parms.participant_id,parms.session_id,...
  parms.proc_rootdir,parms.modality,parms.series_type,parms.task);

% keep the proc outputs only (desc-preproc), not a rawdata-named file left
%   in the proc tree
proc_suffix = sprintf('_desc-preproc_%s.nii.gz',parms.series_type);
if ~isempty(parms.fnames_data)
  idx = find(~cellfun(@isempty,regexp(parms.fnames_data,[regexptranslate('escape',proc_suffix) '$'],'once')));
  if ~isempty(parms.scan_info) && length(parms.scan_info)==length(parms.fnames_data)
    parms.scan_info = parms.scan_info(idx);
  end
  parms.fnames_data = parms.fnames_data(idx);
  parms.fnames_info = parms.fnames_info(idx);
end

% check that there are input files to analyze
if isempty(parms.fnames_data)
  errcode = 1;
  errmsg = 'not processed: no usable scans to analyze';
  return;
end

% report input files found
fprintf('%s: input fMRI data:\n',mfilename);
for f = 1:length(parms.fnames_data)
  fprintf(' %s\n',parms.fnames_data{f});
end

% get fnames_motion based on fnames_data
parms.fnames_motion = strrep(parms.fnames_data,proc_suffix,'_desc-motion_timeseries.tsv');

% check that motion files exist
for f = 1:length(parms.fnames_motion)
  if ~exist(parms.fnames_motion{f},'file')
    [~,fn,fx] = fileparts(parms.fnames_motion{f});
    errcode = 1;
    errmsg = sprintf('motion file not found: %s',[fn fx]);
    return;
  end
end

if ~events_flag, return; end

% get fnames_events based on fnames_data
parms.fnames_events = strrep(parms.fnames_data,proc_suffix,'_events.tsv');

% check that events files exist; warn if missing
for f = 1:length(parms.fnames_events)
  if ~exist(parms.fnames_events{f},'file')
    [~,fn,fx] = fileparts(parms.fnames_events{f});
    fprintf('%s: WARNING: events file not found: %s\n',mfilename,[fn fx]);
    parms.fnames_events{f} = [];
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
