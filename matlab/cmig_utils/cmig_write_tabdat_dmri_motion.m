function cmig_write_tabdat_dmri_motion(varargin)
%function cmig_write_tabdat_dmri_motion([options])
%
% Purpose: write tabulated motion QC data for dMRI (session-mean motion)
%   Reads per-scan *_desc-motion_timeseries.tsv files (written by the dMRI
%   proc pipeline), or _motion.tsv beside a plain _dwi data name (released
%   mproc derivatives), and writes session-average mot/rotat/trans to a TSV.
%   A missing motion file logs a WARNING and that scan is excluded from the
%   session average, rather than aborting.
%
% Optional Input ('key',value):
%   'fnames_data': cell array of DWI data file names (nii/nii.gz)
%     {default = {}}
%   'series_type': series type used to derive the motion tsv name
%     with or without the desc-preproc_ prefix
%     {default = 'dwi'}
%   'motion_radius': head radius (mm) for rotation-to-distance conversion
%     {default = 50}
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem ({outstem}_motionqc.tsv)
%     {default = ''}
%   'colstem': prefix attached to column names
%     {default = 'mr_y_qc__mot__dmri'}
%   'participant_id': participant ID string
%     {default = []}
%   'session_id': session ID string
%     {default = []}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  07/22/2026 by Don Hagler
% Prev Mod: 10/02/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'fnames_data',{},[],...
  'series_type','dwi',[],...
  'motion_radius',50,[1,100],...
  'outdir',pwd,[],...
  'outstem','',[],...
  'colstem','mr_y_qc__mot__dmri',[],...
  'participant_id',[],[],...
  'session_id',[],[],...
  'forceflag',false,[false true],...
});

if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

fprintf('%s: writing motion QC tabulated data...\n',mfilename);

fname_tsv = sprintf('%s/%s_motionqc.tsv',parms.outdir,parms.outstem);
if exist(fname_tsv,'file') && ~parms.forceflag, return; end

mmil_mkdir(parms.outdir);

% series_type given with or without desc-preproc
suffix = regexprep(parms.series_type,'^desc-preproc_','');
pat = sprintf('(_desc-preproc)?_%s\\.nii(\\.gz)?$',suffix);

nscans = length(parms.fnames_data);
mot_sum = 0; rotat_sum = 0; trans_sum = 0; n = 0;
for s = 1:nscans
  fname_data = parms.fnames_data{s};
  if isempty(regexp(fname_data,pat,'once'))
    fprintf('%s: WARNING: data name does not end in _[desc-preproc_]%s.nii[.gz] (scan excluded from session average): %s\n',...
      mfilename,suffix,fname_data);
    continue;
  end
  fname_motion = regexprep(fname_data,pat,'_desc-motion_timeseries.tsv');
  % plain _dwi name (released mproc derivatives): motion file is _motion.tsv
  if ~exist(fname_motion,'file') && isempty(regexp(fname_data,'_desc-preproc_','once'))
    fname_motion_old = regexprep(fname_data,pat,'_motion.tsv');
    if exist(fname_motion_old,'file'), fname_motion = fname_motion_old; end
  end
  if ~exist(fname_motion,'file')
    fprintf('%s: WARNING: motion timeseries file not found (scan excluded from session average): %s\n',...
      mfilename,fname_motion);
    continue;
  end
  motion_tseries = mmil_load_motion_tsv(fname_motion);
  if isempty(motion_tseries)
    fprintf('%s: WARNING: empty motion timeseries: %s\n',mfilename,fname_motion);
    continue;
  end
  [mean_motion,mean_trans,mean_rot] = mmil_mean_motion(motion_tseries,parms.motion_radius);
  mot_sum = mot_sum + mean_motion;
  rotat_sum = rotat_sum + mean_rot;
  trans_sum = trans_sum + mean_trans;
  n = n + 1;
end

if n >= 1
  avgdata = [mot_sum/n, rotat_sum/n, trans_sum/n];
else
  avgdata = [NaN, NaN, NaN];
end

colnames = {sprintf('%s__mot_mean',parms.colstem),...
            sprintf('%s__rotat_mean',parms.colstem),...
            sprintf('%s__trans_mean',parms.colstem)};
header = cat(2,{'participant_id','session_id'},colnames);
datarow = cat(2,{parms.participant_id,parms.session_id},num2cell(avgdata));
datarow(cat(2,false(1,2),isnan(avgdata))) = {'n/a'};
all_data = cat(1,header,datarow);
mmil_write_csv(fname_tsv,all_data,'separator','\t','quote_flag',0);

return;
