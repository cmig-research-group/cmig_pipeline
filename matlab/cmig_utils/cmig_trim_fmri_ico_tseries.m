function cmig_trim_fmri_ico_tseries(varargin)
%function cmig_trim_fmri_ico_tseries([options])
%
% Purpose: trim preprocessed icosahedral surface fMRI time series
%          to a fixed number of valid (uncensored) frames per visit
%
% Optional Parameters ('key',value pairs):
%   'indir': input directory
%     {default = [pwd '/input']}
%   'instem': input file stem of preprocessed ico time series
%     {default = 'ico_tseries_prep'}
%   'outdir': output directory
%     {default = [pwd '/output']}
%   'outstem': output file stem
%     {default = 'ico_tseries_trim'}
%   'trim_ntpoints': number of valid frames to retain
%     if Inf, function returns without output
%     {default = Inf}
%   'trim_rand_flag': [0|1] randomly select trim_ntpoints frames
%     if 0, use first N frames after censoring
%     {default = false}
%   'trim_seed': seed of the random number generator for trim_rand_flag=1
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Output has shape [1, ndirs, trim_ntpoints, nverts] with all frames
%   marked as uncensored (censvec = 0). Visits with fewer than
%   trim_ntpoints valid frames are left as NaN.
%
% Created:  06/22/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Prev Mod: 07/02/2026 by Don Hagler
% Prev Mod: 10/13/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',[pwd '/input'],[],...
  'instem','ico_tseries_prep',[],...
  'outdir',[pwd '/output'],[],...
  'outstem','ico_tseries_trim',[],...
  'trim_ntpoints',Inf,[],...
  'trim_rand_flag',false,[false true],...
  'trim_seed',0,[0,2^32-1],...
  'forceflag',false,[false true],...
});

% parse string inputs as number, string, cell array, or numeric array
if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

if ~isfinite(parms.trim_ntpoints)
  fprintf('%s: trim_ntpoints is Inf, nothing to do\n',mfilename);
  return;
end

% check input file
if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end
fname_in = sprintf('%s/%s.mat',parms.indir,parms.instem);
if ~exist(fname_in,'file')
  error('input file %s not found',fname_in);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = sprintf('%s/%s.mat',parms.outdir,parms.outstem);
if ~exist(fname_out,'file') || parms.forceflag

  fprintf('%s: loading %s...\n',mfilename,fname_in);
  tic; data = load(fname_in); toc;

  nverts    = data.nverts;
  ndirs     = data.ndirs;
  trim_n    = parms.trim_ntpoints;
  nmotparms = size(data.datamat_motion,4);

  % build output struct: copy all scalar/metadata fields, replace array fields
  out = struct();
  skip_fields = {'datamat_tsdata','datamat_motion','censvec','nruns','ntpoints'};
  fns = fieldnames(data);
  for f=1:length(fns)
    fn = fns{f};
    if ~ismember(fn,skip_fields)
      out.(fn) = data.(fn);
    end
  end
  out.nruns    = 1;
  out.ntpoints = trim_n;
  % NaN-initialize; censvec=1 means censored (visit has no data by default)
  out.datamat_tsdata = NaN(1,ndirs,trim_n,nverts,'single');
  out.datamat_motion = NaN(1,ndirs,trim_n,nmotparms,'single');
  out.censvec        = ones(1,ndirs,trim_n,'single');

  % seed rng once before the visit loop (only relevant when trim_rand_flag=1)
  if parms.trim_rand_flag
    rng(parms.trim_seed,'twister');
  end

  % which vertices have no NaN timepoints: [nruns,ndirs,1,nverts]
  defmat = isfinite(sum(data.datamat_tsdata,3));

  for d=1:ndirs

    if d==1 || mod(d,100)==0
      fprintf('%s: diri = %d of %d (now=%s)\n',...
        mfilename,d,ndirs,datestr(now,'HH:MM:SS'));
    end

    % skip visit if no valid data at all
    if ~nnz(mmil_rowvec(defmat(:,d,:,:)))
      continue;
    end

    % concatenate uncensored frames across runs
    D = [];
    M = [];
    for s=1:data.nruns
      if ~nnz(squeeze(defmat(s,d,:,:))), continue; end
      tm = squeeze(data.datamat_tsdata(s,d,:,:));  % [ntpoints, nverts]
      mm = squeeze(data.datamat_motion(s,d,:,:));  % [ntpoints, nmotparms]
      cv = squeeze(data.censvec(s,d,:));            % [ntpoints, 1]
      idx_valid = find(~cv);
      D = cat(1,D,tm(idx_valid,:));
      M = cat(1,M,mm(idx_valid,:));
    end
    ntpoints = size(D,1);

    if ~ntpoints
      continue;
    end

    if ntpoints < trim_n
      fprintf('%s: WARNING: %d/%d : skipping because fewer than %d valid data points...\n',...
        mfilename,d,ndirs,trim_n);
      continue;
    end

    % select trim_n frames
    if parms.trim_rand_flag
      idx_trim = randperm(ntpoints);
    else
      idx_trim = 1:ntpoints;
    end
    D = D(idx_trim(1:trim_n),:);
    M = M(idx_trim(1:trim_n),:);

    % store selected frames; mark all trim_n frames as uncensored
    out.datamat_tsdata(1,d,:,:) = reshape(D,[1,1,trim_n,nverts]);
    out.datamat_motion(1,d,:,:) = reshape(M,[1,1,trim_n,nmotparms]);
    out.censvec(1,d,:)          = 0;

  end

  fprintf('%s: saving results to %s...\n',mfilename,fname_out);
  save(fname_out,'-struct','out','-v7');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
