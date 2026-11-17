function cmig_calc_fmri_network_corr(varargin)
%function cmig_calc_fmri_network_corr([options])
%
% Purpose: concatenate within- and between- network resting-state fMRI correlations
%          calculated from surface-based and volumetric ROI correlations
%
% Optional Parameters ('key',value pairs):
%   'indir': input directory
%     {default = [pwd '/input']}
%   'instem': input file fstem of concatenated time series data
%     {default = 'corrmat'}
%   'infostem': input file stem of info file
%     {default = 'vol_info'}
%   'outdir': output directory
%     {default = [pwd '/output']}
%   'outstem': output file stem of concatenated correlation matrices
%     {default = 'network_corrmat'}
%   'network_flag': [1|2|3] average correlation within networks of ROIs
%     1: mean correlation between each network and each ROI
%     2: mean correlation within each network
%     3: mean correlation within and between each network
%     {default = 3}
%   'roinames': cell array of ROIs to be included in output
%     applies if network_flag = 1
%     if empty, use all ROIs
%     {default = []}
%   'forceflag': [0|1] overwrite existing output files
%     {default = 0}
%
% Created:  06/22/2026 by Don Hagler
% Prev Mod: 06/22/2026 by Don Hagler
% Last Mod: 11/11/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on ABCD_concat_roi_rsfMRI_network_corr
%   Created:  05/31/2024 by Don Hagler
%   Last Mod: 04/23/2026 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms(varargin,{...
  'indir',[pwd '/input'],[],...
  'instem','corrmat',[],...
  'infostem','vol_info',[],...
  'outdir',[pwd '/output'],[],...
  'outstem','network_corrmat',[],...
  'network_flag',3,[],... % [0:3],...
  'roinames',[],[],...
  'forceflag',false,[],... % [false true],...
...
  'networks',[],[],...
  'roi_type','gordon',[],... % {'gordon','aparc'},...
  'save_corr_flag',false,[],... % [false true],...
  'netcorr_tags',{'targ_roinames','network_flag', 'networks',...
                  'roi_type','roinames','strict_roinames_flag','save_corr_flag'},[],...
});

% parse string inputs
%   as number, string, cell array, or numeric array
if isdeployed
  tags = fieldnames(parms);
  for t=1:length(tags)
    parms.(tags{t}) = cmig_parseval(parms.(tags{t}));
  end
end

% check input file
if ~exist(parms.indir,'dir')
  error('input directory %s not found',parms.indir);
end
fname_in = sprintf('%s/%s.mat',parms.indir,parms.instem);
if ~exist(fname_in,'file')
  error('input file %s not found',fname_in);
end

% define networks
if isempty(parms.networks)
  parms.networks = abcd_define_rsBOLD_networks(parms.roi_type);
end
network_names = fieldnames(parms.networks)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = sprintf('%s/%s.mat',parms.outdir,parms.outstem);
if ~exist(fname_out,'file') || parms.forceflag

  % load concatenated rsfMRI ROI correlations
  fprintf('%s: loading %s...\n',mfilename,fname_in);
  tic
  data = load(fname_in);
  toc

  %% NOTE: fields in data struct:
  %  'corrmat','varmat','ntpointvec',...
  %  'nvolvec','nvolfdsubvec','nvolfdsubcontvec',...
  %  'meanfdvec','meanfdmat','meanfdcvec','meanfdcmat',...
  %  'meanrotvec','meantransvec','maxrotvec','maxtransvec','maxfdvec',...
  %  'meanrotcvec','meantranscvec','maxrotcvec','maxtranscvec','maxfdcvec',...
  %  'nroi','nnodes','roinames','roi1vec','roi2vec','ndirs'

  % save value for later
  roinames = data.roinames;
  nroi = data.nroi;
  ndirs = data.ndirs;
  ntpointvec = data.ntpointvec;
  nvolvec = data.nvolvec;
  nvolfdsubvec = data.nvolfdsubvec;
  nvolfdsubcontvec = data.nvolfdsubcontvec;
  meanfdvec = data.meanfdvec;
  meanfdmat = data.meanfdmat;
  meanfdcvec = data.meanfdcvec;
  meanfdcmat = data.meanfdcmat;
  meanrotvec = data.meanrotvec;
  meantransvec = data.meantransvec;
  maxrotvec = data.maxrotvec;
  maxtransvec = data.maxtransvec;
  maxfdvec = data.maxfdvec;
  meanrotcvec = data.meanrotcvec;
  meantranscvec = data.meantranscvec;
  maxrotcvec = data.maxrotcvec;
  maxtranscvec = data.maxtranscvec;
  maxfdcvec = data.maxfdcvec;
  nnet = length(network_names);

  % restore full correlation matrix
  fprintf('%s: restoring full correlation matrix...\n',mfilename);
  tic
  idx_triu = find(triu(ones(nroi,nroi)));
  idx_tril = find(tril(ones(nroi,nroi)));
  R = nan(nroi*nroi,ndirs,'single');
  R(idx_triu,:) = data.corrmat';
  R_lower = reshape(permute(reshape(R,[nroi,nroi,ndirs]),[2,1,3]),[nroi*nroi,ndirs]);
  R(idx_tril,:) = R_lower(idx_tril,:);
  clear R_lower;
  R = reshape(R,[nroi,nroi,ndirs]);
  toc;

  % calculate network correlations
  fprintf('%s: calculating network correlations...\n',mfilename);
  tic
  args = mmil_parms2args(parms,parms.netcorr_tags);
  tmp = abcd_network_corr(R,roinames,args{:});
  toc

  % reformat output depending on network_flag
  switch parms.network_flag
    case 1 % network to ROI
      if ~isempty(parms.roinames)
        nroi = length(parms.roinames);
        roinames = parms.roinames;
      end
      corrmat = reshape(single(tmp.network_corr),[nnet*nroi,ndirs])';
      nnodes = size(corrmat,2);
      net1mat = repmat([1:nnet]',[1,nroi]);
      roi2mat = repmat([1:nroi],[nnet,1]);
      net1vec = reshape(net1mat,[nnet*nroi,1]);
      roi2vec = reshape(roi2mat,[nnet*nroi,1]);
    case 2 % within network
      corrmat = reshape(single(tmp.network_corr),[nnet,ndirs])';
      nnodes = size(corrmat,2);
      net1vec = [1:nnet];
      net2vec = [1:nnet];
    case 3 % network to network (between and within)
      corrmat = reshape(single(tmp.network_corr),[nnet*nnet,ndirs])';
      % upper triangle of nnet x nnet correlation matrix
      idx_triu = find(triu(ones(nnet,nnet)));
      nnodes = length(idx_triu);
      corrmat = corrmat(:,idx_triu);
      % network indices for each node
      net1mat = repmat([1:nnet]',[1,nnet]);
      net2mat = repmat([1:nnet],[nnet,1]);
      net1vec = net1mat(idx_triu);
      net2vec = net2mat(idx_triu);
  end

  % save output file
  fprintf('%s: saving results to %s...\n',mfilename,fname_out);
  if parms.network_flag==1
    save(fname_out,'corrmat','ntpointvec',...
                   'nvolvec','nvolfdsubvec','nvolfdsubcontvec',...
                   'meanfdvec','meanfdmat','meanfdcvec','meanfdcmat',...
                   'meanrotvec','meantransvec','maxrotvec','maxtransvec','maxfdvec',...
                   'meanrotcvec','meantranscvec','maxrotcvec','maxtranscvec','maxfdcvec',...
                   'nnet','network_names','net1vec',...
                   'nroi','roinames','roi2vec',...
                   'nnodes','ndirs',...
                   '-v7');
  else
    save(fname_out,'corrmat','ntpointvec',...
                   'nvolvec','nvolfdsubvec','nvolfdsubcontvec',...
                   'meanfdvec','meanfdmat','meanfdcvec','meanfdcmat',...
                   'meanrotvec','meantransvec','maxrotvec','maxtransvec','maxfdvec',...
                   'meanrotcvec','meantranscvec','maxrotcvec','maxtranscvec','maxfdcvec',...
                   'nnet','network_names','net1vec','net2vec',...
                   'nnodes','ndirs',...
                   '-v7');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

