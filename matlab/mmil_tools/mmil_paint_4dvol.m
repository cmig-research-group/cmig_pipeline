function [surfmat,normvecs,surf] = mmil_paint_4dvol(subj,vol,M,varargin);
%function [surfmat,normvecs,surf] = mmil_paint_4dvol(subj,vol,M,[options]);
%
% Purpose: resample multi-frame volume data to cortical surface
%
% Usage:
%  [surfmat,normvecs,surf] = mmil_paint_4dvol(subj,vol,M,'key1', value1,...);
%
% Required parameters:
%  subj: FreeSurfer subject name
%  vol: multiframe (4D) input volume
%  M: vox2ras matrix for input volume
%
% Optional parameters:
%  'hemi': cortical hemisphere ('lh' or 'rh')
%    {default = 'lh'}
%  'regmat': registration matrix
%    if empty, assumes input volume is in register with FreeSurfer volumes
%    {default = []}
%  'pmin': minimum distance (mm) to project along surface vertex normal
%    {default = 1}
%  'pmax': maximum distance (mm) to project along surface vertex normal
%    {default = 1}
%  'pstep': step distance (mm) along surface vertex normal
%    {default = 1}
%  'interpmethod': interpolation method ('nearest','linear','spline','cubic')
%    {default = 'linear'}
%  'subjdir': subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'surfname': surface to paint onto
%    {default = white}
%  'surfname_outer': surface to check orientation of normal vectors
%    {default = 'pial'}
%  'surf': surface structure
%    if supplied, will skip loading from file
%    {default = []}
%  'normvecs': matrix of normal vectors with size = [nverts,3]
%    if supplied, will skip calculating surf
%    {default = []}
%
% Output:
%   surfmat: surface data matrix with size = [nverts,nframes,nsteps]
%   normvecs: surface data matrix with size = [nverts,3]
%   surf: surface structure
%
% Created:  03/12/21 by Don Hagler
% Last Mod: 03/20/21 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surfmat = []; normvecs = []; surf = [];
if ~mmil_check_nargs(nargin,3), return; end;

% parse input parameters
parms = check_input(subj,vol,M,varargin);

% load surface file and compute normal vectors
[surf,normvecs] = load_surf(parms,vol);

% resample data from volume to surface
surfmat = resample_vol2surf(parms,vol,surf,normvecs);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(subj,vol,M,options)
  parms = mmil_args2parms(options,{...
    'subj',subj,[],...
    'M',M,[],...
  ...
    'hemi','lh',{'lh','rh'},...
    'regmat',eye(4),[],...
    'pmin',1,[-10,10],...
    'pmax',1,[-10,10],...
    'pstep',1,[0.001,10],...
    'nsmooth',0,[0,10000],...
    'interpmethod','linear',{'nearest','linear','spline','cubic'},...
    'subjdir',[],[],...
    'surfname','white',[],...
    'surfname_outer','pial',[],...
    'surf',[],[],...
    'normvecs',[],[],...
  ...
    'verbose',false,[false true],...
  });

  if isempty(parms.subjdir)
    parms.subjdir = getenv('SUBJECTS_DIR');
    if isempty(parms.subjdir)
      error('SUBJECTS_DIR not defined as an environment variable');
    end;
  else
    setenv('SUBJECTS_DIR',parms.subjdir);
  end;

  parms.fname_surf = sprintf('%s/%s/surf/%s.%s',...
    parms.subjdir,parms.subj,parms.hemi,parms.surfname);
  if ~exist(parms.fname_surf,'file')
    error('surface file %s not found',parms.fname_surf);
  end
  parms.fname_surf_outer = sprintf('%s/%s/surf/%s.%s',...
    parms.subjdir,parms.subj,parms.hemi,parms.surfname_outer);
  if ~exist(parms.fname_surf_outer,'file')
    error('surface file %s not found',parms.fname_surf_outer);
  end

  if any(size(parms.regmat)~=[4 4])
    error('regmat is wrong size');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [surf,normvecs] = load_surf(parms,vol)
  surf = []; normvecs = [];
  % for center coordinates
  ctx_vol = ctx_mgh2ctx(vol(:,:,:,1),parms.M);
  if isempty(parms.surf)
    % load surface file
    if parms.verbose
      fprintf('%s: loading surface file %s...\n',...
        mfilename,parms.fname_surf);
    end;
    surf = fs_read_surf(parms.fname_surf);
    surf = prep_surf(surf,ctx_vol.lphcent,parms.regmat);
  else
    surf = parms.surf;
  end;
  if isempty(parms.normvecs)
    % load outer surface file
    if parms.verbose
      fprintf('%s: loading surface file %s...\n',...
        mfilename,parms.fname_surf_outer);
    end;
    surf_outer = fs_read_surf(parms.fname_surf_outer);
    surf_outer = prep_surf(surf_outer,ctx_vol.lphcent,parms.regmat);
    % compute normal vectors
    normvecs = calc_normvecs(surf,surf_outer,parms.verbose);
  else
    normvecs = parms.normvecs;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function surf = prep_surf(surf,lphcent,M)
  ras2lph = [-1 -1 1];
  lph2ras = ras2lph;
  % add center coordinates
  surf.vertices = bsxfun(@plus,surf.vertices,lph2ras.*mmil_rowvec(lphcent));
  % convert from RAS to LPH
  surf.vertices = bsxfun(@times,surf.vertices,ras2lph);
  % transform surface coordinates with M
  V = cat(2,surf.vertices,ones(surf.nverts,1));
  V = (M * V')';
  surf.vertices = V(:,1:3);
  % preprocess surface into ctx format
  surf = preprocessQ(surf);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function normvecs = calc_normvecs(surf,surf_outer,verbose)
  % compute surface normal vectors
  normvecs = ts_calc_surf_normals(surf); %  normvecs = computeNormals(surf);
  % check normal vector orientation is pointing "out"
  outvecs = surf_outer.vertices - surf.vertices;
  outvecs = bsxfun(@rdivide,outvecs,sqrt(sum(outvecs.^2, 2)));
  dp = dot(normvecs,outvecs,2);
  % flip orientation of normal vectors if pointing in
  if median(dp(~isnan(dp))) < 0
    if verbose
      fprintf('%s: NOTE: flipping sign of normal vectors\n',mfilename);
    end;
    normvecs = -normvecs;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function surfmat = resample_vol2surf(parms,vol,surf,normvecs)
  % initialize output matrix
  surfmat = [];
  nverts = size(surf.vertices,1);
  nframes = size(vol,4);
  nsteps = length([parms.pmin:parms.pstep:parms.pmax]);
  surfmat = nan(nverts,nframes,nsteps);
  % sample each frame to surface
  if parms.verbose
    fprintf('%s: sampling values for %d frames and %d steps...\n',mfilename,nframes,nsteps);
  end;
  for f=1:nframes
    ctx_vol = ctx_mgh2ctx(vol(:,:,:,f),parms.M);
    % sample values along normal vectors onto surface
    vals = getValsOnRay(ctx_vol,surf,normvecs,...
                        parms.pmin,parms.pmax,parms.pstep,parms.nsmooth,parms.interpmethod);
    surfmat(:,f,:) = reshape(vals',[nverts,1,nsteps]);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

