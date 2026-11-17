function [A,Ainv,V,VV,F] = rsi_create_matrix(qmat,varargin)
%function [A,Ainv,V,VV,F] = rsi_create_matrix(qmat,[options])
%
% Required Parameters:
%   qmat: matrix of diffusion direction vectors
%
% Optional Parameters:
%   'bvals': vector of b values
%     one for all, or one for each diffusion direction
%     {default = 1000}
%   'lambda': regularization constant
%     {default = 0.1}
%   'ADC_long': longitudinal ADC
%     {default = 1e-3}
%   'ADC_trans_min': minimum transverse ADC
%     {default = 0}
%   'ADC_trans_max': maximum transverse ADC
%     {default = 0.9e-3}
%   'num_ADC_trans': number of transverse ADC size scales
%     {default = 5}
%   'SH_order': spherical harmonic order -- must be even
%     {default = 4}
%   'iso_restricted_flag': [0|1] model isotropic diffusion of restricted water
%     {default = 1}
%   'iso_hindered_flag': [0|1] model isotropic diffusion of hindered water
%     {default = 1}
%   'iso_free_flag': [0|1] model isotropic diffusion of free water
%     {default = 1}
%   'ADC_hindered': ADC of isotropic hindered water (e.g. edema)
%     {default = 1.5e-3}
%   'ADC_free': apparent diffusion coefficient (ADC) of
%               isotropic free water (e.g. CSF)
%     {default = 3e-3}
%   'ADC_iso_min': minimum isotropic ADC
%     {default = 0}
%   'ADC_iso_max': maximum isotropic ADC
%     {default = 3e-3}
%   'num_ADC_iso': number of isotropic ADC size scales
%     {default = 0}
%   'ADC_iso_vals': vector of isotropic ADC values
%     if empty, will be set according to
%       ADC_iso_min, ADC_iso_max, and num_ADC_iso
%     {default = []}
%
% Created:  09/12/2019 DH
% Last Mod: 09/12/2019 DH
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = []; Ainv = []; V = []; VV = []; F = [];
if ~mmil_check_nargs(nargin,2), return; end;

parms = check_input(qmat,varargin);

VV = icosaSphere(3);
V = VV.vertices;
Q = parms.qmat.*repmat(sqrt(parms.bvals),1,3);
F = rsi_SH_matrix(V,parms.SH_order);
F0 = rsi_SH_matrix(V,0);
A = [];

% series of FODs for varying size scales
for t=1:parms.num_ADC_trans
  R = rsi_FOD_matrix(Q,V,parms.ADC_long,parms.ADC_trans_vals(t));
  A = [A R*F];
end

% series of isotropic FODs for varying size scales
for t=1:parms.num_ADC_iso
  R = rsi_FOD_matrix(Q,V,parms.ADC_iso_vals(t),parms.ADC_iso_vals(t));
  A = [A R*F0];
end

% isotropic restricted
if parms.iso_restricted_flag
  R0 = rsi_FOD_matrix(Q,V,0,0);
  A = [A R0*F0];
end;

% isotropic hindered water fractions
if parms.iso_hindered_flag
  R_hindered = rsi_FOD_matrix(Q,V,parms.ADC_hindered,parms.ADC_hindered);
  A = [A R_hindered*F0];
end;

% isotropic free and hindered water fractions
if parms.iso_free_flag
  R_free = rsi_FOD_matrix(Q,V,parms.ADC_free,parms.ADC_free);
  A = [A R_free*F0];
end;

% compute regularized inverse
AtA = A'*A;
Ainv = inv(AtA+parms.lambda*mean(diag(AtA))*eye(size(AtA)))*A';

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(qmat,options)
  parms = mmil_args2parms(options,{...
    'qmat',qmat,[],...
    ...
    'bvals',1000,[],...
    'lambda',0.1,[],...
    'ADC_trans_min',0,[],...
    'ADC_trans_max',0.9e-3,[],...
    'num_ADC_trans',5,[],...
    'SH_order',4,[2:2:10],...
    'iso_free_flag',true,[false true],...
    'iso_hindered_flag',true,[false true],...
    'iso_restricted_flag',true,[false true],...
    'ADC_free',3e-3,[],...
    'ADC_hindered',1.5e-3,[],...
    'ADC_long',1e-3,[],...
    'ADC_iso_min',0,[],...
    'ADC_iso_max',3e-3,[],...
    'num_ADC_iso',0,[],...
    'ADC_iso_vals',[],[],...
  });
  
  % check qmat is correct size
  if size(parms.qmat,2)~=3
    error('qmat must have 3 columns');
  end;
  % ensure that qmat contains unit vectors
  qlength = max(sqrt(sum(parms.qmat.^2,2)),eps);
  parms.qmat_orig = parms.qmat;
  parms.qmat = parms.qmat ./ repmat(qlength,1,3);
  % number of frames according to qmat
  parms.nf = size(parms.qmat,1);
  % extend bvals to nf elements if only 1 supplied
  if numel(parms.bvals)==1
    parms.bvals = parms.bvals*ones(parms.nf,1);
  elseif numel(parms.bvals)==parms.nf
    parms.bvals = reshape(parms.bvals,[parms.nf,1]);
  else
    error('number of bvals must match frames in qmat');
  end;
  % set ADC iso vals
  if isempty(parms.ADC_iso_vals) && parms.num_ADC_iso>0
    parms.ADC_iso_vals = ...
      linspace(parms.ADC_iso_min,parms.ADC_iso_max,parms.num_ADC_iso);
  else
    parms.num_ADC_iso = length(parms.ADC_iso_vals);
  end;
  % set ADC trans vals
  parms.ADC_trans_vals =...
    linspace(parms.ADC_trans_min,parms.ADC_trans_max,parms.num_ADC_trans);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


