function [M_af regStruct moving_af] = register_volumes_affine(moving, fixed, wvol, regtype, niter)

if ~exist('wvol','var') || ~isequal(size(wvol),size(fixed.imgs))
  wvol = true(size(fixed.imgs));
end

if ~exist('regtype','var') || isempty(regtype)
  regtype = 'icorr';
end

if ~exist('niter','var') || isempty(niter)
  niter = 1;
end

fixed_sub = subsample_volume(subsample_volume(fixed));
wvol_sub = subsample_volume(subsample_volume(wvol));
moving_sub = subsample_volume(subsample_volume(moving));

ivec_sub = find(wvol_sub>0.5);
sel = @(x)x(ivec_sub);
switch lower(regtype)
  case {'icorr'} % This appears to be broken
    sel = @(x,el)x(ivec_sub); 
    costfun = @(x)1-nancorr(sel(colvec(getfield(vol_resample_amd(moving_sub,fixed_sub,pars2M_af(x),1),'imgs'))),sel(colvec(fixed_sub.imgs))); 
  case {'mse'} % Should make this normalized by target variance
    costfun = @(x)nanmean(colvec((getfield(vol_resample_amd(moving_sub,fixed_sub,pars2M_af(x),1),'imgs')-fixed_sub.imgs).^2.*wvol_sub)); % Should perhaps use different cost function for binary volumes?
  case {'mi'}
    costfun = @(x)1-mutualinfo(sel(colvec(getfield(vol_resample_amd(moving_sub,fixed_sub,pars2M_af(x),1),'imgs'))),sel(colvec(fixed_sub.imgs))); % Should perhaps use different cost function for binary volumes?
  otherwise
    error(sprintf('%s: regtype=%s not implemented\n',mfilename,regtype));
end

parvec_init = [0 0 0 1 1 1 0 0 0 0 0 0]; % Should initialize translation based on center of mass
%disp(costfun(parvec_init))
%moving_init0 = vol_resample_amd(moving,fixed,pars2M_af(parvec_init),1);

parvec_init(1:3) = moving.Mvxl2lph(1:3,:)*([size(moving.imgs)/2 1])';
%disp(costfun(parvec_init))
%moving_init = vol_resample_amd(moving,fixed,pars2M_af(parvec_init),1);
%moving_init_sub = vol_resample_amd(moving_sub,fixed,pars2M_af(parvec_init),1);

M_tr = pars2M_af(parvec_init);
moving_tr = vol_resample_amd(moving,fixed,M_tr,1);

% reset random number generator so that results are reproducible
rng('default');

cost_opt = inf;
for iter = 1:niter % Random starts
  fprintf(1,'iter=%d/%d cost_opt=%f\n',iter,niter,cost_opt);
  parvec_opt_tmp = parvec_init;
  for phasei = 1:2
    if phasei==1 % Start with coarse rb-reg
      options.randscale = [3 3 3 0.0 0.0 0.0 5*pi/180 5*pi/180 5*pi/180 0 0 0];
      [parvec_opt_tmp cost_opt_tmp] = fminsearch_stochastic(costfun,parvec_opt_tmp,options);
%      disp(cost_opt_tmp)
      [M_af_tmp Mss Ms Mrl Mrp Mrh Mt] = pars2M_af(parvec_opt_tmp);
      M_rb_tmp = Mrl*Mrp*Mrh*Mt; 
      [parvec_opt_tmp cost_opt_tmp] = fminsearch(costfun,parvec_opt_tmp,options);
%      disp(cost_opt_tmp)
%      disp(parvec_opt_tmp);
      [M_af_tmp Mss Ms Mrl Mrp Mrh Mt] = pars2M_af(parvec_opt_tmp);
%      M_rb_tmp = Mrl*Mrp*Mrh*Mt; % Rigid-body component of affine -- not neccessarily least-squares optimal -- need to back out the scaling and shear component of affine;  M_af = Mss*Ms*Mrl*Mrp*Mrh*Mt*M_reg
    end
  end
  if cost_opt_tmp<cost_opt
    cost_opt = cost_opt_tmp;
    parvec_opt = parvec_opt_tmp;
    M_rb = M_rb_tmp;
    M_af = M_af_tmp;
  end
end

moving_rb = vol_resample_amd(moving,fixed,M_rb,1);
moving_af = vol_resample_amd(moving,fixed,M_af,1);

if 0 % isinteractive
  showVol(moving_tr,moving_rb,moving_af,fixed,wvol);  
  keyboard
end

regStruct.min_cost_af = cost_opt;
regStruct.parvec_af = parvec_opt;
regStruct.M_atl_to_vol_tr = M_tr;
regStruct.M_atl_to_vol_rb = M_rb;
regStruct.M_atl_to_vol_af = M_af;

% ToDos
%   Initialize with actual rb reg


