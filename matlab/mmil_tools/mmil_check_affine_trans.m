function [results,errcode] = mmil_check_affine_trans(vol_t2,vol_t1,vol_mask,varargin)
%function [results,errcode] = mmil_check_affine_trans(vol_t2,vol_t1,vol_mask,varargin)
%
% Purpose: estimate rigid and affine transforms between T1w and T2w scans
%   Rigid transform (rotation + translation) reflects registration accuracy
%   Warp measure reflects relative stretching (e.g., B0 warping)
%   Estimates warp measure seperate from Rigid transform for B0 Warp assesment
%    
% Usage: results = mmil_check_affine_trans(vol_t2,vol_t1,vol_mask,varargin)
%
% Required Input:
%   vol_t2 : T2(*)-weighted volume; e.g., 3D-T2, dMRI, fMRI
%   vol_t1 : T1-weighted volume volume
%   vol_mask : binary mask region of interest in which to check affine transform
%
% Optional Input:
%   'fname_plot': full path of output image file (jpg or tif)
%     {default = []}
%   'visible_flag': [0|1] show figure on screen
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Output:
%   results: Structure that holds the following output parameters
%     VolumeCorr   : corr( vol_t2(mask==1), vol_t1(mask==1) );
%     Delta_1      : |img_T2 -  img_T1|
%     Delta_2      : |Rectified(img_T2) -  img_T1|
%     MI_1         : MutualInfo(img_T2,imgT1)
%     MI_2         : MutualInfo(Rectified(img_T2),img_T1)
%     ImgCorr1     : Corr(img_T2,imgT1) within bounding box,
%                    where Rectified(img_T2) has data
%     ImgCorr2     : Corr(Rectified(img_T2),img_T1) within bounding box,
%                    where Rectified(img_T2) has data
%     ImgDist3     : 1-Corr(img_T2,imgT1) within mask, smoothed,
%                    where Rectified(img_T2) has data
%     ImgDist4     : 1-Corr(Rectified(img_T2),img_T1) within mask, smoothed,
%                    where Rectified(img_T2) has data
%     qAffineT     : argmin(|Rectify(img_T2,qAffineT*)-img_T1|)
%     RigidT       : argmin(|Rectify(img_T2, RigidT*)-img_T1|)
%     RawRegT1Norm : mean(| Rigid_T*unit_Circle_xy-unit_Circle_xy |)
%     TotRegT1Norm : RawRegT1Norm*(Delta_1/Delta_2)
%     WarpNorm     : Norm of Warp |[dx-1,dy-1]|, [dx,dy]=getScaleFactor(qAffineT);
%                    times (MI_2/MI_1)*10^tanh((MI_2-MI_1)/MI_1/0.05),
%                    times 1/(1+exp(-(ImgDist3-0.13)/0.02))
%     Figure=imgw  : Image (2D bitmap) of Matlab figure
%   errcode: returned as 1 if there was an error computing metric, 0 otherwise
%
% Notes:
%    | v | = sqrt(sum( vi^2)  Eculidian L2 Norm of vector v
%    RigidT   = (3 x 3) Transformation Matrix, Translation and Rot homogenous coords
%    qAffineT = (3 x 3) Transformation Matrix, Translation, Scale and Rot homogenous coords
%
% Created:  04/28/20 by Payam Saisan
% Prev Mod: 05/28/20 by Payam Saisan
% Prev Mod: 05/30/20 by Don Hagler
% Prev Mod: 08/03/23 by Don Hagler
% Prev Mod: 06/15/2026 by Don Hagler
% Prev Mod: 10/01/2026 by Don Hagler
% Prev Mod: 10/02/2026 by Don Hagler
% Prev Mod: 10/05/2026 by Don Hagler
% Prev Mod: 11/02/2026 by Don Hagler
% Prev Mod: 11/03/2026 by Don Hagler
% Last Mod: 11/07/2026 by Don Hagler
%

% initialize results
results = [];
results.VolumeCorr=NaN;
results.Delta_1=NaN;
results.Delta_2=NaN;
results.MI_1=NaN;
results.MI_2=NaN;
results.ImgCorr1=NaN;
results.ImgCorr2=NaN;
results.ImgDist3=NaN;
results.ImgDist4=NaN;
results.qAffineT=NaN;
results.RigidT=NaN;
results.RawRegT1Norm=NaN;
results.TotRegT1Norm=NaN;
results.WarpNorm=NaN;
results.Figure=NaN;
errcode = 0;

if ~mmil_check_nargs(nargin,3), return; end;

parms = mmil_args2parms(varargin, { ...
  'fname_plot',[],[],...
  'visible_flag',false,[false true],...
  'forceflag',false,[false true],...
  ...
  'title','Check Registration to T1',[],...
  'fig_size',[7,5],[],...
  'img_dpi',300,[],...
});

try
  % get bounding box
  vol_mask = 1.0*(vol_mask>eps);
  vmask = squeeze(vol_mask(round(end/2),:,:));
  rp = regionprops(vmask,'BoundingBox');
  xy_idx = [round(rp.BoundingBox(1:2)) floor(rp.BoundingBox(1:2)+rp.BoundingBox(3:4))];  
  xy_idx(4) = round(0.75*(xy_idx(4)-xy_idx(2))+xy_idx(2));
  %  
  vdti = vol_t2;
  vdti(~(vol_mask==1)) = 0;
  uti = squeeze(mean(vol_t2(round(end/2)+(-4:4),xy_idx(2):xy_idx(4),xy_idx(1):xy_idx(3)),1));
  %
  vmask_box = squeeze(vmask(xy_idx(2):xy_idx(4),xy_idx(1):xy_idx(3)));
  vmask_box = conv2(vmask_box,fspecial('gaussian',35,15),'same');
  vmask_box = (vmask_box/sum(vmask_box(:)));
  %
  vt1 = vol_t1;
  vt1(~(vol_mask==1)) = 0;
  ut1 = squeeze(mean(vol_t1(round(end/2)+(-4:4),xy_idx(2):xy_idx(4),xy_idx(1):xy_idx(3)),1));
  % normalize image intensities
  urgb(:,:,1) = normalize(ut1);
  urgb(:,:,2) = 1-normalize(uti);
  urgb(:,:,3) = ut1*0;
  %
  img1 = double( (255-histeq(uint8(round(normalize(uti)*255))) ));
  img2 = double( (histeq(uint8(round(normalize(ut1)*255))) ));
  
  img22 = normalize(conv2(img2,fspecial('gaussian',3,2),'same').*vmask_box);
  img11 = normalize(conv2(img1,fspecial('gaussian',3,2),'same').*vmask_box);
  
  % initialize optimizer
  [optimizer, metric] = imregconfig('monomodal');
  optimizer.MaximumIterations=400;
  
  % for over-all rigid transform norm
  tForm = imregtform(img11,img22,'rigid',optimizer,metric);
  % for affine warp norm
  tForm1 = imregtform(img11,img22,'similarity',optimizer,metric);
  % image reconstruction is based on rigid transform model
  % fill pixels shifted in from outside with NaN, so the mismatch is
  %   measured where the warped image has data
  img1_hat = imwarp(img1,tForm,'OutputView',imref2d(size(img2)),'FillValues',NaN);
  valid = ~isnan(img1_hat);
  img1_hat(~valid) = 0;
  % adjusted warp measure
  img1_hat2=imwarp(img1,tForm1,'OutputView',imref2d(size(img2)),'FillValues',NaN);
  valid2 = ~isnan(img1_hat2);
  % normalize image intensities
  clear urgb;
  urgb(:,:,1) = normalize(img2);
  urgb(:,:,2) = normalize(img1_hat);
  urgb(:,:,3) = ut1*0;
  % calculate metrics
  dist = corr( vol_t2(vol_mask==1), vol_t1(vol_mask==1) );
  % before and after the rigid fit, over the same pixels
  dist1 = calc_delta(img1(valid),img2(valid));
  dist2 = calc_delta(img1_hat(valid),img2(valid));
  % mutual information before and after the similarity fit, over the same pixels
  dist3 = mmil_information(img1(valid2)',img2(valid2)');
  % adjusted warp measure
  dist4 = mmil_information(img1_hat2(valid2)',img2(valid2)');
  dist5 = corr(img1(valid),img2(valid));
  dist6 = corr(img1_hat(valid),img2(valid));
  % normalize images
  urgb(:,:,2) = normalize(img1_hat);
  u_img_hat = urgb;
  urgb(:,:,2) = normalize(img1);
  u_img=urgb;
  img1(:,:,1) = normalize(img1);img1(:,:,2:3)=img1(:,:,[1 1]);
  img2(:,:,1) = normalize(img2);img2(:,:,2:3)=img2(:,:,[1 1]);
  % extract quasi affine warp (sx sy) scaling
  qt = TSR_decomp(tForm1.T') ;
  % calculate rigid transform norm
  tN = T_Norm(tForm);
  nN = tN.*((dist1+0.0001)./(dist2+0.0001));
  nN(find(isnan(nN)))=0;
  % adjusted warp measure
  wN=min((100*norm( [qt.Vals(4)-1, qt.Vals(3)-1])),100);
  %   weighted by 10^tanh of the relative MI change, in place of
  %   10^sign, so near-zero changes do not jump by a factor of 100
  wN2=norm( [qt.Vals(4)-1, qt.Vals(3)-1])* (((dist4)./(dist3)).*(10.^(tanh(((dist4-dist3)./dist3)/0.05))) );
  wN2=min(100*wN2/5,100);
  % correlation measures based on ROI mask rather than complete bounding box
  %   rigid-warped image smoothed and weighted as img11, over mask pixels
  %   the warped image covers, away from the smoothed fill border
  img11_hat = conv2(img1_hat,fspecial('gaussian',3,2),'same').*vmask_box;
  vmask_box = normalize(vmask_box)>0.5;
  pix = vmask_box & imerode(valid,true(3));
  dist7 = abs(1-corr(img11(pix),img22(pix))); % original dist
  dist8 = abs(1-corr(img11_hat(pix),img22(pix))); % dist for rectified
  % update wN2 (Warp) measure
  %   weighted by a logistic ramp in the mask-pixel dist7, centred on 0.13,
  %   in place of zeroing below 0.13, so scans near 0.13 do not jump to 0
  if ~isnan(dist7), wN2=wN2/(1+exp(-(dist7-0.13)/0.02)); end

  % plot summary figure
  imgw.cdata = [ [img1 u_img(:,1:5,:,:)*0+max(u_img(:)) img2 ]; [u_img u_img(:,1:5,:,:)*0+max(u_img(:)) u_img_hat]];
  if parms.visible_flag ||...
    ~isempty(parms.fname_plot) && (~exist(parms.fname_plot) || parms.forceflag)
    % create figure
    figure; clf;
    if ~parms.visible_flag, set(gcf,'Visible','off'); end;
    subplot(2,1,1);
    imagesc([ img1 u_img(:,1:5,:,:)*0+max(u_img(:)) img2]);
    title(parms.title); axis off;
    subplot(2,1,2);
    imagesc([u_img u_img(:,1:5,:,:)*0+max(u_img(:)) u_img_hat]);
    axis off;
    title([['|T_{D,Rot}| = ' , num2str(nN) '    |Warp2_{xy}| = ' num2str(wN2) '    |Warp_{xy}| = ' num2str(wN) '     D_{xy}=' num2str(norm(tForm.T(3,1:2)))  '  Rot= '  num2str(atan2(tForm.T(1,2),tForm.T(1,1))*180/pi)],10,...
        ['Delta_{1}=' num2str(dist1)  '  MI_{1}=' num2str(dist3)  '   --   Delta_{2}=' num2str(dist2)  '  MI_{2}=' num2str(dist4)]]);    
    set(gcf,'Position', [[560 367 707 550]]);
    % save figure
    imgw = getframe(gcf);
    drawnow;
    shg;
    if ~isempty(parms.fname_plot) && (~exist(parms.fname_plot) || parms.forceflag)
      mmil_save_fig(gcf,parms.fname_plot,...
        parms.fig_size,parms.img_dpi,0,parms.visible_flag);
    end;
  end
  % assign outputs
  results.VolumeCorr = dist;
  results.Delta_1 = dist1;
  results.Delta_2 = dist2;
  results.MI_1 = dist3;
  results.MI_2 = dist4;
  results.ImgCorr1 = dist5;
  results.ImgCorr2 = dist6;
  results.ImgDist3 = dist7;
  results.ImgDist4 = dist8;
  results.qAffineT = tForm1.T;
  results.RigidT = tForm.T;
  results.RawRegT1Norm = tN;
  results.TotRegT1Norm = nN;
  results.WarpNorm = wN2;
  results.Figure = imgw.cdata;
catch me
  fprintf('%s: WARNING: error estimating affine transform:\n%s\n',...
    mfilename,me.message);
  errcode = 1;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% norm of transform tForm
function tN=T_Norm(tForm)
  pts=[cos( (0:0.01:1)*2*pi)' sin((0:0.01:1)*2*pi)' sin((0:0.01:1)*2*pi)'*0+1]';
  for k=1:length(pts), px(:,k)= tForm.T'*pts(:,k); tNorm(k)=norm(pts(:,k)-px(:,k)/px(3,k));end
  tN=mean(tNorm);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% factor out trans, scale, shear, rot from affine transform
function qt=TSR_decomp(affineTrans)
  % get transform parameters
  a=affineTrans(1,1);
  b=affineTrans(1,2);
  c=affineTrans(1,3);
  d=affineTrans(2,1);
  e=affineTrans(2,2);
  f=affineTrans(2,3);
  A=[a b;d e];
  %
  p=sqrt(a^2+b^2);
  %
  r=det(A)/p;
  %
  q=(a*d+b*e)/det(A);
  th=atan2(b,a);
  %
  t = [1 0 c;0 1 f; 0 0 1];
  sh = [1 0 0;q 1 0; 0 0 1];
  sc = [p 0 0;0 r 0;0 0 1];
  rot_M = [cos(th) sin(th) 0;-sin(th) cos(th) 0; 0 0 1];
  %
  dx=c; dy=f; q_shear=q; sx=p; sy=r; theta=th;
  qt.T=[t  sc sh rot_M];
  qt.Vals=[dx dy sx sy q_shear theta*180/pi];
  qt.key='qt.Vals=[dx dy sx sy q_shear theta] 10 qt.T=[TransM, Scale_M, Shear_M, Rotation_M]';
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = calc_delta(a,b)
  a = normalize(a); b = normalize(b);
  d = sum((a-b).^2)/sum((a+b).^2);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imgg=normalize(imgg)
  for i=1:size(imgg,3)
    img=imgg(:,:,i);
    y=img-min(img(:));
    if ~(max(y(:))==0)
      y=y/max(y(:));
    end
    imgg(:,:,i)=y;
  end
return;

