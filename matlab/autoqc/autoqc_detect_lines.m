function [line_score,line_count] = autoqc_detect_lines(vol,nslices,outdir,outstem,verbose,plotflag)
%function [line_score,line_count] = autoqc_detect_lines(vol,nslices,outdir,outstem,verbose,plotflag)
%
% Purpose: Determine if scan has line artifacts
%
% Required Input:
%   vol: 3D scan volume
%
% Optional Input:
%   nslices: total number of slices to be analyzed
%     if 0, will analyze all slices
%     {default = 1}
%   outdir: output directory
%     {default = pwd}
%   outstem: output file stem
%     {default = 'vol'}
%   verbose: [0|1] display messages about lines detected
%     {default = 0}
%   plotflag: [0|1] plot image if lines detected
%     {default = 1}
%
% Output:
%   line_score: scalar measure of line artifact (averaged across slices tested)
%   line_count: number of lines detected (averaged across slices tested)
%
% Created:   08/29/21 by Payam Saisan
% Prev Mod:  11/29/21 by Payam Saisan
% Prev Mod:  02/01/22 by Don Hagler
% Last Mod:  06/26/2026 by Don Hagler
% 

line_score = 0; line_count = 0;

if ~mmil_check_nargs(nargin,1), return; end
if ~exist('nslices','var') || isempty(nslices), nslices = 1; end
if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end
if ~exist('outstem','var') || isempty(outstem), outstem = 'vol'; end
if ~exist('verbose','var') || isempty(verbose), verbose = 0; end
if ~exist('plotflag','var') || isempty(plotflag), plotflag = 1; end

if length(size(vol))>2
  dslices = size(vol,3);
else
  % allow single slice as input
  tmpvol(1,:,:) = squeeze(vol);
  vol = tmpvol;
  nslices = 1;
  dslices = 1;
end
if nslices<1 || nslices>dslices
  nslices = dslices;
end

% if multi-frame, use middle frame
if length(size(vol))==4
  vol = vol(:,:,:,round(end/2));
end

mmil_mkdir(outdir);

line_scores = zeros(nslices,1);
line_counts = zeros(nslices,1);
for i=1:nslices
  % set to actual slice number for cross-plane
  slicenum = round((size(vol,1)/(nslices+1))*i);
  imgo=normalize(rot90(squeeze(vol(slicenum,:,:))));
  % filter and smooth
  img=medfilt2(imgo,[1 14]);
  img=conv2(img,fspecial('gauss',3,1),'same');
  % create mask / erode
  msk=double(img>max(img(:))*0.15);
  se = strel('disk',3);
  msk = imerode(msk,se);
  msk=conv2(double(msk),fspecial( 'gauss',4,2),'same');
  % find edges
  [xe,ye] = find(edge(msk,'canny'));
  msk(msk==0) = NaN;
  img = double(normalize(imgo)).*msk;
  if length(size(img))>2
    rot = double(rgb2gray(img));
  else
    rot = img;
  end
  rotI = double(rot);
  dims = max(size(rotI));
  mn = abs(diff(nanmean(normalize(img'))));
  rotI(rotI<0.1*max(rotI(:))) = min(rotI(:));
  rotI = rotI.*msk;
  rf = sum(msk)/prod(size(msk));
  rotI = sigmoid(rotI,nanmean(rotI(:)),10);
  mn = diff(nanmean(rotI'));
  mn = [mn 0].*nanmean(rotI'>0.1*max(rotI(:)));
  mn = diff(mn);
  [pks,idx] = findpeaks(mn);
  thresh = 0.2;
  idx2 = find(pks>thresh);
  line_scores(i) = sum(mn(idx(idx2)));
  line_counts(i) = length(idx2);
  
  % plot image if line detected
  if plotflag && ~isempty(idx2)
    figure; clf;
    colormap(bone);
    imagesc(normalize(img)); hold on;
    plot(ye,xe,'r.');
    title(sprintf('%s: slice=%d score=%0.2f lines=%d',...
      regexprep(outstem,'_',' '),slicenum,line_scores(i),line_counts(i)),...
      'FontSize',8);
    for p=1:line_counts(i)
      xd=1:1:size(img,2);
      yd=0*(1:1:size(img,2))+idx(idx2(p))+1;
      plot( xd([1:10  (end-10):end]),yd([1:10  (end-10):end]) ,'y.');
      plot( xd,yd ,'b.','MarkerSize',0.01);
    end
    imgw=img;
    outfile = sprintf('%s/%s_lines_slice%03d.jpg',outdir,outstem,slicenum);
    mmil_save_fig(gcf,outfile);
  end
end

line_score = mean(line_scores);
line_count = mean(line_counts);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imgg = normalize(imgg)
  for i=1:size(imgg,3)
    img=imgg(:,:,i);
    y=img-min(img(:));
    if ~(max(y(:))==0)
      y=y/max(y(:));
    end
    imgg(:,:,i)=y;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = sigmoid(z,c1,c2)
  if nargin==1
    c1=1;
    c2=2;
  end
  a = 1.0 ./ (1+  exp(-c2*(z-c1)));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


