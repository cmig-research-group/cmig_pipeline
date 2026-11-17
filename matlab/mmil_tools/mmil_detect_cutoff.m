function [cutoff,dorsal,ventral]=mmil_detect_cutoff(fname_in,varargin)
%function [cutoff,dorsal,ventral]=mmil_detect_cutoff(fname_in,[options])
%
% Purpose: determine the measure of frame cutoff in scans 
%          assumes masks passed to it are derived from "raw" data
%
% Usage: [cutoff,dorsal,ventral] = mmil_detect_cutoff(fname_in,'key1', value1,...);
%
% Required Input:
%   fname_in: full path of input mask volume file
%         
% Optional Input:
%   'fname_plot': full path of output image file (jpg or tif)
%
% Output:
%   cutoff: total cutoff score (sum of dorsal and ventral)
%   dorsal: cutoff score for dorsal edge
%   ventral: cutoff score for ventral edge
%
% Created:  02/07/20 by Payam Saisan
% Prev Mod: 02/13/20 by Don Hagler
% Prev Mod: 02/01/22 by Don Hagler
% Prev Mod: 02/01/23 by Don Hagler
% Last Mod: 07/16/24 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cutoff = 0; dorsal = 0; ventral = 0;
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
  'fname_plot',[],[],...
  'forceflag',false,[false true],...
  ...
  'fig_size',[5,5],[],...
  'img_dpi',300,[],...
  'visible_flag',false,[false true],...
});

[vol,M] = fs_load_mgh(fname_in);

% reorient volume to standard
[vol_PRS,M_PRS] = fs_reorient(vol,M,'PRS');

% determine mask and frame border intersetctions
[b_data img_data c_img]=mask_border_intersection((rot90(squeeze((mean(vol_PRS(:,:,:,:),4))))));

% create scores based on intersections
[dd1 ss_top1] = feat1d(b_data{1}.top) ;
[dd2 ss_bot1] = feat1d(b_data{1}.bot) ;
[dd3 ss_top2] = feat1d(b_data{2}.top) ;
[dd4 ss_bot2] = feat1d(b_data{2}.bot) ;

dorsal=mean([ss_top1,ss_top2]);
ventral=mean([ss_bot1,ss_bot2]);
cutoff=dorsal+ventral;

if ~isempty(parms.fname_plot) &&...
  (~exist(parms.fname_plot) || parms.forceflag)
  figure; clf;
  if ~parms.visible_flag, set(gcf,'Visible','off'); end;
  imagesc(c_img,[0,max(1,prctile(c_img(:),99))]);
  set(gca,'XTickLabel','','YTickLabel','');
  colormap hot;
  title(sprintf('dorsal cutoff = %0.1f',dorsal));
  xlabel(sprintf('ventral cutoff = %0.1f',ventral));
  mmil_save_fig(gcf,parms.fname_plot,...
    parms.fig_size,parms.img_dpi,0,parms.visible_flag);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dd score]=feat1d(d)
  dd=d>max(d(:))/1.5;
  score=round(sum(dd)/length(dd)*100);
return;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [bi im cimg]=borderI(imgg)
% analyze scan stack imgg for intersection with frame borders
% bi   = border intersetctions strcut
% im   = image struct data
% cimg = combined 2 views showing intersection of mask/scan with frame

function [bi im cimg]=mask_border_intersection(imgg)
  for k=1:3
    img= normalize(squeeze(sum(imgg,k)));    
    img=double(img>0.1) ;
    img=conv2(img,fspecial('gauss',3,1),'same'); 
    bi{k}.top=img(:,end); 
    bi{k}.bot=img(:,1); 
    im{k}.data=img;
    im{k}.edge=edge(img,'prewitt');
  end
  cimg = [(rot90(im{2}.edge+normalize(squeeze(mean(imgg,2))))) rot90(im{1}.edge+normalize(squeeze(mean(imgg,1))))];
  cimg(:,round(end/2))=max(cimg(:));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=normalize(img)
  y=img-min(img(:));
  y=y/max(y(:));
return;







