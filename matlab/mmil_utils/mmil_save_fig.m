function mmil_save_fig(fig_num,fname_out,fig_size,img_dpi,eps_flag,visible_flag)
%function mmil_save_fig(fig_num,fname_out,fig_size,img_dpi,eps_flag,visible_flag)
%
% created:  04/20/16 by Don Hagler
% last mod: 02/09/20 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('fig_size','var') || isempty(fig_size), fig_size = [4 4]; end;
if ~exist('img_dpi','var') || isempty(img_dpi), img_dpi = 300; end;
if ~exist('eps_flag','var') || isempty(eps_flag), eps_flag = 0; end;
if ~exist('visible_flag','var') || isempty(visible_flag), visible_flag = 0; end;

if ~visible_flag, set(fig_num,'visible','off'); end;
set(fig_num, 'PaperUnits', 'inches');
set(fig_num, 'PaperSize', [fig_size(1) fig_size(2)]);
set(fig_num, 'PaperPositionMode', 'manual');
set(fig_num, 'PaperPosition', [0 0 fig_size(1) fig_size(2)]);
[fpath,fstem,fext] = fileparts(fname_out);
options = sprintf('-r %d',img_dpi);
switch fext
  case '.tif'
    print(fig_num,'-dtiff',fname_out,options);
  case '.jpg'
    print(fig_num,'-djpeg80',fname_out,options);
  otherwise
    eps_flag = 1;
end;
if eps_flag
  fname_eps = regexprep(fname_out,fext,'.eps');
  mmil_printeps(fig_num,fname_eps);
end;
if ~visible_flag, close(fig_num); end;

return;

