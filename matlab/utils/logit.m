function y = logit(x,invflag)
%function y = logit(x,invflag)
%
% logit transformation
%
% Input:
% ------
% x,            input value
% invflag,      logical, inverse?
%
% Return:
% ------
% y,            transformed value
%
% Created:  10/24/2011 by Anders Dale
% Last Mod: 02/25/2014 by Anders Dale
% Last Mod: 09/04/2020 by Don Hagler
%

if exist('invflag','var') & invflag
  y = exp(x)./(1+exp(x));
else
  y = log(x./(1-x));
end

