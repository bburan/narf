% function pal=redblue(expon);
%
% passed to colormap() function to set blue-->white-->red palette
%
% expon is the exponent of coloration. 
% (ie, 1= linear, <1 more white, >1 more color)
% 
function pal=gray2(expon);

if ~exist('expon','var'),
   expon=1.0;
end

pal=[linspace(0,1,64).^expon' linspace(0,1,64).^expon' ...
     linspace(0,1,64).^expon'];

