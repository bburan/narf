% function pal=blackwhite(expon);
%
% passed to colormap() function to set blue-->white-->red palette
%
% expon is the exponent of coloration. 
% (ie, 1= linear, <1 more white, >1 more color)
% 
function pal=blackwhite(expon);

if ~exist('expon','var'),
   expon=0.75;
end

skbottom=[linspace(0,1,32).^expon' linspace(0,1,32).^expon' ...
          linspace(0,1,32).^expon']./2;
pal=[skbottom; 1-flipud(skbottom)];

