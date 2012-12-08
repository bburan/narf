% function pal=redblue(expon);
%
% passed to colormap() function to set blue-->white-->red palette
%
% expon is the exponent of coloration. 
% (ie, 1= linear, <1 more white, >1 more color)
% 
function pal=redblue(expon);

if 1,
   if ~exist('expon','var'),
      expon=0.75;
   end
   skbottom=[(linspace(0,1,32).^expon)' (linspace(0,1,32).^expon)' ones(32,1)];
   sktop=[ones(32,1) linspace(1,0,32).^expon' linspace(1,0,32).^expon'];
   pal=[skbottom;sktop];
else
   if ~exist('expon','var'),
      expon=1.0
   end
   grc=0.5;
   rb=0.5;
   
   skbottom=[linspace(0,grc,32).^expon' linspace(0,grc,32).^expon' ...
             linspace(1,grc,32).^expon'];
   sktop=[ linspace(grc,1,32).^expon' linspace(grc,rb,32).^expon' ...
           linspace(grc,rb,32).^expon'];
   pal=[skbottom;sktop];
end

