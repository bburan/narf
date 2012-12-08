% function [psi,u,v]=gabcurve(ang,sf,ph,sigma0,sigma1,rad,x,y,[N or u,v,s])
%
% by default the window is one unit of space unless (u,v) are
% specified differently. generally for the case where u and v range
% from -0.5 to 0.5, parameters in units of space should be < 0.5
%
% basic gabor parameters:
% ang - angle of center (in radians, from x axis)
% sf - sf (cycles per unit of space)
% ph - phase of gabor (radians)
% sigma0 - radius of envelope along gabor (units of space)
% sigma1 - radius of envelope across gabor (units of space)
% rad - radius of curvature (units of space)
% x,y - position of gabor center (units of space)  (default 0,0)
%
% sampling parameters
% N - sampling (  u=v=linspace(-0.5,0.5,N)  ) (default 16)
% u,v - vectors specifying sampling
% s - scale (for changing wavelet scale)
% 
% created SVD 6/23/04
%
function [psi,u,v]=gabcurve(ang,sf,ph,sigma0,sigma1,rad,x,y,u,v,s)

if ~exist('x','var'),
   x=0;
end
if ~exist('y','var'),
   y=0;
end
if ~exist('u','var'),
   u=16;
end
if ~exist('v','var'),
   if length(u)>1,
      N=length(u);
   else
      N=u;
      u=linspace(-0.5,0.5,N);
   end
   v=u;
   s=1;
elseif ~exist('s','var'),
   s=1;
end

% recenter in window
xs=(x-u)./s;
ys=(y-v)./s;

% recenter at appropriate curvature radius
xs=xs-rad*cos(ang);
ys=ys+rad*sin(ang);

[XX,YY]=meshgrid(xs,ys);

TH=atan2(YY,XX)+pi;
RR=sqrt(XX.^2+YY.^2);

TH=mod(TH+ang+pi,2*pi)-pi;

g=exp(-(RR-rad).^2./(2*sigma0^2) ...
      -(TH.*rad).^2./(2*sigma1^2))./ ...
       (2*pi*sqrt(sigma0*sigma1));
f=exp(-i*(2.*pi*sf.*(RR-rad)+ph));

psi=real((f.*g)./s);


return

% parms for circular gaussian windows
%sigma=sigma0;
%g=exp(-(XX.^2+YY.^2)./(2*sigma))./(2*pi*sigma);
%f=exp(-i*(sf.*(XX.*cos(ang)+YY*sin(ang))+ph));

g=exp(-(XX.*cos(ang)+YY*sin(ang)).^2./(2*sigma0^2) ...
      -(XX.*cos(ang+pi/2)+YY*sin(ang+pi/2)).^2./(2*sigma1^2))./ ...
       (2*pi*sqrt(sigma0*sigma1));
f=exp(-i*(2.*pi*sf.*(XX.*cos(ang)+YY*sin(ang))+ph));

psi=(f.*g)./s;
