% function psi=gabor(ang,sf,ph,sigma0,sigma1,x,y,u,v,s)
%
function psi=gabor(ang,sf,ph,sigma0,sigma1,x,y,u,v,s)

xs=(x-u)./s;
ys=(y-v)./s;

[XX,YY]=meshgrid(xs,ys);

% parms for circular gaussian windows
%sigma=sigma0;
%g=exp(-(XX.^2+YY.^2)./(2*sigma))./(2*pi*sigma);
%f=exp(-i*(sf.*(XX.*cos(ang)+YY*sin(ang))+ph));

g=exp(-(XX.*cos(ang)+YY*sin(ang)).^2./(2*sigma0^2) ...
      -(XX.*cos(ang+pi/2)+YY*sin(ang+pi/2)).^2./(2*sigma1^2))./ ...
       (2*pi*sqrt(sigma0*sigma1));
f=exp(-i*(2.*pi*sf.*(XX.*cos(ang)+YY*sin(ang))+ph));

psi=(f.*g)./s;
