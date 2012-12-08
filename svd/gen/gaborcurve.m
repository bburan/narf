% function [fe,fo]=gaussfp(beta,x)
%
% f=a/(sqrt(2*pi*s^2)) .* exp(-(x-m).^2/(2*s^2)) + d
%
% where beta=[or sf orw sfw a d]
%
% ie, 2D gaussian with mean (or,sf), peak height a (above d),
% offset d, evaluated over the values in x
%
function [fe,fo,r]=gaborcurve(beta,x)

or=beta(1);
sf=beta(2);
orw=beta(3);
sfw=beta(4);
a=beta(5);
d=0;

pixcount=32;
if ~exist('x','var'),
   [xx,yy]=meshgrid(linspace(-16,15,pixcount),linspace(-16,15,pixcount));
   x=[xx(:) yy(:)];
end

rmat=[cos(or) sin(or); -sin(or) cos(or)];
x1=x*rmat;
x1(find(x1(:,1)==0),1)=0.0001;

xx=reshape(x1(:,1),pixcount,pixcount);
yy=reshape(x1(:,2),pixcount,pixcount);

if 1,
   % polar orientation bw
   r=sqrt(xx.^2+yy.^2);
   r(find(r==0))=0.0001;
   
   ang=atan(yy./xx);
   
   f=exp(-(log(r)-log(sf)).^2./(2.*log(sfw)^2)- ...
         (ang).^2./(2.*orw^2)) ./ (2*sfw*orw*pi);
else
   % cartesian orientation bw
   f=exp(-(log(abs(xx))-log(sf)).^2./(2.*log(sfw)^2)- ...
         (yy).^2./(2.*orw^2)) ./ (2*sfw*orw*pi);
end

f=f.*a + d;


fe=real(fftshift(ifft2(fftshift(f))));
f(xx<0)=-i.*f(xx<0);
f(xx>0)=i.*f(xx>0);
fo=real(fftshift(ifft2(fftshift(f))));

