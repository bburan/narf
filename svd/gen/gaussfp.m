% function f=gaussfp(beta,x)
%
% f=a/(sqrt(2*pi*s^2)) .* exp(-(x-m).^2/(2*s^2)) + d
%
% where beta=[or sf orw sfw a d]
%
% ie, 2D gaussian with mean (or,sf), peak height a (above d),
% offset d, evaluated over the values in x
%
% CURRENTLY FORCING d=0
%
% NOTE: hard-coded switch for either polar or cartesian orientation
%       bandwidth (orw). which one is better??
%
function f=gaussfp(beta,x)

%f=exp(-(x-m).^2./(2.*s^2))./(sqrt(2*pi)*s);
%f=f./max(f).*a + d;

or=beta(1);
sf=beta(2);
orw=beta(3);
sfw=beta(4);
a=beta(5);
d=0;

if ~exist('x','var'),
   [xx,yy]=meshgrid(-8:7,-8:7);
   x=[xx(:) yy(:)];
end

rmat=[cos(or) sin(or); -sin(or) cos(or)];
x1=x*rmat;
x1(find(x1(:,1)==0),1)=0.0001;

if 1,
   % polar orientation bw
   r=sqrt(x1(:,1).^2+x1(:,2).^2);
   r(find(r==0))=0.0001;
   
   ang=atan(x1(:,2)./x1(:,1));

   f=exp(-(log(r)-log(sf)).^2./(2.*log(sfw)^2)- ...
         (ang).^2./(2.*orw^2)) ./ (2*sfw*orw*pi);
else
   % cartesian orientation bw
   f=exp(-(log(abs(x1(:,1)))-log(sf)).^2./(2.*log(sfw)^2)- ...
         (x1(:,2)).^2./(2.*orw^2)) ./ (2*sfw*orw*pi);
end

f=f.*a + d;
