% function tf=gaborsf(beta,x,[fitsize]=16)
%
% beta=[or,orw,sf,sfw,A]
function tf=gaborsf(beta,x,fitsize)

if ~exist('fitsize'),
   fitsize=16;
end

or=beta(1);
orw=beta(2);
sf=beta(3);
sfw=beta(4);
A=beta(5);
offset=beta(6);

%[xx,yy]=meshgrid(1:fitsize+1,1:fitsize+1);
%x0=round((fitsize+1)/2);
%y0=x0;
%xx=(xx-x0)./fitsize*16;
%yy=(yy-y0)./fitsize*16;

[xx,yy]=meshgrid(-round(fitsize/sqrt(2)+1):round(fitsize/sqrt(2)+1),...
                 -round(fitsize/sqrt(2)+1):round(fitsize/sqrt(2)+1));
xx=xx./fitsize*16;
yy=yy./fitsize*16;

xc=round(size(xx,1)/2);
x1=xc-fitsize/2;
x2=xc+fitsize/2-1;

orwpix=sf*tan(orw*pi/180);

%if (orw > 60 | sfw > 4 | sf > 7.5)
%   tf=zeros(fitsize*fitsize,1);
%   return
%end

%yy=exp(abs(yy))./exp(sf);
%sfw=sfw./exp(sf);
%sf=1;
yy=yy+(yy==0).*0.01;
yy=log(abs(yy))-log(sf);
sfw=(log(sfw/2)+log(sf)); %sfw-log(sf);
sf=0;

% SVD 5/1/01 : got rid of square of exponent term
tf=A.*exp(-(abs(yy)-sf).^2/(2*sfw.^2)-xx.^2./(2*orwpix.^2)) + offset;

% rotate to appropriate orientation
tf=imrotate(tf,-or,'bilinear','crop');

% crop off edges from rotated image
%tf=reshape(tf(1:fitsize,1:fitsize),fitsize*fitsize,1);
tf=reshape(tf(x1:x2,x1:x2),fitsize*fitsize,1);

