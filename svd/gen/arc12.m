% function [x,y]=arc12(x1,y1,x2,y2,r,n);
%
% return plot points for an arc running ccw from (x1,y1) to (x2,y2)
% with radius r
%
% created SVD 3/7/05 (wrapper for arc.m)
%
function [x,y]=arc12(x1,y1,x2,y2,r,n);

if ~exist('n','var'),
   n=21;
end

if sign(r)<1,
   xt=x1; yt=y1;
   x1=x2; y1=y2;
   x2=xt; y2=yt;
end

xm=(x1+x2)./2;
ym=(y1+y2)./2;
d=sqrt((x1-xm).^2+(y1-ym).^2);
b=sqrt(r.^2-d.^2);

btan=[x2-x1; y2-y1];
btan=btan./norm(btan);
bperp=[-(y2-y1); x2-x1];
bperp=bperp./norm(bperp);

p0=[xm;ym]+bperp*b;
x0=p0(1);
y0=p0(2);

th1=angle(x1-x0 + (y1-y0)*i);
th2=angle(x2-x0 + (y2-y0)*i);

[x,y]=arc(x0,y0,th1.*180/pi,th2.*180/pi,abs(r),n);

if sign(r)<1,
   x=flipud(x);
   y=flipud(y);
end

