% function [x,y]=arc(x0,y0,th1,th2,r,n[=21]);
%
% return plot points for an arc centered at (x0,y0) running from
% angle th1 to th2 (in degrees) at radius r, sampled at n points.
%
% angles in complex plane coordinates (ie, 0= x axis right, 90 = y axis up)
%
% created SVD 2/14/05
%
function [x,y]=arc(x0,y0,th1,th2,r,n);

if not(exist('n')),
   n=21;
end

x=zeros(n,1);
y=x;

for rr=1:n,
   theta=mod(th2-th1,360)/(n-1)*(rr-1)+th1;
   theta=theta/180*pi;
   x(rr)=r*cos(theta)+x0;
   y(rr)=r*sin(theta)+y0;
end

