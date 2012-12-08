% function [x,y]=ellipse(x0,y0,rx,ry,theta0,n);
function [x,y]=ellipse(x0,y0,rx,ry,theta0,n);

if ~exist('theta','var'),
   theta=0;
end
if ~exist('n','var'),
   n=21;
end

x=zeros(n,1);
y=x;

rr=(1:n)';
theta=2*pi/(n-1).*(rr-1)+theta0;
x(rr)=rx*cos(theta)+x0;
y(rr)=ry*sin(theta)+y0;

