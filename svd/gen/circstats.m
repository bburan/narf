% function [cm,cstd,cv]=circstats(r,theta,display);
%
% r - length of vectors
% theta - corresponding angle for each vector(radians)
%         default assume evenly spaced values from 0 to 2pi*(1-1/N)
%         where N is length(r)
%
% outputs:
% cm - circular mean (in radians)
% cstd - circular standard dev (also radians)
% cv - circular variance. 0=no variance, 1=circular symmetry
%
% created SVD 10/04
%
function [cm,cstd,cv]=circstats(r,theta,display);

if ~exist('theta','var'),
   theta=linspace(0,2.*pi,length(r)+1);
   theta=theta(1:end-1);
end

x=r(:).*exp(i * theta(:));
mx=mean(x);
cm=angle(mx);

if sum(abs(r))>0,
   cv=1-abs(sum(x))./sum(abs(r));
else
   cv=0;
end

cstd=sqrt(-2 .* log(1-cv));

if exist('display','var') & display,
   plot(real(x),imag(x),'kx');
   hold on
   plot(real(mx),imag(mx),'ro');
   mlen=abs(mx);
   mpe=[mx.*exp(i*cstd) mx mx.*exp(-i*cstd)];
   plot(real(mpe),imag(mpe),'-');
   hold off
   axis([-1 1 -1 1].*max(abs(r)))
end


return

% demo code:

% aka, demonstrate that circular std is the same as the standard
% deviation (in radians) of a linear 1d gaussian function wrapped
% around the circle ... for std<pi/4 or so.

theta=linspace(-pi,pi-pi/200,400);

Nsig=100;
cvar=zeros(Nsig,1);
cst=zeros(Nsig,1);
sigrange=linspace(pi/30,pi/4,Nsig);
for sigma=1:Nsig,
   ss=sigrange(sigma);
   r=gauss1([0 ss 1 0],theta);
   [cm,cv,cstd]=circstats(r,theta);
   cvar(sigma)=cv;
   cst(sigma)=cstd;
end

plot(cst,sigrange.*180/pi);
xlabel('circular std dev');
ylabel('standard dev of gaussian (deg)');

p = polyfit(cst,sigrange'.*180/pi,1)

title(sprintf('std dev (deg)= %.1f * circular std dev',p(1)));

%hold on
%plot(cst,sigrange.*180/pi,'r');
%%hold off
