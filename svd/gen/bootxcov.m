% function [cxy,exy,tt,p]=bootxcov(x,y,maxlag,N);
%
% cxy is the cross covariance of x and y.  then x is shuffled
% multiple times to get a measure of expected error cross
% covariance, exy, for 95% of shuffled trials.
%
% maxlag-number of bins +/- from 0 in cross-cov
% N-number of shuffles (default 100)
%
% created: SVD 10/17/04 (ripped off randxcov)
%
function [cxy,exy,tt,p]=bootxcov(x,y,maxlag,N);

if ~exist('N','var'),
   N=100;
end
if ~exist('maxlag','var'),
   maxlag=0;
end

gidx=find(~isnan(x+y));
x=x(gidx);
y=y(gidx);

L=length(x);
epsilon=1e-100;

p=zeros(maxlag*2+1,1);

mx=mean(x);
my=mean(y);

if sum(abs(x-mx))<=epsilon | sum(abs(y-my))<=epsilon,
   cxy=zeros(maxlag*2+1,1);
   exy=zeros(maxlag*2+1,1);
   tt=(-maxlag:maxlag)';
   p=-ones(maxlag*2+1,1);
   return
end

[cxy,tt]=xcov(x,y,maxlag,'coeff');

% bootstrap to compute error bars
bexy=zeros(maxlag*2+1,N);
bstep=length(x)./N;
for n=1:N,
   bidx=[1:round((n-1)*bstep) round(n*bstep+1):length(x)];
   bexy(:,n)=xcov(x(bidx),y(bidx),maxlag,'coeff');
end
exy=std(bexy).*sqrt(N+1)./sqrt(N);

return





