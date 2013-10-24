% function [cxy,exy,p]=rfloor(x,y,N);
%
% cxy - r value between x and y, ie, xcov(x,y,0,'coeff')
%
% then x is shuffled multiple times to get a measure of expected error cross
% covariance, exy, for 95% of shuffled trials.
%
% N-number of shuffles (default 100)
%
% created: SVD 2013-10-24, ripped off randxcov
%
function [cxy,exy,p]=rfloor(x,y,N);

if ~exist('maxlag','var'),
   maxlag=0;
end
if ~exist('N','var'),
   N=100;
end
if ~exist('Nerr','var'),
   Nerr=0.05;
end
L=length(x);
epsilon=1e-100;

p=zeros(maxlag*2+1,1);

mx=mean(x);
my=mean(y);

if sum(abs(x-mx))<=epsilon | sum(abs(y-my))<=epsilon,
   cxy=zeros(maxlag*2+1,1);
   exy=zeros(maxlag*2+1,1);
   p=-ones(maxlag*2+1,1);
   return
end

[cxy,tt]=xcov(x,y,maxlag,'coeff');

bexy=zeros(maxlag*2+1,N);
for n=1:N,
   %s=round(rand(3,1)*L)+1;
   %tx=(x([s(1):L 1:s(1)-1])+x([s(2):L 1:s(2)-1])+x([s(3):L 1:s(3)-1]))/3;
   %tx=(x([s(1):L 1:s(1)-1]));
   tx=shuffle(x);
   bexy(:,n)=(xcov(tx,y,maxlag,'coeff'));
end

%exy=sqrt(sum(bexy,2)./N);
bexy=sort([bexy cxy],2);
exy=bexy(:,round((1-Nerr)*N));  %  exy=95% error for each latency
for ii=1:length(cxy),
    p(ii)=sum(cxy(ii)<=bexy(ii,:))./(N+1);
end

%keyboard
return

[cxy,exy,tt]=randxcov(rectpred(trespvalididx),trespvalid,100);
figure(11)
plot(tt,cxy)
hold on
plot(tt,exy,'--');
hold off



