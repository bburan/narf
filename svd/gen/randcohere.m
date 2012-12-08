% function [cxy,exy,tt]=randcohere(x,y,nfft);
%
% cxy is the cross covariance of x and y.  then x is shuffled
% multiple times to get a measure of expected error cross
% covariance, exy
%
function [cxy,exy]=randcohere(x,y,nfft);

N=100;
L=length(x);
epsilon=1e-100;

if sum(abs(x))<=epsilon | sum(abs(y))<=epsilon,
   cxy=zeros(nfft/2+1,1);
   exy=zeros(nfft/2+1,1);
   return
end

cxy=cohere(x,y,nfft,1);

for n=1:N,
   s=round(rand(3,1)*L)+1;
   %tx=(x([s(1):L 1:s(1)-1])+x([s(2):L 1:s(2)-1])+x([s(3):L 1:s(3)-1]))/3;
   tx=x([s(1):L 1:s(1)-1]);
   if n==1,
      exy=cohere(tx,y,nfft,1);
   else
      exy=exy+cohere(tx,y,nfft,1);
   end
end
exy=exy/N;

% zero out coherence at frequency zero (DC), since everything is
% rectified, the value will be high and meaningless
if 1,
   exy(1)=0;
   cxy(1)=0;
end

return




