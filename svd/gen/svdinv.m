% function result=svdinv(A,maxsingratio)
%
% A:		near-singular matrix to invert
% maxsingratio:	use singular values >= max(singular value)*maxsingratio
%			(default=0.01 -- too low????)
%

function result=svdinv(A,maxsingratio)

[U S V] = svd(A);
len=size(S);

%if not specified, choose singular value cutoff at 1/100 of maximum
if nargin<2,
   maxsingratio=1/100
end
maxsing=S(1,1)*maxsingratio;

%Sinv=zeros(len(1),len(2));
for ii=1:len(1),
   if S(ii,ii) < maxsing | S(ii,ii)==0,
%      Sinv(ii,ii)=0;
      S(ii,ii)=0;
   else
%      Sinv(ii,ii)=1/S(ii,ii);
      S(ii,ii)=1/S(ii,ii);
   end
end
result=V*S*U';

