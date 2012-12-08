% function result=svdinv(A,maxsing)
%
% A:		near-singular matrix to invert
% maxsingratio:	use singular values >= max(singular value)*maxsingratio
%			(default=0.01 -- too low????)
%

function result=svdinv(A,maxsing)

[U S V] = svd(A);
len=size(S);

%if not specified, choose singular value cutoff at 1/100 of maximum
%if nargin<2,
%   maxsing=1/100
%end
%maxsing=S(1,1)*maxsingratio;

if exist('maxsing'),
   Emin=max(find(diag(S) > maxsing));
else
   E=zeros(len(1),1);
   pixperside=sqrt(len(1));
   [XX,YY]=meshgrid(-ceil((pixperside-1)/2):floor((pixperside-1)/2),-ceil((pixperside-1)/2):floor((pixperside-1)/2));
   Ws=1./(sqrt(XX.^2 + YY.^2 +1));
   Ws=reshape(Ws,len(1),1);
   W=(Ws*Ws');

   Sinv=zeros(len(1),len(2));

   for cutoff=1:len(1),
      if S(cutoff,cutoff)==0,
         if cutoff>1,
            E(cutoff)=E(cutoff-1);
         end
      else
         Sinv(cutoff,cutoff)=1/S(cutoff,cutoff);
         result=V*Sinv*U';

         E(cutoff)=sum(sum(W.*abs((result*A-eye(len(1))).^2)));
      end
   end
   Emin=min(find(E==min(E)));
   S(Emin,Emin)./S(1,1);
end

Sinv=zeros(len(1),len(2));
for cutoff=1:Emin,
   Sinv(cutoff,cutoff)=1/S(cutoff,cutoff);
end
result=V*Sinv*U';








