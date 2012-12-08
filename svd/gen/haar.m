% function psi=haar(t,u,s)
%
function psi=haar(t,u,j)

ts=(t-u)./2.^j;

psi=zeros(length(t),1);
psi(find(ts>=0 & ts<0.5))=1;
psi(find(ts>=0.5 & ts<1))=-1;
psi=psi.*2^(-j);

