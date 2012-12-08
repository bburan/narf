% function S=sparseness(r)
%
function S=sparseness(r)

m=mean(r);
sig=std(r);
n=length(r);

S=(1-m./(m+sig)) ./ (1 - 1./n);