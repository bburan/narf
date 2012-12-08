% function ypred=invlineerr(beta,T);
%
% ypred=1./ (beta(1).*(1./T) + 1./beta(2))
%
%
function ypred=invline(beta,T);

ypred=1./ (beta(1).*(1./T) + beta(2));

