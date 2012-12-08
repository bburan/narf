% function err=invlineerr(beta,T,y,yerr);
%
% ypred=1./ (beta(1).*(1./T) + beta(2))
%
% err=sum( abs(y-ypred)./yerr );
%
function err=invlineerr(beta,T,y,yerr);

ypred=1./ (beta(1).*(1./T) + beta(2));

err=sqrt(mean( (abs(y-ypred)./y).^2 ));
