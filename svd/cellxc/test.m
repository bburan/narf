function ppbetter=test(pp,T,prset,prerr);

ppbetter=fsolve('invlineerr',pp,optimset('Display','iter'),T,prset,prerr);

return


function F = myfun(x)
F = sin(x);
 

%function err=invlineerr(x,p1,p2,p3);
%
%ypred=1./ (x(1).*(1./p1) + x(2));
%
%err=sqrt(mean( (abs(p2-ypred)./p3).^2 ));

