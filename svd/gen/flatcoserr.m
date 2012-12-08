function chi2=flatcoserr(beta,x,y,e)

ya=flatcos(x,beta);
chi2=sum((ya-y).^2 ./ e.^2);

