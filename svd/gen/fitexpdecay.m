% function [beta,beta0]=fitexpdecay(x,y,verbose);
%
function [beta,beta0]=fitexpdecay(x,y,verbose);

if ~exist('e'),
   e=[];
end
if ~exist('beta0'),
   beta0=[];
end
if ~exist('verbose','var'),
   verbose=0;
end

opt=optimset('Display','off');

% beta=[tau,A,offset]
[x,xi]=sort(x);
y=y(xi);
beta0=[mean(diff(x)).*3 mean(y(1:2))-mean(y((end-1):end)) ...
       mean(y((end-1):end))];
betalb=[min(diff(x))./100 -Inf -Inf];
betaub=[max(diff(x)).*100 Inf Inf];

beta=lsqcurvefit('expdecay',...
                 beta0,x,y,betalb,betaub,opt);

beta=beta(:);
beta0=beta0(:);

if verbose,
   % display fit
   figure(1)
   clf
   plot(x,y);
   hold on
   plot(x,expdecay(beta,x),'r--');
   plot(x,expdecay(beta0,x),'k:');
   hold off
   legend('actual','fit','fit0');
end


return

