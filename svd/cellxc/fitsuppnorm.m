% function [parms]=fitexpthresh(pred,resp,showfit[=0]);
%
% parms=[thresh expon scf];
%
function [parms]=fitsuppnorm(pred,resp,showfit);

if ~exist('showfit','var'),
   showfit=0;
end

beta0=[1 1-min(pred(:,2)) 0];
lb=   [-inf 0 -inf];
ub=   [inf inf inf];

fitopt=optimset('Display','off');

if find(abs(ub-lb)==0),
   parms=[0 0];
   return
end

parms=lsqcurvefit('suppnorm',beta0,pred,resp,lb,ub,fitopt);


if showfit,
   [xcov(pred(:,1)-pred(:,2),resp,0,'coeff'),...
    xcov(suppnorm(parms,pred),resp,0,'coeff')]
   clf
   
   scatter(pred(1:100),resp(1:100));
   xx=linspace(min(pred),max(pred),50);
   hold on
   plot(xx,expthresh(beta0,xx),'r--');
   plot(xx,expthresh(parms,xx));
   hold off
   drawnow
end

