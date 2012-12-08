% function [parms]=fitexpthresh(pred,resp,showfit[=0]);
%
% parms=[thresh expon scf];
%
function [parms]=fitexpthresh(pred,resp,showfit);

if ~exist('showfit','var'),
   showfit=0;
end

beta0=[min(pred)+std(pred)   1.3  1.0];
lb=   [min(pred)             0.8  0.01];
ub=   [max(pred)-std(pred)   3    inf];

fitopt=optimset('Display','off');

if find(abs(ub-lb)==0),
   parms=[0 0];
   return
end

parms=lsqcurvefit('expthresh',beta0,pred,resp,lb,ub,fitopt);

if showfit,
   clf
   
   scatter(pred(1:100),resp(1:100));
   xx=linspace(min(pred),max(pred),50);
   hold on
   plot(xx,expthresh(beta0,xx),'r--');
   plot(xx,expthresh(parms,xx));
   hold off
   drawnow
end

