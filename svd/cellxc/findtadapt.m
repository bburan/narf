%function [parms,res]=findtadapt(pred,resp,showfit);
function [parms,res]=findtadapt(pred,resp,showfit);

if ~exist('showfit'),
   showfit=0;
end


% parms=[adaptstr tdecay]
beta0=[0.5  1 mean(pred)];
lb=   [0.01 0 min(pred)];
ub=   [0.99 2 max(pred)-std(pred)/10];

fitopt=optimset('Display','off');

parms=lsqcurvefit('tadapt',beta0,pred,resp,lb,ub,fitopt);
res=[];

if showfit,
   clf
   y=tadapt(parms,pred);
   
   plot(pred(1:100),resp(1:100),'gx');
   hold on
   plot(y(1:100),resp(1:100),'k.');
   hold off
   drawnow;
end


return

