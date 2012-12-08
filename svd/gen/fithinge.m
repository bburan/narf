% function [hingeparms]=fithinge(predpsth,resp,showfit[=0],forcevalues);
%
% hingparms=[x0 y0 slope]
%
% y= (x-x0) * slope  * (x>x0) + y0
%
function [hingeparms]=fithinge(predpsth,resp,showfit,forcevalues);

if ~exist('showfit','var'),
   showfit=0;
end
if ~exist('forcevalues','var'),
   forcevalues=[nan nan nan];
end

beta0=[median(predpsth) min(resp)+(max(resp)-min(resp))/10 ...
       (max(resp)-min(resp))./(max(predpsth)-min(predpsth))];
fitopt=optimset('Display','off');
lb=[min(predpsth);                    min(resp);               0  ];
ub=[max(predpsth)/2+min(predpsth)/2;  max(resp)/2+min(resp)/2; inf];

if find(abs(ub-lb)==0),
   hingeparms=[0 0 0];
   return
end
   
hingeparms=lsqcurvefit('hinge',beta0,predpsth,resp,lb,ub,fitopt,forcevalues);
hingeparms(find(~isnan(forcevalues)))=forcevalues(find(~isnan(forcevalues)));

if showfit,
   clf
   scatter(predpsth(round(linspace(1,length(predpsth)))),resp(round(linspace(1,length(resp)))));
   xx=linspace(min(predpsth),max(predpsth),50);
   hold on
   plot(xx,hinge(beta0,xx),'r--');
   plot(xx,hinge(hingeparms,xx));
   hold off
   drawnow
end

