%function [parms,res]=findthresh(pred,resp,showfit,thcount,rescount);
function [parms,res]=findthresh(pred,resp,showfit,thcount,rescount);

if ~exist('showfit','var'),
   showfit=0;
end
if ~exist('thcount','var'),
   thcount=20;
end
if ~exist('rescount','var'),
   rescount=2;
end

if std(resp)==0 || std(pred(:))==0,
   parms=0;
   res(1).thrange=0;
   res(1).xc=0;
   res(1).bestthxidx=1;
   return
end

pmin=min(pred(:));
pmax=max(pred(:));

pmax=(pmax-pmin)*1/2 + pmin;
%pmax=(pmax-pmin)*3/5 + pmin;
%pmax=(pmax-pmin)*4/5 + pmin;

for residx=1:rescount,
   
   thrange=linspace(pmin,pmax,thcount);
   
   xc=zeros(thcount,1);
   rt=pred;
   for thidx=1:length(thrange),
      rt(rt<thrange(thidx))=0;
      %rt=thresh(thrange(thidx),pred);
      if std(rt)>0,
         xc(thidx)=xcov(rt,resp,0,'coeff');
      end
   end
   
   bestthidx=min(find(xc==max(xc)));
   if isempty(bestthidx),
      bestthidx=1;
   end
   
   %fprintf('iteration %d: pmin=%.3f pmax=%.3f  thbest=%d (%.3f)\n',...
   %        residx,pmin,pmax,bestthidx,thrange(bestthidx));
   
   res(residx).thrange=thrange;
   res(residx).xc=xc;
   res(residx).bestthidx=bestthidx;
   
   %figure(1);
   %clf
   %plot(thrange,xc);
   %drawnow
   
   if bestthidx>1,
      pmin=thrange(bestthidx-1);
   else
      pmin=thrange(1);
   end
   if bestthidx<thcount,
      pmax=thrange(bestthidx+1);
   else
      pmax=thrange(thcount);
   end
   
end

parms=thrange(bestthidx);

return

beta0=[min(pred)+std(pred)];
lb=   [min(pred)          ];
ub=   [max(pred)-2*std(pred)];

fitopt=optimset('Display','off');

if find(abs(ub-lb)==0),
   parms=[0];
   res=[];
   return
end

parms=lsqcurvefit('thresh',beta0,pred,resp,lb,ub,fitopt);
res=[];

if showfit,
   clf
   
   scatter(pred(1:100),resp(1:100));
   xx=linspace(min(pred),max(pred),50);
   hold on
   plot(xx,thresh(beta0,xx),'r--');
   plot(xx,thresh(parms,xx));
   hold off
   drawnow
end


return

