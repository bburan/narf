%function [parms,res]=findthreshsat(pred,resp,showfit,thcount,rescount);
function [parms,res]=findthreshsat(pred,resp,showfit,thcount,rescount);

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
smax=max(pred(:));

smin=(smax-pmin)*0.7 + pmin;
pmax=(smax-pmin)*0.3 + pmin;

%pmax=(pmax-pmin)*3/5 + pmin;
%pmax=(pmax-pmin)*4/5 + pmin;

for residx=1:rescount,
   
   thrange=linspace(pmin,pmax,thcount);
   sarange=linspace(smin,smax,thcount);
   
   xc=zeros(thcount,2);
   rt=pred;
   rs=pred;
   for thidx=1:length(thrange),
      rt(rt<thrange(thidx))=0;
      if std(rt)>0,
         xc(thidx,1)=xcov(rt,resp,0,'coeff');
      end
      
      rs(rs>sarange(thidx))=sarange(thidx);
      if std(rs)>0,
         xc(thidx,2)=xcov(rs,resp,0,'coeff');
      end
   end
   
   bestthidx=min(find(xc(:,1)==max(xc(:,1))));
   if isempty(bestthidx),
      bestthidx=1;
   end
   bestsaidx=min(find(xc(:,2)==max(xc(:,2))));
   if isempty(bestsaidx),
      bestsaidx=1;
   end
   
   %fprintf('iteration %d: pmin=%.3f pmax=%.3f  thbest=%d (%.3f)\n',...
   %        residx,pmin,pmax,bestthidx,thrange(bestthidx));
   
   res(residx).xc=xc;
   res(residx).thrange=thrange;
   res(residx).bestthidx=bestthidx;
   res(residx).sarange=sarange;
   res(residx).bestsaidx=bestsaidx;
   
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
   
   if bestsaidx>1,
      smin=sarange(bestsaidx-1);
   else
      smin=sarange(1);
   end
   if bestsaidx<thcount,
      smax=sarange(bestsaidx+1);
   else
      smax=sarange(thcount);
   end
   
end

parms=[thrange(bestthidx);sarange(bestsaidx)];

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

