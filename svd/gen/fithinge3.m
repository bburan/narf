% function [hingeparms,beta0]=fithinge3(predpsth,resp,showfit[=0],phasecount[=1],forcevalues,beta0);
%
% hingparms=[x01 slope1 slope2 ...  slope<phasecount>]
%
% y= (x-x0) * slopeP  * (x>x0)  summed over each phase channel
%
function [hingeparms,beta0]=fithinge3(predpsth,resp,showfit,phasecount,forcevalues);

if ~exist('showfit','var'),
   showfit=0;
end
if ~exist('phasecount','var'),
   phasecount=1;
end
if ~exist('forcevalues','var'),
   forcevalues=[nan];
end
if ~exist('beta0','var'),
   beta0=zeros(1,phasecount+1);
   b0given=0;
else
   b0given=1;
end

spacecount=size(predpsth,2);
if spacecount>phasecount,
   tpred=zeros(size(predpsth,1),phasecount);
   stepsize=ceil(spacecount/phasecount);
   for phaseidx=1:phasecount,
      
      pidx=(1:stepsize)+stepsize*(phaseidx-1);
      pidx=pidx(find(pidx<=spacecount));
      
      tpred(:,phaseidx)=sum(predpsth(:,pidx),2);
   end
   predpsth=tpred;
end

lb=-inf*ones(1,phasecount+1);
ub=inf*ones(1,phasecount+1);

for phaseidx=1:phasecount,
   tx=predpsth(:,phaseidx);
   r0=resp-mean(resp);
   r1=tx-mean(tx);
   d1=sum(r1.^2);
   if d1>0,
      scf=sum(r0.*r1)./d1;
   else
      scf=(max(resp)-min(resp))./(max(tx)-min(tx));
   end
   
   if length(forcevalues)>=phaseidx+1 & ~isnan(forcevalues(phaseidx+1)),
      beta0(phaseidx+1)=forcevalues(phaseidx+1);
   else
      beta0(phaseidx+1)=scf;
   end
end

% threshold - choose after
tx=sum(predpsth,2);
beta0(1)=min(resp)+std(resp);


fitopt=optimset('Display','off');

hingeparms=lsqcurvefit('hinge3',beta0,predpsth,resp,lb,ub,fitopt,forcevalues);
hingeparms(find(~isnan(forcevalues)))=forcevalues(find(~isnan(forcevalues)));

%[beta0; lb; ub; hingeparms ]

if showfit,
   for phaseidx=1:phasecount,
      subplot(phasecount,1,phaseidx);
      cla
      scatter(predpsth(:,phaseidx),resp);
      xx=linspace(min(predpsth(:)),max(predpsth(:)),50);
      hold on
      plot(xx,hinge2(beta0(phaseidx*2+(-1:0)),xx),'r--');
      plot(xx,hinge2(hingeparms(phaseidx*2+(-1:0)),xx));
      hold off
      drawnow
   end
end


