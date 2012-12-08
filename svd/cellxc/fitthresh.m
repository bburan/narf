% function [threshparms]=fitthresh(predpsth,resp,showfit[=0],phasecount[=1]);
%
%
function [threshparms]=fitthresh(predpsth,resp,showfit,phasecount);

if ~exist('showfit','var'),
   showfit=0;
end
if ~exist('phasecount','var'),
   phasecount=1;
end

spacecount=size(predpsth,2);
if spacecount<phasecount,
   disp('spacecount<phasecount, reducing phasecount to spacecount');
   phasecount=spacecount;
end

beta0=zeros(1,phasecount);
lb=zeros(1,phasecount);
ub=zeros(1,phasecount);
stepsize=ceil(spacecount/phasecount);
for phaseidx=1:phasecount,
   pidx=(1:stepsize)+stepsize*(phaseidx-1);
   pidx=pidx(find(pidx<=spacecount));
   tx=sum(predpsth(:,pidx),2);
   beta0(phaseidx)=min(tx)+std(tx);
   lb(phaseidx)=0;
   ub(phaseidx)=median(tx)+std(tx)*1.5;
   if ub(phaseidx)-lb(phaseidx)==0,
      ub(phaseidx)=lb(phaseidx)+0.0000001;
   end
end


fitopt=optimset('Display','off');

threshparms=lsqcurvefit('thresh',beta0,predpsth,resp,lb,ub,fitopt);

%[beta0; threshparms; lb; ub]

if showfit,
   for phaseidx=1:phasecount,
      subplot(floor(sqrt(phasecount)),...
              ceil(phasecount/floor(sqrt(phasecount))),phaseidx,'.');
      cla
      pidx=(1:stepsize)+stepsize*(phaseidx-1);
      pidx=pidx(find(pidx<=spacecount));
      scatter(sum(predpsth(:,pidx),2),resp);
   end
   drawnow
end
