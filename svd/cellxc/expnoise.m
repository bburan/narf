% function [asympxc,meanxc,pprset,bfitgood]=expnoise(prset,totlen,frac,datainv[=1],PMIN[=0.05]);
function [asympxc,meanxc,pprset,bfitgood]=expnoise(prset,totlen,frac,datainv,PMIN);

if ~exist('datainv','var'),
   datainv=1;
end

if ~exist('PMIN','var'),
   PMIN=0.05
end

DOMEDIAN=0;

goodidx=find(~isnan(prset(:,1)));
prset=prset(goodidx,:);
lens=totlen(goodidx);
ll=lens*frac;

% find asymptote for each cell.
pp=zeros(length(goodidx),2);
ppe=zeros(length(goodidx),1);
pprset=zeros(size(prset));
ppasymp=zeros(length(goodidx),1);
bfitgood=zeros(size(prset,1),1);
for cc=1:length(prset),
   f=1./(frac.*lens(cc));  % ie, x-axis is 1/T
   
   if sum(prset(cc,:)==0)>0,
      pp(cc,:)=polyfit(f,prset(cc,:),1);   % y axis is predxc
      pprset(cc,:)=(pp(cc,2)+pp(cc,1).*f);
   else
      pp(cc,:)=polyfit(f,1./prset(cc,:),1);   % y axis is predxc
      pprset(cc,:)=1./(pp(cc,2)+pp(cc,1).*f);
   end
   
   % error of linear fit, to determine if this cell qualifies
   if mean(prset(cc,:).^2)==0,
      ppe(cc)=1;
   else
      ppe(cc)=sqrt(mean((prset(cc,:)-pprset(cc,:)).^2) ./ ...
                   mean(prset(cc,:).^2));
   end
   ppasymp(cc)=pp(cc,2);
end
bfitgood(goodidx)=ppe<PMIN;
fprintf('%d\n',sum(bfitgood));
fitidx=find(bfitgood(goodidx));

% clear current axes for plotting results
hold off;
cla;

if datainv,
   % find asymptotic pred for 1./mean across all cells.
   
   f=1./(frac.*mean(lens));
   pp=polyfit(f,1./mean(prset(fitidx,:),1),1);
   ppmean=1./pp;
   
   pm=polyfit(f,1./median(prset(fitidx,:),1),1);
   ppmed=1./pm;
   
   plot([f 0],pp(2)+pp(1).*[f 0],'r--');
   hold on
   plot(f,1./mean(prset(fitidx,:)),'bx');
   errorbar(f,1./mean(prset(fitidx,:),1),...
            1./std(prset(fitidx,:),0,1)./sqrt(length(fitidx)),'bx');
   plot(0,pp(2),'rx');
   hold off
   
elseif datainv==0
   % find asymptotic pred for mean across all cell
   
   f=1./(frac.*mean(lens));
   pp=polyfit(f,mean(prset(fitidx,:),1),1);
   ppmean=pp;
   
   pm=polyfit(f,median(prset(fitidx,:),1),1);
   ppmed=pm;
   
   plot([f 0],pp(2)+pp(1).*[f 0],'r--');
   hold on
   plot(f,mean(prset(fitidx,:)),'bx');
   errorbar(f,mean(prset(fitidx,:),1),...
            std(prset(fitidx,:),0,1)./sqrt(length(fitidx)),'bx');
   plot(0,pp(2),'rx');
   hold off
end


if DOMEDIAN,
   asympxc=sqrt(ppmed(2));
   meanxc=sqrt(median(prset(fitidx,end)));
else
   asympxc=sqrt(ppmean(2));
   meanxc=sqrt(mean(prset(fitidx,end)));
end
   
   
