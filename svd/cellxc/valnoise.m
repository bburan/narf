% function [ccr,ccrM,yinf,n,y,p]=valnoise(pred,r2);
%
% to get expected asymptotic cc^2 from single trial cc:
%
% ccexp=ccp./ccr;
%
% (ccp = single trial cc^2)
% (ccexp = asymp cc^2)
%
% ccr - ratio of single trial cc^2 to infinite trial cc^2
% ccrM - ratio of single trial cc^2 to M trial cc^2
%
function [ccr,ccrM,yinf,n,y,p]=valnoise(pred,r2);

gidx=find(~isnan(pred) & ~isnan(r2(:,1)));
pred=pred(gidx);
r2=r2(gidx,:);

M=size(r2,2);

fracs=[0 0.05 0.2 0.75];
edges=ceil(M.*cumsum(fracs));
if edges(3)==2,
   edges(3)=3;
end
ccr0=zeros(M,length(fracs)-1);
xc1=zeros(M,1);
if var(pred)>0,
   xcM=xcov(mean(r2,2),pred,0,'coeff').^2;
else
   xcM=0;
end

ra=zeros(size(r2,1),length(fracs)-1);

for mm=1:M
   for ii=1:length(fracs)-1,
      if edges(ii+1)<edges(ii)+1,
         warning('not enough trials to measure ceiling');
         [M edges]
      end
      ra(:,ii)=mean(r2(:,edges(ii)+1:edges(ii+1)),2);
      if std(ra(:,ii))>0 & std(pred)>0,
         ccr0(mm,ii)=xcov(ra(:,ii),pred,0,'coeff').^2;
      end
   end
   r2=shift(r2,[0 1]);
   
   if std(r2(:,1))>0 & std(pred)>0,
      xc1(mm)=xcov(r2(:,1),pred,0,'coeff').^2;
   end
end

y=ccr0;
x=diff(edges);

if 0,
   [yinf,yinferr,p]=fitceiling(x,y);
else
   yinf=zeros(size(y,1),1).*nan;
   yinferr=zeros(size(y,1),1).*nan;
   p=zeros(size(y,1),2).*nan;
   for ii=1:size(y,1),
      if ~sum(y(ii,:)==0),
         
         [yinf(ii),yinferr(ii),p(ii,:)]=fitceiling(x,y(ii,:));
      end
   end
   yinf=nanmean(yinf);
   yinferr=nanstd(yinf).*sqrt(sum(~isnan(yinf))-1);
end


ccr=mean(xc1)./yinf;
ccrM=xcM./yinf;
ccre=mean(xc1)./(yinf-yinferr) - ccr;

n=diff(edges);
%y=mean(ccr0,1);
ye=std(ccr0).*sqrt(M-1);


if 0,
   clf
   errorbar(x,y,ye);
   hold on
   for ii=1:size(p,1);
   aa=1:20;
   bb=1./(p(ii,2)+p(ii,1)./aa);
   plot(aa,bb,'k:')
   end
   hold off
end



return


ccrM=(1+sqrt(1./mean(ccr0))) ./ (-M + M * sqrt(1./mean(ccr0)) + 2);
ccr=2./(2-M+M*sqrt(1./ccr0));

ccre=std(ccr,1).*sqrt(M-1);
ccr=mean(ccr);



