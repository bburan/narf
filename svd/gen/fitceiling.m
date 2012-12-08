% function [yinf,yinferr,p]=fitceiling(x,y);
function [yinf,yinferr,p]=fitceiling(x,y);

y(find(y==0))=nan;

aa=1./x;
if size(y,1)>1,
   bb=nanmean(1./y);
else
   bb=1./y;
end

p0 = polyfit(aa,bb,1);
%1./p0(2)

x=repmat(x,size(y,1),1);
gidx=find(~isnan(y));

pmin=[0 -0.001];
p0(p0<pmin)=pmin(p0<pmin);

fitopt=optimset('Display','off');
p=lsqcurvefit('ceilfun',p0,x(gidx),y(gidx),...
              pmin,[inf inf],fitopt);

if p(2)<0.25, % bad fit!
   p(:)=nan;
end

yinf=1./p(2);
yinferr=0;

if 0,
   figure(1);
   clf
   subplot(1,2,1);
   plot(aa,bb);
   aa0=linspace(0.05,1.5);
   bb0=mean((p(2)*ones(size(aa0))+p(1)*(1./aa0)),1);
   hold on
   plot(1./aa0,1./ceilfun(p,aa0),'--');
   hold off
   
   subplot(1,2,2);
   plot(x,y);
   hold on
   plot(aa0,ceilfun(p,aa0),'--');
   %plot(aa0,ceilfun(p0,aa0),':');
   hold off
end

%fprintf('p=(%.3f, %.3f) y: %.3f --> %.3f',p,mean(y(:,end)),1./p(2));

return


y(find(y==0))=nan;

aa=1./x;
bb=1./y;

p=zeros(size(y,1),2);
for ii=1:size(y,1),
   
   if sum(isnan(bb(ii,:)))==0,
      
      fitopt=optimset('Display','off');
      p0 = polyfit(aa,bb(ii,:),1);
      if p0(1)<0,
         p0(1)=-p0(1);
      end
      
      p(ii,:)=lsqcurvefit('ceilfun',p0,x,y(ii,:),...
                          [-inf -inf],[inf inf],fitopt);
      
      if 1,
         figure(1);
         clf
         subplot(1,2,1);
         plot(aa,bb);
         aa0=linspace(0.05,1.5);
         bb0=mean((p(ii,2)*ones(size(aa0))+p(ii,1)*(1./aa0)),1);
         hold on
         plot(1./aa0,1./ceilfun(p(ii,:),aa0),'--');
         hold off
         
         subplot(1,2,2);
         plot(x,y);
         hold on
         plot(aa0,ceilfun(p(ii,:),aa0),'--');
         %plot(aa0,ceilfun(p0,aa0),':');
         hold off
      end
      
      if p(ii,2)<0.0, % bad fit!
         p(ii,:)=nan;
      end
   end
end

nbad=sum(isnan(p(:,2)));
if nbad<size(p,1)/2,
   [dd,pidx]=sort(-p(:,2));
   p(pidx(1:nbad),:)=nan;
end

yinf=nanmean(1./p(:,2));
yinferr=nanstd(1./p(:,2));

fprintf('p=(%.3f, %.3f) y: %.3f --> %.3f\n',p(1,:),mean(y(:,end)),yinf);
pause

