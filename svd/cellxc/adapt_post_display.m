function adapt_post_display(resp,pred);

p1=pred;
p2=conv2(pred(:),[0 0 0 0 1/3 1/3 1/3]','same');
%p2=conv2(pred(:),[0 0 0 1/2 1/2]','same');

ggidx=find(~isnan(resp(:)+p1(:)+p2(:)));
resp=resp(ggidx);
p1=p1(ggidx);
p2=p2(ggidx);
pred=pred(ggidx);

idx=round(linspace(1,length(resp),1000));
scatter(p1(idx),p2(idx),round(resp(idx).*30+1));
drawnow


beta=fit_adapt_post(pred,resp);
zz0=adapt_post(beta,pred);


[xx,yy]=meshgrid(linspace(min(p1),max(p1),50),...
                 linspace(min(p2),max(p2),50));

dd=yy(2)-yy(1);
zz=zeros(size(xx));
zz2=zeros(size(xx));
for ii=1:length(zz(:)),
   xi=xx(ii);
   yi=yy(ii);
   di=sqrt((p1-xi).^2 + (p2-yi).^2);
   di(di<dd./4)=dd./4;
   didx=find(di<dd.*2);
   di=1./di(didx);
   if length(di)>1,
      zz(ii)=sum(resp(didx).*di)./sum(di);
      zz2(ii)=sum(zz0(didx).*di)./sum(di);
   end
end

imagesc(xx(1,:)',yy(:,1),zz);
hold on
contour(xx,yy,zz,[0.01 0.01],'w:');
plot([min(p1) max(p1)],[min(p1),max(p1)],'w--');
contour(xx,yy,zz2,[0.1 0.1].*max(zz0),'b');
contour(xx,yy,zz2,[0.3 0.3].*max(zz0),'g');
contour(xx,yy,zz2,[0.5 0.5].*max(zz0),'r');
contour(xx,yy,zz2,[0.7 0.7].*max(zz0),'y');
hold off

axis xy
colormap(gray)
drawnow

[xcov(pred,resp,0,'coeff') xcov(zz0,resp,0,'coeff')]




keyboard

return


idx=round(linspace(1,length(resp),400));

figure(1);
clf
scatter(p1(idx),p2(idx),round(resp(idx).*30+1));
