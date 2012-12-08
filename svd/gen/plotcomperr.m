% function p=plotcomperr(d1,d2,e1,e2,n1,n2,axisrange);
%
% axisrange=[mina maxa minb maxb]
%
function p=plotcomperr(d1,d2,e1,e2,n1,n2,axisrange);

TREP=2000;
if ~exist('axisrange'),
   minx=min([0; d1(:); d2(:)]);
   maxx=max([0.5; d1(:); d2(:)]);
   axisrange=[minx maxx minx maxx];
end

cla
hold on
if sum(e1)==0,
   scatter(d1,d2,20,'k','filled');
else
   for ii=1:length(d1),
      line([d1(ii)-e1(ii) d1(ii)+e1(ii)],[d2(ii) d2(ii)]);
      line([d1(ii) d1(ii)],[d2(ii)-e2(ii) d2(ii)+e2(ii)]);
   end
end

hl=plot([axisrange(1) axisrange(2)],[axisrange(3) axisrange(4)],'k--');
set(hl,'LineWidth',1);

hold off;
axis(axisrange);
axis square
box on

[p,m]=randpairtest(d1,d2,TREP);
ht=title(sprintf('%s v %s m=%.2f p<%.3f n=%d',n1,n2,m,p,length(d1)));
hx=xlabel(sprintf('%s (m=%.2f)',n1,mean(d1)));
hy=ylabel(sprintf('%s (m=%.2f)',n2,mean(d2)));

set(ht,'FontSize',8);
set(hx,'FontSize',8);
set(hy,'FontSize',8);
set(gca','FontSize',8);






