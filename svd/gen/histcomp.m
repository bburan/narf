% function [p1,p2]=histcomp(d1,d2,n1a,n1b,sunits,axisrange,DOMEDIAN);
%
% axisrange=[mina maxa minb maxb]
%
function [p,m]=histcomp(d1,d2,n1a,n1b,sunits,axisrange,DOMEDIAN);

if ~exist('axisrange','var'),
   diffmax=max(abs([d1(:); d2(:)]));
   axisrange=[-diffmax diffmax];
end
if ~exist('sunits','var'),
   sunits=[];
end
if ~exist('DOMEDIAN','var'),
   DOMEDIAN=0;
end


cla
[p,m]=randpairtest(d1,zeros(size(d1)),5000);
if DOMEDIAN,
   mall=median([d1(:);d2(:)]);
   m1=median(d1(:));
   if length(d2(:))>0
      m2=median(d2(:));
   else
      m2=0;
   end
   
else
   mall=mean([d1(:);d2(:)]);
   m1=mean(d1(:));
   if length(d2(:))>0
      m2=mean(d2(:));
   else
      m2=0;
   end
end

diffrange=linspace(axisrange(1),axisrange(2),13);
data(1).d=d1;
data(2).d=d2;
pophist(data,diffrange);

ht=title(sprintf('%s v %s (n=%d,%d m=%.2f/%.2f p<%.3f)',n1a,n1b,...
                 length(d1),length(d2),mall,m1,p));
hx=xlabel(sprintf('diff (%s)',sunits));

set(ht,'FontSize',8);
set(hx,'FontSize',8);
%set(hy,'FontSize',8);
set(gca','FontSize',8);

if length(axisrange)==4,
   axis(axisrange);
end
axis square

