% function [m]=histfracs(d1,d2,n1,sunits,axisrange,n);
%
% axisrange=[mina maxa minb maxb]
%
function [m]=histfracs(d1,d2,n1,sunits,axisrange,n);

if ~exist('n1','var'),
   n1='data';
end
if ~exist('sunits','var'),
   sunits='';
end
if ~exist('n','var'),
   n=10;
end

cla
%[p,m]=randpairtest(d1,zeros(size(d1)),5000);
mall=median([d1(:);d2(:)]);
m1=median(d1(:));
if length(d2(:))>0
   m2=median(d2(:));
else
   m2=0;
end

if nargout>0,
   m=[mall;m1;m2];
end

if ~exist('axisrange','var'),
   [xx,histrange]=hist([d1(:);d2(:)],n);
   axisrange=[histrange(1) histrange(end)];
else
   histrange=linspace(axisrange(1),axisrange(2),n);
end
if isempty(d1),
   data=zeros(n,1);
else
   data=reshape(hist(d1,histrange),length(histrange),1);
end
if isempty(d2),
   data=[data zeros(n,1)];
else
   data=[data reshape(hist(d2,histrange),length(histrange),1)];
end

hb=bar(histrange,data,'stacked');

ediff=(histrange(2)-histrange(1))/2;
if length(axisrange)==4,
   axis([histrange(1)-ediff histrange(end)+ediff axisrange(3:4)]);
else
   a=axis;
   axis([histrange(1)-ediff histrange(end)+ediff a(3) a(4)]);
end
axis square

ht=title(sprintf('%s n=%d,%d m=%.2f/%.2f',n1,...
                 length(d1)+length(d2),length(d2),mall,m1));
hx=xlabel(sprintf('%s',sunits));

set(ht,'FontSize',8);
set(hx,'FontSize',8);
%set(hy,'FontSize',8);
set(gca','FontSize',8);


