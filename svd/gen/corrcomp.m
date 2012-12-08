% function [p,m]=corrcomp(da,db,na,nb,axisrange,catidx,colorset,fillset);
%
% da, db= vectors, same size
% catidx = category index of each pair in (da,db) (catidx==0 skipped)
%
% na={name_a1,name_a2,...} array of strings, same size diff entries
% in catidx
%
% axisrange=[mina maxa minb maxb]
%
% created SVD 6/01
% modified SVD 2/02 - added multiple category support
%
function [p,m]=corrcomp(da,db,na,nb,axisrange,catidx,colorset,fillset);

if ~exist('na','var'),
   na='set 1';
end
if ~exist('nb','var'),
   nb='set 1';
end
if ~exist('axisrange','var') | isempty(axisrange),
   axisrange=[min(da(:)) max(da(:)) min(db(:)) max(db(:))];
end

if ~exist('colorset','var'),
   colorset=[0 0 0; 0 0 0; 0 0 0; 0 0 0];
end
if ~exist('fillset','var'),
   fillset=[0 0 0; 1 1 1; 0.67 0.67 0.67; 0.33 0.33 0.33];
end

if ~exist('catidx'),
   catidx=ones(size(da));
end

catidx(find(isnan(da) | isnan(db)))=0;

TTAIL=1;
TREP=2000;
RAD=25;
GCOUNT=max(catidx);
r=zeros(GCOUNT+1,1);
p=zeros(GCOUNT+1,1);
m=zeros(GCOUNT+1,1);
lf=zeros(GCOUNT+1,2);

% test all
[r(GCOUNT+1),exy,tt,p(GCOUNT+1)]=...
    randxcov(da(find(catidx)),db(find(catidx)),0,400);
%[p(GCOUNT+1),m(GCOUNT+1)]=randpairtest(da(find(catidx)),db(find(catidx)),TREP,TTAIL);
stitle=sprintf('%s v %s r=%.2f(%.3f)/',na,nb,r(GCOUNT+1),p(GCOUNT+1));
sx=sprintf('%s (%.2f/',na,mean(da(find(catidx))));
sy=sprintf('%s (%.2f/',nb,mean(db(find(catidx))));

cla
hold on;
for ii=1:GCOUNT,
   
   hp=scatter(da(find(catidx==ii)),db(find(catidx==ii)),RAD);
   set(hp,'MarkerFaceColor',fillset(ii,:),'MarkerEdgeColor',colorset(ii,:));
   
   [r(ii),exy,tt,p(ii)]=...
       randxcov(da(find(catidx==ii)),db(find(catidx==ii)),0,400);
   lf(ii,:)=polyfit(da(find(catidx==ii)),db(find(catidx==ii)),1);
   %[p(ii),m(ii)]=randpairtest(da(find(catidx==ii)),db(find(catidx==ii)),TREP,TTAIL);
   
   if ii<=1,
      stitle=[stitle,sprintf('%.2f(%.3f) ',r(ii),p(ii))];
   end
   if length(find(catidx==ii))>0,
      sx=[sx,sprintf('%.2f ',mean(da(find(catidx==ii))))];
      sy=[sy,sprintf('%.2f ',mean(db(find(catidx==ii))))];
   else
      sx=[sx,'0.00 '];
      sy=[sy,'0.00 '];
   end
end
sx(end)=')';
sy(end)=')';

if ~exist('axisrange'),
   minx=min([0; da(:); db(:)]);
   maxx=max([0.5; da(:); db(:)]);
   axisrange=[minx maxx minx maxx];
end

xx=[axisrange(1) axisrange(2)];
yy=xx*lf(1,1)+lf(1,2);
%hl=plot([axisrange(1) axisrange(2)],[axisrange(3) axisrange(4)],'k--');
hl=plot(xx,yy,'k--');
set(hl,'LineWidth',1);
hold off;
axis(axisrange);
axis square
box on

ht=title(stitle);
hx=xlabel(sx);
hy=ylabel(sy);
set(ht,'FontSize',8);
set(hx,'FontSize',8);
set(hy,'FontSize',8);

return
if exist('d2a') & ~isempty(d2a),
   p2=randpairtest(d2a,d2b,TREP);
   ht=title(sprintf(['%s v %s p<%.3f n=%d (red %s,%s p<%.3f)'],...
                    n1a,n1b,p1,length(d1a),n2a,n2b,p2));
   hx=xlabel(sprintf('%s m=%.2f (red %s %.2f)',...
                     n1a,mean(d1a),n2a,mean(d2a)));
   hy=ylabel(sprintf('%s m=%.2f (red %s %.2f)',...
                     n1b,mean(d1b),n2b,mean(d2b)));
else
   p2=1;
   m2=0;
   ht=title(sprintf('%s v %s m=%.2f p<%.3f n=%d',n1a,n1b,m1,p1,length(d1a)));
   hx=xlabel(sprintf('%s (m=%.2f)',n1a,mean(d1a)));
   hy=ylabel(sprintf('%s (m=%.2f)',n1b,mean(d1b)));
end

set(ht,'FontSize',8);
set(hx,'FontSize',8);
set(hy,'FontSize',8);
set(gca','FontSize',8);



