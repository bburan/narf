% function [p,m]=plotcomp(da,db,na,nb,axisrange,catidx,stat,colorset,fillset);
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
function [p,m]=plotcomp(da,db,na,nb,axisrange,catidx,stat,colorset,fillset);

if ~exist('na','var'),
   na='set 1';
end
if ~exist('nb','var'),
   nb='set 1';
end
if ~exist('axisrange','var') | isempty(axisrange),
   axisrange=[min([da(:);db(:)]) max([da(:);db(:)]) ...
              min([da(:);db(:)]) max([da(:);db(:)])];
end
if ~exist('catidx'),
   catidx=ones(size(da));
end
if ~exist('stat','var'),
   %stat='median';
   stat='mean';
end
if ~exist('colorset','var'),
   colorset=[0 0 0; 0 0 0; 0 0 0; 0 0 0];
end
if ~exist('fillset','var'),
   fillset=[0 0 0; 1 1 1; 0.67 0.67 0.67; 0.33 0.33 0.33];
end

% skip any nan entries
catidx(find(isnan(da) | isnan(db)))=0;

TTAIL=1;
TREP=2000;
RAD=20;
GCOUNT=max(catidx);
p=zeros(GCOUNT+1,1);
m=zeros(GCOUNT+1,1);

% test all
[p(GCOUNT+1),m(GCOUNT+1)]=...
    randpairtest(da(find(catidx)),db(find(catidx)),TREP,TTAIL,stat);
stitle=sprintf('%s v %s %.3f/',na,nb,p(GCOUNT+1));

sx=sprintf('%s (%.2f/',na,feval(stat,da(find(catidx))));
sy=sprintf('%s (%.2f/',nb,feval(stat,db(find(catidx))));

cla
hold on;
for ii=GCOUNT:-1:1,
   
   hp=scatter(da(find(catidx==ii)),db(find(catidx==ii)),RAD);
   set(hp,'MarkerFaceColor',fillset(ii,:),'MarkerEdgeColor',colorset(ii,:));
   
   [p(ii),m(ii)]=randpairtest(da(find(catidx==ii)),...
                              db(find(catidx==ii)),TREP,TTAIL,stat);
   
   if ii<=1,
      stitle=[stitle,sprintf('%.3f ',p(ii))];
   end
   if length(find(catidx==ii))>0,
      sx=[sx,sprintf('%.2f ',feval(stat,da(find(catidx==ii))))];
      sy=[sy,sprintf('%.2f ',feval(stat,db(find(catidx==ii))))];
   else
      sx=[sx,'0.00 '];
      sy=[sy,'0.00 '];
   end
end
sx(end)=')';
sy(end)=')';

if ~exist('axisrange') || isempty(axisrange),
   minx=min([0; da(:); db(:)]);
   maxx=max([0.5; da(:); db(:)]);
   axisrange=[minx maxx minx maxx];
end

hl=plot([axisrange(1) axisrange(2)],[axisrange(3) axisrange(4)],'k--');
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




