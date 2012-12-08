% function pophist(data,edges,p);
%
% data(1...n).d = vector of point to hist in each group
%
function pophist(data,edges,p);

n=[];
ediff=(edges(2)-edges(1))/2;

for ii=1:length(data),
   if isempty(data(ii).d),
      n=[n zeros(length(edges),1)];
   else
      n=[n reshape(histc(data(ii).d,edges),length(edges),1)];
   end
end

%hb=bar(edges+ediff,n,'stacked');
hb=bar(edges+ediff,n);
set(hb,'barwidth',1);

if exist('p','var'),
   a=axis;
   nmax=max(sum(n,2));
   m=mean(data(1).d);
   r=max(edges)-min(edges);
   hold on
   plot(m,nmax+1,'x');
   ht=text(m+r/20,nmax+1,sprintf('p<%.3f',max([p 0.001])));
   hold off
end

