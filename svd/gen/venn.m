% function venn(d,showlabels)
%
% plot a 2- or 3 circle venn diagram
%
% d is a 2 x 2 or 2 x 2 x 2 matrix, telling the number of points
% falling in each combination of out (1) and in (2) along each
% of the dimcount dimensions. plots on currently active axes
%
% created SVD 1/2/04
%
function venn(d,showlabels)

if ~exist('showlabels','var'),
   showlabels=0;
end

dimcount=length(size(d));

marg=zeros(dimcount,1);
olap=zeros(dimcount);
for ii=1:dimcount,
   notii=[1:ii-1 ii+1:dimcount];
   m=permute(d,[ii notii]);
   marg(ii)=sum(m(2,:));
   
   for jj=1:dimcount,
      if jj==ii,
         olap(ii,jj)=marg(ii);
      else
         m=permute(d,[ii jj setdiff(1:dimcount,[ii jj])]);
         olap(ii,jj)=sum(m(2,2,:));
      end
   end
end


rad=sqrt(marg./pi);
dist=zeros(dimcount);
optset=optimset('display','off');

for ii=1:dimcount-1,
   for jj=ii+1:dimcount,
      dist0=mean([rad(ii) rad(jj)]);
      r1=rad(ii);
      r2=rad(jj);
      dist(ii,jj)=lsqcurvefit('olaparea',dist0,0,olap(ii,jj),...
                              0,rad(ii)+rad(jj),optset,r1,r2);
      %[ r1 r2 olap(ii,jj) dist0 ...
      %  dist(ii,jj) olaparea(dist(ii,jj),0,r1,r2)]
   end
end

if dimcount==3,
   
   s=(dist(1,2)+dist(1,3)+dist(2,3))/2;
   a=sqrt(s*(s-dist(1,2))*(s-dist(1,3))*(s-dist(2,3)));
   
   %[dist(1,2) dist(1,3) dist(2,3)]
   y0=[0 0 a./dist(1,2).*2];
   x0=[0 dist(1,2) sqrt(dist(1,3).^2-y0(3).^2)];
   solox=x0+[-rad(1).*0.9 rad(2).*0.9 0];
   soloy=y0+[0 0 rad(3).*0.9];
   
   pairx=[mean(x0(2:3))-x0(1)/3 mean(x0([1 3]))-x0(2)/3 ...
          mean(x0(1:2))-x0(3)/3];
   pairy=[mean(y0(2:3))-y0(1)/3 mean(y0([1 3]))-y0(2)/3 ...
          mean(y0(1:2))-y0(3)/3];
else
   y0=[0 0];
   x0=[0 dist(1,2)];
   solox=x0+[-rad(1).*0.9 rad(2).*0.9];
   soloy=y0+[0 0];
end

pcol={'k','g','r'};
for ii=1:dimcount,
   [x,y]=circle(x0(ii),y0(ii),rad(ii),50);
   hold on
   if showlabels,
      text(solox(ii),soloy(ii),num2str(d(1+(ii==1),1+(ii==2),1+(ii==3))));
      if dimcount==3,
         text(pairx(ii),pairy(ii),num2str(d(2-(ii==1),2-(ii==2),2-(ii==3))));
      end
   end
   
   plot(x,y,pcol{ii},'LineWidth',2);
end
text(mean(x0),mean(y0),num2str(d(end,end,end)));
text(0.9*(x0(2)+rad(2)),0.9*(y0(1)-rad(1)),num2str(d(1,1,1)));

hold off

%axis equal
%axis off
