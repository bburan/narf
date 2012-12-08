function [outnl,xx,sigmoidparms]=findoutnl(predpsth,resp,showfit);

if ~exist('showfit'),
   showfit=1;
end

BINENTRIES=30;
if length(predpsth)./BINENTRIES<5,
   BINENTRIES=floor(length(predpsth)./5);
end


rcount=length(predpsth);
binstart=round(linspace(1,rcount+1,rcount/BINENTRIES));
BINCOUNT=length(binstart)-1;
binval=zeros(BINCOUNT,2);

if 0,
   BINENTRIES=20;
   BINCOUNT=10;
   
   rcount=length(predpsth);
   binstart=round(linspace(1,rcount+1,rcount/BINENTRIES));
   BINCOUNT=length(binstart)-1;
   %binstart=round(linspace(1,rcount+1,BINCOUNT+1));
   binval=zeros(BINCOUNT,2);
end

rmat=[predpsth resp];
rmat=sortrows(rmat);

for ii=1:BINCOUNT,
   binstop=binstart(ii+1)-1;
   binval(ii,:)=mean(rmat(binstart(ii):binstop,:),1);
end

outnl=binval(:,2);
xx=binval(:,1);

if sum(abs(binval(:,1)-mean(binval(:,1))))==0,
   sigmoidparms=zeros(4,1);
   disp('no dev from mean!');
   return
end

fitopt=optimset('Display','off');
%fitopt=zeros(14,1);
%fitopt(14)=2000;
%fitopt=foptions(fitopt);
x0=[mean(xx);(max(outnl)-min(outnl))./(max(xx)-min(xx));...
    max(outnl)-min(outnl);min(outnl)];

lb=[min(xx);  0;   0;   min(outnl)-abs(min(outnl))];
ub=[max(xx);  inf; inf; max(outnl)];

sigmoidparms=lsqcurvefit('sigmoid',x0,xx,outnl,[],[],fitopt);

if showfit,
   figure(1)
   clf
   plot(binval(:,1),binval(:,2));
   hold on
   plot(binval(:,1),sigmoid(sigmoidparms,binval(:,1)),'g');
   plot(binval(:,1),sigmoid(x0,binval(:,1)),'r');
   hold off
   
   legend('data','fit');
   drawnow
end
%keyboard



return

rmat=[predpsth resp];
rmat=sortrows(rmat);
size(rmat)
xx=rmat(:,1);
outnl=rmat(:,2);

x0=[mean(xx);max(outnl)./(max(xx)-min(xx));max(outnl)-min(outnl); min(outnl)];
sigmoidparms=lsqcurvefit('sigmoid',x0,predpsth,resp);


if 1,
   figure(1)
   clf
   plot(xx,outnl);
   hold on
   plot(xx,sigmoid(sigmoidparms,xx),'g');
   plot(xx,sigmoid(x0,xx),'r');
   hold off
   
   legend('data','fit');
   drawnow
end


