% function [beta,beta0]=fitgaussfp(fftkern,N,sfrange);
%
% beta=[or sf orw sfw a]  repeated for each N
%
function [beta,beta0]=fitgaussfp(fftkern,N,sfrange);

if ~exist('N','var'),
   N=1;
end
if ~exist('sfrange','var'),
   Xmax=size(fftkern,1);
   xc=round(Xmax/2);
   sfrange=(1:Xmax)-xc-1;
end
[xx,yy]=meshgrid(sfrange,sfrange);
x=[xx(:) yy(:)];

opt=optimset('Display','off');

% fit spatial frequency tuning
xx=-1:1;
g=exp(-(xx*1.5).^2);
g=g./sum(g(:));

tkern=fftkern;
beta0=[];
lb=[-inf 0.0001 pi/32 1.01 -inf];
ub=[inf max(sfrange).*sqrt(2) pi./2 max(sfrange)/2 inf];
%ub=[inf max(sfrange).*sqrt(2) max(sfrange) max(sfrange)/2 inf];

for nn=1:N,
   
   sfftkern=conv2(tkern,g,'same');
   sfftkern=conv2(sfftkern,g','same');
   %smaxidx=min(find(sfftkern==max(sfftkern(:))));
   smaxidx=min(find(abs(sfftkern)==max(abs(sfftkern(:)))));
   
   [w1,w2]=ind2sub(size(fftkern),smaxidx);
   w1=sfrange(w1);
   w2=sfrange(w2);
   
   if abs(w2)>0
      oo=-atan(w1./w2);
   else
      oo=pi/2;
   end
   
   ss=sqrt(w1.^2+w2.^2)+lb(2);
   
   tbeta0=[oo ss pi/20 1.25 tkern(smaxidx)];
   tbeta=lsqcurvefit('gaussfpN',tbeta0,x,tkern(:),lb,ub,opt);
   
   tkern=tkern-reshape(gaussfp(tbeta,x),Xmax,Xmax);
   beta0=[beta0 tbeta];
end
lb=repmat(lb,[1 N]);
ub=repmat(ub,[1 N]);


beta=lsqcurvefit('gaussfpN',beta0,x,fftkern(:),lb,ub,opt);

%fprintf('   sf0: or=%.1f sf=%.1f orw=%.3f sfw=%.3f a=%.3f\n',...
%        beta0(1),beta0(2),beta0(3),beta0(4),beta0(5));
%fprintf('   sf : or=%.1f sf=%.1f orw=%.3f sfw=%.3f a=%.3f\n',...
%        beta(1),beta(2),beta(3),beta0(4),beta(5));

beta=beta(:);
beta0=beta0(:);

return
