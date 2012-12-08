% function [beta,beta0]=fitgaussfp(fftkern,mS,bgpix,stimfilterparms);
%
% reasonable polar specs:
% beta(1) phase = [0 270];
% beta(2) fc = [0.005 0.01 0.02 0.03 0.04 0.06 0.08]
% beta(3) fr = [1 3 5 7]
% beta(4) pol = [-1 1]
% beta(5) A = [0 255]
%
function [beta,beta0]=fitpolfp(fftkern,mS,bgpix,stimfilterparms);

opt=optimset('Display','off');

% f=polfp(beta,r,mS,bgpix,stimfilterparms)

lb=[-inf 0.005 1 -1 -inf];
ub=[inf 0.08 7 1 inf];
beta0=[90 0.02 4 1 std(fftkern)];
r=sqrt(length(fftkern)*2)/2 .* ones(size(fftkern));

beta=lsqcurvefit('polfp',beta0,r,fftkern,lb,ub,opt,mS,bgpix,stimfilterparms);

return


% fit spatial frequency tuning
xx=-1:1;
g=exp(-(xx*1.5).^2);
g=g./sum(g(:));

tkern=fftkern;
beta0=[];
lb=[-inf 0.0001 pi/16 0.5 -inf];
ub=[inf max(sfrange) pi max(sfrange)/2 inf];

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
   
   tbeta0=[oo ss pi/6 1.5 tkern(smaxidx)];
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
