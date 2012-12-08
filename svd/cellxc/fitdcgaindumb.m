% function [dcgparms,beta0]=fitdcgaindumb(predpsth,resp,attcode);
%
% dcgparms=[dc g]
% y= dc + (g * x)
%
% beta0 = initial guess at dcgparms
% attcode allows certain parameters to float with attention state.
%
% NOTE: assumes attcode ranges from 1...max for each column
%
% doesn't use lsqcurvefit.
%
function [dcgparms,beta0]=fitdcgaindumb(predpsth,resp,attcode);

if ~exist('attcode','var'),
   attcode=ones(length(predpsth),2);
   dcattcount=1;
   gattcount=1;
else
   % assume attcode ranges from 1...max for each column
   dcattcount=max(attcode(:,1));
   gattcount=max(attcode(:,2));
end

if length(resp)==0,
   dcgparms=[0;1];
   beta0=[0;1];
   return
end

% fit dc and gain terms to jackknife subsets of the data
resampcount=20;
dc=zeros(dcattcount,resampcount);
gn=zeros(gattcount,resampcount);

meanall=nanmean(resp);

for residx=1:resampcount,
   for dcidx=1:dcattcount,
      attidx=find(attcode(:,1)==dcidx);
      if resampcount>1,
         attidx=attidx([1:round((residx-1)/resampcount*length(attidx)) ...
                        round(residx/resampcount*length(attidx)+1):end]);
      end
      
      dc(dcidx,residx)=mean(resp(attidx));
   end
end

mdc=mean(dc,2);
edc=std(dc,1,2).*sqrt(resampcount-1);
dc=shrinkage(mdc,edc,1);
%[mdc-meanall edc].*60

tresp=resp;
for dcidx=1:dcattcount,
   attidx=find(attcode(:,1)==dcidx);
   tresp(attidx)=tresp(attidx)-dc(dcidx);
end

for residx=1:resampcount,
   for gnidx=1:gattcount,
      attidx=find(attcode(:,2)==gnidx);
      if resampcount>1,
         attidx=attidx([1:round((residx-1)/resampcount*length(attidx)) ...
                        round(residx/resampcount*length(attidx)+1):end]);
      end
      
      r0=tresp(attidx);
      r1=predpsth(attidx);
      d1=sum(r1.^2);
      
      if d1==0,
         % no modulation in predicted response
         dcgparms=[ones(dcattcount,1) .* mean(resp); ones(gattcount,1)];
         beta0=dcgparms;
      else
         gn(gnidx,residx)=sum(r0.*r1)./d1;
      end
   end
end

sgn=sign(sum(gn,2));
gn=gn.*repmat(sgn,[1 resampcount]);
gn(find(gn<=0))=1;

mgn=mean(log(gn),2);
egn=std(log(gn),1,2).*sqrt(resampcount-1);
gn=exp(shrinkage(mgn,egn,1)) .* sgn;


beta0=[dc;gn];
dcgparms=beta0;


