% function [dcgparms,beta0]=fitdcgain(predpsth,resp,attcode,resampcount);
%
% dcgparms=[dc g]
% y= dc + (g * x)
%
% beta0 = initial guess at dcgparms
% attcode allows certain parameters to float with attention state.
%
% NOTE: assumes attcode ranges from 1...max for each column
%       resampcount>1 (default=1) jackknife to find gain. doesn't
%       seem to help?
%
function [dcgparms,beta0]=fitdcgain(predpsth,resp,attcode,resampcount);

if ~exist('attcode','var'),
   attcode=ones(length(predpsth),2);
   dcattcount=1;
   gattcount=1;
   
else
   % assume attcode ranges from 1...max for each column
   dcattcount=max(attcode(:,1));
   gattcount=max(attcode(:,2));
end

if ~exist('resampcount','var'),
   resampcount=1; %min([length(resp),20]);
end

if length(resp)==0,
   dcgparms=[0;0];
   beta0=[0;0];
   return
end

% shift reference to align with position in beta vector
attcode(:,2)=attcode(:,2)+dcattcount;

r0=resp-mean(resp);
r1=predpsth-mean(predpsth);
d1=sum(r1.^2);

if d1==0,
   % no modulation in predicted response
   dcgparms=[ones(dcattcount,1) .* mean(resp); zeros(gattcount,1)];
   beta0=dcgparms;
   
   return
end

scf=sum(r0.*r1)./d1;

beta0=[ones(dcattcount,1) .* mean(resp); ...
       ones(gattcount,1)  .* scf];
lb=[ones(dcattcount,1) .* -inf; ...
    ones(gattcount,1)  .* 0.0];
ub=[ones(dcattcount,1) .* inf; ...
    ones(gattcount,1)  .* inf];

fitopt=optimset('Display','off');

if resampcount<=1,
   dcgparms=lsqcurvefit('dcgain',beta0,predpsth,resp,lb,ub,fitopt,attcode);
else
   dcgp=zeros(length(beta0),resampcount);
   for residx=1:resampcount,
      attidx=[1:round((residx-1)/resampcount*length(resp)) ...
              round(residx/resampcount*length(resp)+1):length(resp)]';
      
      r0=resp(attidx);
      r1=predpsth(attidx);
      d1=sum(r1.^2);
      
      if d1==0,
         % no modulation in predicted response
         dcgp(:,residx)=[ones(dcattcount,1) .* mean(resp); 
                         ones(gattcount,1)];
      elseif 1,
         dcgp(:,residx)=lsqcurvefit('dcgain',beta0,r1,r0-r1,lb,ub,...
                                    fitopt,attcode(attidx,:));
         % add 1 to gains
         dcgp(dccount+1:end,residx)=dcgp(dccount+1:end,residx)+1;
      else
         dcgp(:,residx)=lsqcurvefit('dcgain',beta0,r1,r0,lb,ub,...
                                    fitopt,attcode(attidx,:));
      end
   end
   
   %mm=mean(dcgp(1:dcattcount,:),2)-mean(resp);
   %sm=std(dcgp(1:dcattcount,:),1,2).*sqrt(resampcount-1);
   %m=shrinkage(mm,sm,1)+mean(resp);
   
   % don't use shrinkage on mean since SNR seems to be reasonably high
   m=mean(dcgp(1:dcattcount,:),2);
   
   mg=mean(log(dcgp(dcattcount+(1:gattcount),:)),2);
   sg=std(log(dcgp(dcattcount+(1:gattcount),:)),1,2).*sqrt(resampcount-1);
   s=exp(shrinkage(mg,sg,1));
   s(find(isnan(s)))=1;
   
   dcgparms=[m;s];
end

%[beta0 dcgparms]

