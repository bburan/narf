% function [sigparms,beta0]=fitsigmoid2(predpsth,resp,showfit,attcode);
%
% dcgparms=[dc g]
% y= dc + (g * x)
%
% beta0 = initial guess at dcgparms
% attcode allows certain parameters to float with attention state.
%
% NOTE: assumes attcode ranges from 1...max for each column
%
function [sigparms,beta0]=fitsigmoid2(predpsth,resp,showfit,attcode);

if ~exist('attcode','var'),
   attcode=ones(length(predpsth),2);
   attcounts=[1 1 1 1];
else
   % assume attcode ranges from 1...max for each column
   attcounts=max(attcode,[],1);
end

% shift reference to align with position in beta vector
attcode(:,2)=attcode(:,2)+attcounts;

r0=resp-mean(resp);
r1=predpsth-mean(predpsth);
d1=sum(r1.^2);

if d1==0,
   % no modulation in predicted response
   sigparms=[ones(dcattcount,1) .* mean(resp); zeros(gattcount,1)];
   beta0=sigparms;
   
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

sigparms=lsqcurvefit('dcgain',beta0,predpsth,resp,lb,ub,fitopt,attcode);


%[beta0 sigparms]