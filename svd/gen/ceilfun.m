% function y=ceilfun(parms,x);
%
% parms=[m b]
% y=1./ (m * (1./x) + b)
%
function y=ceilfun(parms,x);

y=1./(parms(2)+parms(1)./x);

