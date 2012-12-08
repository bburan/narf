% function rt=threshsat(parms,x);
%
% apply threshold and saturation 
% threshold x<parms(1) and saturate x>parms(2)
%
function rt=threshsat(parms,x);

rt=x;
rt(rt>parms(2))=parms(2);
rt=(rt-parms(1)).*(rt>=parms(1));

