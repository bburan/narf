% function s=expthresh(parms,x);
%
% parms=[thresh expon scf];
%

function s=expthresh(parms,x);

if sum(parms)==0,
   s=x;
   return
end

thresh=parms(1);
expon=parms(2);
scf=parms(3);

s=((x-thresh).*(x>thresh)).^expon .* scf;


