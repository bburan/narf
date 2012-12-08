% function y=hinge(parms,x,oparms);
%
% parms=[x0 y0 slope];
%

function y=hinge(parms,x,oparms);

if exist('oparms','var'),
   parms(find(~isnan(oparms)))=oparms(find(~isnan(oparms)));
end

x0=parms(1);
y0=parms(2);
slope=parms(3);

y=y0+ (x-x0) .* slope .* (x>x0);
