% function y=dcgain(beta,x,attcode)
% 
% y=beta(1) + beta(2) .* x;
%
function y=dcgain(beta,x,attcode)

if exist('attcode','var'),
   y=beta(attcode(:,1)) + beta(attcode(:,2)) .* x;
else
   y=beta(1) + beta(2) .* x;
end

