% function y=hinge4(beta,x,attcode)
% 
%   y=beta(1) + beta(2) .* x;       % dcg
%   y=(y-beta(3)) .* (x > beta(3)); % threshold
%
% created SVD 1/04
%
function y=hinge4(beta,x,attcode)

if exist('attcode','var'),
   y=beta(attcode(:,1)) + beta(attcode(:,2)) .* x;
   thr=beta(attcode(:,3));
   y=(y-thr) .* (y > thr);
else
   y=beta(1) + beta(2) .* x;
   y=(y-beta(3)) .* (x > beta(3));
end

