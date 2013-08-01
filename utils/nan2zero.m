% function y=nan2zero(x)
%
% convert NaN entries of x to zero
%
function y=nan2zero(x)
    
    y=x;
    y(isnan(y))=0;
    