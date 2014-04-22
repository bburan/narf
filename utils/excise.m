function x = excise(x)
% x = excise(x)
% Excises all rows of matrix x which contain a NaN at any position. 
% For higher dimensional matrices, excises a hyper-row (dimensions 2 to N)
% if there is a NAN anywhere.

if ndims(x) < 3
    x(any(isnan(x)'),:) = [];
elseif ndims(x) == 3    
    [a,b,c] = size(x);
    y = reshape(x, a, b*c); 
    y(any(isnan(y)'),:) = [];   
    x = reshape(y, [], b, c);
else
    error('Not sure how to excise that!');
end

