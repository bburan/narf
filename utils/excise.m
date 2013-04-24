function x = excise(x)
% x = excise(x)
% Excises all rows of matrix x which contain a NaN at any position. 
x(any(isnan(x)'),:) = [];
