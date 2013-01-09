function [x] = excise(x)
% Exrcises all rows of matrix x containing NaNs
x(any(isnan(x)'),:) = [];
