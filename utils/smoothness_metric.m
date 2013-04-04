function smoothness = smoothness_metric(coefs)
% smoothness = smoothness_metric(coefs)
% 
% ARGUMENTS: 
%  coefs      Assumed to be FIR coefficients from the fir_filter module.
%             Unlike most data structures time is contained in the 2nd
%             dimension and not the first.
%
% RETURNS: 
%  smoothness The root of the sum of squares of the first derivative of
%             coefs along dimension 2 (assumed to be the time dimension).

smoothness = filter([1,-1], 1, coefs, [], 2);
smoothness = sqrt(sum(smoothness(:).^2));
