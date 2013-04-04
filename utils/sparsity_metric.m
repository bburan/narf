function sparsity = sparsity_metric(coefs)
% sparsity = sparsity_metric(coefs)
%
% A metric of the sparsity of the coefficients of an FIR filter. 
% Mathematically, it is the ratio of the L1 norm to the L2 norm:
%         sparsity = L1 / L2
% 
% The metric therefore has a minimum value of one, which occurs when only 
% one coefficient is nonzero. The maximum value it can take is the total
% number of coefficients there are, which occurs when all coefficients are
% nonzero and of equal magnitude. 
%
% It is a comparable metric regardless of how many coefficients
% there are. If the sparsity number is 3.1, then it has a dimensionality or
% complexity of about 3 strong coefficients (and everything else is
% assumed to be nearly zero).
%
% The bounds can be demonstrated by the following:
%  for ii = 1:100,
%      X(ii) = ii;
%      nosignal = ones(1, ii);
%      Y1(ii) = sparsity_metric(nosignal);
%      perfect = zeros(1, ii);
%      perfect(1) = 1;
%      Y2(ii) = sparsity_metric(perfect);
%  end
% 
% plot(X, Y1, X, Y2);

ncoefs = numel(coefs);

if ncoefs < 2
    sparsity = 1;
    return
end

c_bar = sqrt(sum(coefs(:).^2)); % L2 norm
sparsity = sum(abs(coefs(:)) / c_bar)^2;
