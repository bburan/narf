function v = channel_variance(x)
% Compute the variance of the signal across each row. 

v = var(x, 0, 2); % Variance across each row normalized by N-1
