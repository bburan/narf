function xs = jackknifed_prediction(x, n_jackknifes)

% Compute how much to shrink each parameter
m = length(phi_init);
n = n_jacks;
mu = mean(phi_jack, 2);
sigma_sq = var(phi_jack, [], 2);

% James-Stein Estimator that shrinks towards the prior.
if m < 3
    fprintf('WARNING: A James Stein Estimator is only better if there are 3 or more params. Using simple mean.');
    phi_best = mu;
else
    phi_best = prior + (mu - prior) .* (  1 - (((m-2) * (sigma_sq / n)) / sum((mu - prior).^2))); 
end

% Stephen's shrinkage eqn is essentially the same, but without prior.
% phi_best = real(mu .* sqrt(1 - ((sqrt(sigma_sq) / sqrt(n)) ./ mu).^2));

% Recalc all the way through
unpack_fittables(phi_best);
XXX{2} = cache;
calc_xxx(2);
termcond = NaN;

end
