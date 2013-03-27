function ret = fit_sboost()
global STACK;
n_params = length(pack_fittables(STACK));
n_jacks = 10;
n_iterations = max(30, n_params * 3);

function quickboost()
    fit_boost('score', n_iterations);
end

prior = zeros(size(pack_fittables(STACK)));
jackshrink(n_jacks, @quickboost, prior);
ret = NaN;
end