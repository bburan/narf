function ret = fit_sboost()
n_params = pack_fittables(STACK);
n_jacks = 10;
n_iterations = n_params * 10;
fitter = @(x) fit_boost('score', n_iterations);
prior = zeros(size(pack_fittables(STACK)));
jackshrink(n_jacks, fitter, prior);
ret = NaN;