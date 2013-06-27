function ret = fit_sgene()
global STACK;
n_jacks = 10;
fitter = @fit_genetic;
prior = zeros(size(pack_fittables(STACK)));
jackshrink(n_jacks, fitter, prior);
ret = NaN;