function ret = fit_slsq()
n_jacks = 10;
fitter = @fit_lsq;
prior = zeros(size(pack_fittables(STACK)));
jackshrink(n_jacks, fitter, prior);
ret = NaN;