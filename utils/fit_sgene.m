function ret = fit_sgene()
global STACK;
n_jacks = 10;
fitter = @fit_genetic;
jackshrink(n_jacks, fitter, prior);
ret = NaN;