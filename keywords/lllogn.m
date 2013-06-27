function lllogn()

global MODULES META;

append_module(MODULES.bayesian_likelihood.mdl(struct('probdist', 'Lognormal')));
append_module(MODULES.correlation);

META.perf_metric = @pm_nlogl;
