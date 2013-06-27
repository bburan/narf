function llexp()

global MODULES META;

append_module(MODULES.bayesian_likelihood.mdl(struct('probdist', 'exponential')));
append_module(MODULES.correlation);

META.perf_metric = @pm_nlogl;