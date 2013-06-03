function llexp3()

global MODULES META;

append_module(MODULES.bayesian_likelihood.mdl(struct('probdist', 'exponential', ...
                                                     'probcutoff', 0.003)));
append_module(MODULES.correlation);

META.perf_metric = @pm_nlogl;