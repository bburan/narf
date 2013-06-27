function bic()

global MODULES META;

append_module(MODULES.bayesian_likelihood);

META.perf_metric = @pm_bic;