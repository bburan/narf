function mse()

global MODULES META;

append_module(MODULES.mean_squared_error.mdl(struct('output', 'score')));
append_module(MODULES.correlation);

META.perf_metric = @pm_mse;