function nmse()

global MODULES META;

append_module(MODULES.mean_squared_error);
append_module(MODULES.correlation);

META.perf_metric = @pm_nmse;