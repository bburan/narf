function err20()

global MODULES META;

append_module(MODULES.error_norm.mdl(struct('pnorm', 2.0, ...
                                            'output', 'score')));
append_module(MODULES.correlation);

META.perf_metric = @pm_norm;