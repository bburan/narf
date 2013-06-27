function err10()

global MODULES META;

append_module(MODULES.error_norm.mdl(struct('pnorm', 1.0, ...
                                            'output', 'score')));
append_module(MODULES.correlation);

META.perf_metric = @pm_norm;