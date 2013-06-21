function logfit()

global MODULES STACK;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], ...
                                              'nlfn', @nl_log, ...
                                              'fit_fields', {{'phi'}})));
append_module(MODULES.correlation);
META.perf_metric = @pm_corr;
fit_fminlsq();
pop_module();
STACK{end}{1} = rmfield(STACK{end}{1}, 'fit_fields');