function logfreep()

global MODULES STACK;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(10^-2)], ...
                                              'nlfn', @nl_log, ...
                                              'fit_fields', {{'phi'}})));
STACK{end}{1}.auto_plot=[];
