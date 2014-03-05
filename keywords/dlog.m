function dlog()

global MODULES STACK;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2], ...
                                              'nlfn', @nl_dlog, ...
                                              'fit_fields', {{'phi'}})));
STACK{end}{1}.auto_plot=[];
