function log2n()

global MODULES STACK;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(10^-2)], ...
                                              'nlfn', @nl_log)));
STACK{end}{1}.auto_plot=[];
