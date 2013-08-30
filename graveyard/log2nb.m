function log2nb()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
                                          
append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], ...
                                              'nlfn', @nl_log)));