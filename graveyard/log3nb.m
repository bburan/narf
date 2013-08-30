function log3nb()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
                                          
append_module(MODULES.nonlinearity.mdl(struct('phi', [-3 -log(0+10^-3)], ...
                                              'nlfn', @nl_log)));