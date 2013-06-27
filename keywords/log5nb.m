function log5nb()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
                                          
append_module(MODULES.nonlinearity.mdl(struct('phi', [-5 -log(0+10^-5)], ...
                                              'nlfn', @nl_log)));