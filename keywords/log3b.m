function log3b()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-3 -log(0+10^-3)], ...
                                              'nlfn', @nl_log)));