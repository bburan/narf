function log4b()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-4 -log(0+10^-4)], ...
                                              'nlfn', @nl_log)));