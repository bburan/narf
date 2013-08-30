function log1b()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-1 -log(0+10^-1)], ...
                                              'nlfn', @nl_log)));