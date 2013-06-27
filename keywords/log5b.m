function log5b()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-5 -log(0+10^-5)], ...
                                              'nlfn', @nl_log)));