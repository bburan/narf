function log2()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2], ...
                                              'nlfn', @nl_log)));