function logfreeb()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 0], ...
                                              'nlfn', @nl_log, ...
                                              'fit_fields', {{'phi'}})));
