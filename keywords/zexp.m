function zexp()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [1 10 -2 0], ...
                                              'nlfn', @nl_linexp, ...
                                              'fit_fields', {{'phi'}})));
