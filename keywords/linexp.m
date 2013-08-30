function cp1free()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [1 10 -2], ...
                                              'nlfn', @nl_linexp, ...
                                              'fit_fields', {{'phi'}})));
