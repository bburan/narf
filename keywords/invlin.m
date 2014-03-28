function invlin()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [1 1], ...
                                              'nlfn', @nl_invlin, ...
                                              'fit_fields', {{'phi'}})));
