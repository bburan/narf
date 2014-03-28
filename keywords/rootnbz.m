function rootnbz()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [2 0 0], ...
                                              'nlfn', @nl_root, ...
                                              'fit_fields', {{'phi'}})));
