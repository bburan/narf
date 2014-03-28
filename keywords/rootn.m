function rootn()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [2], ...
                                              'nlfn', @nl_root, ...
                                              'fit_fields', {{'phi'}})));
