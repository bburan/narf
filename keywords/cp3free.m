function cp3free()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [1 1 1], ...
                                              'nlfn', @nl_invpoly, ...
                                              'fit_fields', {{'phi'}})));
