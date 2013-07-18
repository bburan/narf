function poly2()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [1 1 1], ...
                                                  'nlfn', @polyval)));