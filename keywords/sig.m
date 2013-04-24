function sig()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [0 1 1 0], ...
                                                  'nlfn', @nl_sigmoid)));
