function expo()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [1 1 0], ...
                                                  'nlfn', @nl_exponential)));