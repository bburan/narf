function sigumber()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [0 1 1 1 1], ...
                                                  'nlfn', @nl_sig_gumber)));
