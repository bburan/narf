function tune()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 1, ...
                                            'fit_fields', {{'coefs'}})));