function firno16()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 16, ...
                                            'fit_fields', {{'coefs', 'baseline'}})));