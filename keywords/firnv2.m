function firnv2()

global MODULES;

append_module(MODULES.normalize_channels);

append_module(MODULES.add_nth_order_terms);

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                            'fit_fields', {{'coefs'}})));

append_module(MODULES.normalize_channels);
