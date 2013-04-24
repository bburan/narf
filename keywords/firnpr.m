function firnpr()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('splitter', @split_by_respfile,...
                                            'unifier', @unify_respfiles, ...
                                            'num_coefs', 12, ...
                                            'fit_fields', {{'coefs'}})));

append_module(MODULES.normalize_channels);