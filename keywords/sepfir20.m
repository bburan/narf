function sepfir20()

global MODULES XXX;

fir_num_coefs=20;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_separable_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'spec_weights','time_weights','baseline', 'v'}})));

fitSubstack([],10^-2);
