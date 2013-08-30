function fir()

global MODULES XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'coefs','baseline'}})));

fitSubstack();
