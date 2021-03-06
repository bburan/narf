function sepfir15()

global MODULES XXX;

fir_num_coefs=15;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_separable_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'spec_weights','time_weights','baseline', 'v'}}, ...
              'fit_constraints', {{struct( ...
                                         'var', 'v', ...
                                         'lower', -5, ...
                                         'upper', 5), ...
                                   struct( ...
                                         'var', 'spec_weights', ...
                                         'lower', -50, ...
                                         'upper', 50), ...
                                   struct( ...
                                         'var', 'time_weights', ...
                                         'lower', -50, ...
                                         'upper', 50), ...
                                   struct( ...
                                         'var', 'baseline', ...
                                         'lower', -50, ...
                                         'upper', 50)}} ...
              )));

fitSubstack([],10^-2);
