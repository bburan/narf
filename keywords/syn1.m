function syn1()

global MODULES XXX STACK;

ff=XXX{end}.training_set{1};
resp = XXX{end}.dat.(ff).respavg(:);
stim = XXX{end}.dat.(ff).stim(:);
phi_init = polyfit(stim(1:length(resp)), resp, 2);
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', phi_init, ...
                                                  'nlfn', @polyval)));





append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0.1], ...
                           'tau',      [5], ...
                           'tau_norm',100,...
                           'num_channels', 1, ...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

fitSubstack([],10^-2);
fitSubstack(length(STACK)-2,10^-2.5);
