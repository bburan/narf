function depfree2w()

global MODULES XXX STACK;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0.01 0.2], ...
                           'tau',      [10 200], ...
                           'num_channels', 2, ...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

%init10();
fitSubstack(size(STACK,1),10^-1);
fitSubstack(size(STACK,1)-2,10^-2);
