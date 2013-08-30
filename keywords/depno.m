function depno()

global MODULES XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0,  0.2   0.2    0.2    0.2], ...
                           'tau',      [0,  25,   75,    150,   400], ...
                           'num_channels', 5)));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

%init10();
fitSubstack([],10^-2);
