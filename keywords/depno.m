function depno()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0,  0.2   0.2    0.2    0.2], ...
                           'tau',      [0,  25,   75,    150,   400], ...
                           'num_channels', 5)));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                            'fit_fields', {{'coefs','baseline'}})));
