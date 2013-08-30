function depfreewnn3()

global MODULES;

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0.1], ...
                           'tau',      [20], ...
                           'num_channels', 1, ...
                           'fit_fields', {{'strength', 'tau'}})));
                       

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                            'fit_fields', {{'coefs','baseline'}})));

init10();
