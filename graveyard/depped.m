function depped()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0, 0.15 0.35 0.15 0.35], ...
                           'tau',      [0, 200,    200,   400,     400], ...
                           'pedestal', 0.1,...
                           'num_channels', 5)));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 10, ...
                                               'fit_fields', {{'coefs'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive',true)));