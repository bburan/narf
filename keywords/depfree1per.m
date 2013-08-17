function depfree1per()

global MODULES;
global XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

sf=XXX{end}.training_set{1};
chan_count=size(XXX{end}.dat.(sf).stim,3);
init_str=repmat(0.1,[1 chan_count]);
init_tau=repmat(30,[1 chan_count]);

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', init_str, ...
                           'tau',      init_tau, ...
                           'per_channel', 1, ...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                            'fit_fields', {{'coefs','baseline'}})));

init10();
