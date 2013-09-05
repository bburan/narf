function dep1perchan()

global MODULES;
global XXX STACK;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

sf=XXX{end}.training_set{1};
chan_count=size(XXX{end}.dat.(sf).stim,3);
init_str=repmat(0.1,[1 chan_count]);
init_tau=repmat(3,[1 chan_count]);

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', init_str, ...
                           'tau',      init_tau, ...
                           'tau_norm',  100,...
                           'per_channel', 1, ...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

fitSubstack();
fitSubstack(length(STACK)-2,10^-2.5);
