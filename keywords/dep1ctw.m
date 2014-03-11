% single, per-channel depression, with cross talk and weak ICs. no
% FIR filter (must be added separately);
function dep1ctw()

global MODULES;
global XXX STACK;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

sf=XXX{end}.training_set{1};
chan_count=size(XXX{end}.dat.(sf).stim,3);
init_str=repmat(1,[1 chan_count]);
init_tau=repmat(20,[1 chan_count]);

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', init_str, ...
                           'tau',      init_tau, ...
                           'tau_norm',  100,...
                           'crosstalk', 10,...
                           'per_channel', 1, ...
                           'fit_fields', {{'strength', 'tau','crosstalk'}})));
