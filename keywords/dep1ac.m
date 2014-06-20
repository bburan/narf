function dep1ac()
% single depression synapse, same synapse on all channels, with
% normal intial conditions.

global MODULES XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

sf=XXX{end}.training_set{1};
chan_count=size(XXX{end}.dat.(sf).stim,3);
init_str=5;
init_tau=20;
init_offset=0;

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', init_str, ...
                           'tau',      init_tau, ...
                           'offset_in',init_offset, ...
                           'tau_norm',  100,...
                           'per_channel', 0, ...
                           'fit_fields', {{'strength', 'tau', 'offset_in'}})));
