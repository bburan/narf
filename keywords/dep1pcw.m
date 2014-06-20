function dep1pcw()
% single depression synapse, independent per stimulus channel, with
% WEAK intial conditions.

global MODULES XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

sf=XXX{end}.training_set{1};
chan_count=size(XXX{end}.dat.(sf).stim,3);
init_str=repmat(1,[1 chan_count]);
init_tau=repmat(20,[1 chan_count]);
init_offset=zeros(1,chan_count);

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', init_str, ...
                           'tau',      init_tau, ...
                           'tau_norm',  100,...
                           'offset_in',init_offset,...
                           'per_channel', 1, ...
                           'fit_fields', {{'strength', 'tau', 'offset_in'}})));
