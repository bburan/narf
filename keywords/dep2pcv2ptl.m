function dep2pcv2ptl()
% single depression synapse, independent per stimulus channel, with
% normal intial conditions.

global MODULES XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

sf=XXX{end}.training_set{1};
chan_count=size(XXX{end}.dat.(sf).stim,3);
init_str=repmat([1 5],[1 chan_count]);
init_tau=repmat([10 20],[1 chan_count]);

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', init_str, ...
                           'tau',      init_tau, ...
                           'tau_norm',  100,...
                           'per_channel', 1, ...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.add_nth_order_terms);

fns = fieldnames(XXX{end}.dat);
[T, S, C] = size(XXX{end}.dat.(fns{1}).stim);
if C==3,
    keep_channels=[1 2 3];
elseif C==10,
    keep_channels=[1:4 6 9];
    %keep_channels=[6 9];
elseif C==21,
    keep_channels=[1:6 9 14 18];
end

append_module(MODULES.subsample_channels.mdl(...
    struct('keep_channels',keep_channels)));

