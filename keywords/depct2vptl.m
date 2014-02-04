function depct2vptl()

global MODULES XXX STACK;
global ANESTHETIZED_MOUSE

append_module(MODULES.normalize_channels.mdl(struct('force_positive',true)));

if ANESTHETIZED_MOUSE==1,
    disp('dep2: Using longer depression time constant IC for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
    stop_exp=3.0;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
    stop_exp=2.0;
end

append_module(MODULES.depression_filter_bank.mdl(...
    struct('strength', [1.5 10], ...
           'tau',      [10 20], ...
           'tau_norm',dep_tau_norm,...
           'crosstalk',20,...
           'num_channels', 1, ...
           'fit_fields', {{'strength', 'tau', 'crosstalk'}})));

append_module(MODULES.add_nth_order_terms);

fns = fieldnames(XXX{end}.dat);
[T, S, C] = size(XXX{end}.dat.(fns{1}).stim);
if C==3,
    keep_channels=[1 2 3];
elseif C==10,
    keep_channels=[1:4 6 9];
elseif C==21,
    keep_channels=[1:6 9 14 18];
end

append_module(MODULES.subsample_channels.mdl(...
    struct('keep_channels',keep_channels)));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

fitSubstack([],10^-(stop_exp));
fitSubstack(length(STACK)-3,10^-(stop_exp+0.5));
