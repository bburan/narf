function dep1()

global MODULES XXX STACK;
global ANESTHETIZED_MOUSE

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

if ANESTHETIZED_MOUSE,
    disp('dep1: Using longer depression time constant IC for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
end

append_module(MODULES.depression_filter_bank.mdl(...
    struct('strength', [0.1], ...
           'tau',      [5], ...
           'tau_norm',dep_tau_norm,...
           'num_channels', 1, ...
           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

fitSubstack([],10^-2);
fitSubstack(length(STACK)-2,10^-2.5);
