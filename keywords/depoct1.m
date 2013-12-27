function depoct1()

global MODULES XXX STACK;
global ANESTHETIZED_MOUSE

append_module(MODULES.normalize_channels.mdl(struct('force_positive',true)));

if ANESTHETIZED_MOUSE==1,
    disp('dep1: Using longer depression time constant IC for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
    stop_exp=3.0;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
    stop_exp=2.0;
end

append_module(MODULES.depression_filter_bank.mdl(...
    struct('strength', [5], ...
           'tau',      [20], ...
           'tau_norm',dep_tau_norm,...
           'num_channels', 1, ...
           'fit_fields', {{'strength', 'tau'}})));
