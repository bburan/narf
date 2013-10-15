function dep1nn()

global MODULES XXX STACK;
global ANESTHETIZED_MOUSE

if ANESTHETIZED_MOUSE==1,
    disp('dep1nn: Using longer FIR for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
    stop_exp=3.0;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
    stop_exp=2.0;
end

meanstim = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'stim'));
str0=10 ./ meanstim;

append_module(MODULES.depression_filter_bank_nonorm.mdl(...
    struct('strength', str0, ...
           'tau',      [20], ...
           'gain',     1,...
           'tau_norm',dep_tau_norm,...
           'num_channels', 1, ...
           'fit_fields', {{'strength', 'tau', 'gain'}})));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                            'baseline',meanresp,...
                            'fit_fields', {{'coefs','baseline'}})));

fitSubstack([],10^-(stop_exp));
fitSubstack(length(STACK)-1,10^-(stop_exp+0.5));
