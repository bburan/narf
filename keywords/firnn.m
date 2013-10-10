function firnn()

global MODULES XXX;
global ANESTHETIZED_MOUSE

if ANESTHETIZED_MOUSE==1,
    disp('fir: Using longer num_coefs for mouse data.');
    dep_tau_norm=100;
    fir_num_coefs=20;
    stop_exp=3.0;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
    stop_exp=2.0;
end

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'coefs','baseline'}})));

fitSubstack([],10^-(stop_exp));
