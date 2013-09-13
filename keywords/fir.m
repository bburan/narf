function fir()

global MODULES XXX;
global ANESTHETIZED_MOUSE

if ANESTHETIZED_MOUSE,
    disp('fir: Using longer num_coefs for mouse data.');
    dep_tau_norm=50;
    fir_num_coefs=20;
else
    dep_tau_norm=100;
    fir_num_coefs=12;
end

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'coefs','baseline'}})));

fitSubstack();
