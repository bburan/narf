function firpupil()

global MODULES STACK XXX;

dep_tau_norm=100;
fir_num_coefs=20;
stop_exp=3.0;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', fir_num_coefs, ...
                     'baseline',meanresp,...
                     'fit_fields', {{'coefs','baseline'}})));

fit05();
pop_module();  % trim the mean_squared_error module from the stack

npg();
fit05();
% nmse();
