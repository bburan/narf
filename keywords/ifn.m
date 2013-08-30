function ifn()

global MODULES;
global STACK XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'baseline',meanresp,...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim')));
fitSubstack();
STACK{end}{1}.output='stim1';
calc_xxx(2);

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim2')));


% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
fitSubStack();
