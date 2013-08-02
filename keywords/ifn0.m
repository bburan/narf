function ifn0()

global MODULES;
global STACK XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim1')));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim2')));

init10();

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
