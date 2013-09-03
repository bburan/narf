function dep1pifn()  % single, independent depression channel for
                     % each 

global MODULES;
global STACK XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', false)));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim1')));

fitSubstack();
savefilt=STACK{end}{1}.coefs;

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0.1], ...
                           'tau',      [40], ...
                           'per_channel', 1, ...
                           'input','stim1',...
                           'output','stim1',...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim2')));
STACK{end}{1}.coefs=savefilt./8;

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0.1], ...
                           'tau',      [40], ...
                           'per_channel', 1, ...
                           'input','stim2',...
                           'output','stim2',...
                           'fit_fields', {{'strength', 'tau'}})));

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
