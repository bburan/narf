function dep1pifn()  % single, independent depression channel for
                     % each 

global MODULES STACK XXX;

fir();

% re-route output to 'stim1', which means excitatory synapse
STACK{end}{1}.output='stim1';
calc_xxx(2);

savefilt=STACK{end}{1}.coefs;

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [5], ...
                           'tau',      [20], ...
                           'tau_norm',100,...
                           'per_channel', 1, ...
                           'input','stim1',...
                           'fit_fields', {{'strength', 'tau'}})));
fitSubstack();
STACK{end}{1}.output='stim1';
calc_xxx(2);

%keyboard

append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                     'fit_fields', {{'coefs','baseline'}},...
                                     'output','stim2')));
STACK{end}{1}.coefs=savefilt./8;

append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [5], ...
                           'tau',      [20], ...
                           'tau_norm',100,...
                           'per_channel', 1, ...
                           'input','stim2',...
                           'output','stim2',...
                           'fit_fields', {{'strength', 'tau'}})));

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
fitSubstack();
