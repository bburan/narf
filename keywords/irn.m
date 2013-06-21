function irn()

global MODULES;

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct('fit_fields', {{'Vrest','V0','gL'}})));
