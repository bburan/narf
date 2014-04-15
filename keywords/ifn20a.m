function ifn20()

global MODULES STACK XXX;

% start with simple linear filter
fir20();
fitSubstack([],10^-3.5);

baseline=STACK{end}{1}.baseline;
STACK{end}{1}.sum_channels=0;
STACK{end}{1}.baseline=0;

signal = 'stim';

% Compute number of input channels
fns = fieldnames(XXX{end}.dat);
n_input_chans = size(STACK{end}{1}.coefs,1);
n_half=ceil(n_input_chans/2);
weights1=[ones(n_half,1); zeros(n_input_chans-n_half,1)];
weights2=[zeros(n_half,1); ones(n_input_chans-n_half,1)];

append_module(MODULES.weight_channels.mdl(...
       struct('weights', weights1, ...
              'input','stim_filtered','output', 'stim1', ...
              'y_offset', baseline, ...
              'fit_fields', {{'y_offset'}})));
append_module(MODULES.weight_channels.mdl(...
       struct('weights', weights2, ...
              'input','stim_filtered','output', 'stim2', ...
              'y_offset', baseline, ...
              'fit_fields', {{'y_offset'}})));

calc_xxx(2);

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.int_fire_neuron.mdl(struct(...
    'fit_fields', {{'Vrest','V0','gL'}},...
    'rectify_inputs',1)));
fitSubstack();
