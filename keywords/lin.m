function lin()
% Linear weighting of all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES STACK XXX;

signal = 'stim';
n_output_chans = 1;

% Compute number of input channels
fns = fieldnames(XXX{end}.dat);
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);

append_module(MODULES.weight_channels.mdl(...
       struct('weights', ones(n_input_chans, n_output_chans), ...
              'y_offset', 0.0, ...
              'fit_fields', {{'weights'}})));
