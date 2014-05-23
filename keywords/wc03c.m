function wc03c()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES STACK XXX;

signal = 'stim';
n_output_chans = 3;

% Compute number of input channels
fns = fieldnames(XXX{end}.dat);
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);

SF=1;

if n_input_chans>n_output_chans,
    w_init=SF*ones(n_input_chans, n_output_chans);
else
    w_init=SF*[eye(n_input_chans) ...
               ones(n_input_chans,n_output_chans-n_input_chans)];
end

append_module(MODULES.weight_channels.mdl(...
       struct('weights', w_init, ...
              'y_offset', zeros(n_output_chans, 1), ...
              'fit_fields', {{'weights'}})));