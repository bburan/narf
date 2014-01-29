function wc03()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES STACK XXX;

signal = 'stim';
n_output_chans = 3;

% Compute number of input channels
fns = fieldnames(XXX{end}.dat);
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);

append_module(MODULES.weight_channels.mdl(...
       struct('weights', zeros(n_input_chans, n_output_chans), ...
              'y_offset', zeros(n_output_chans, 1), ...
              'fit_fields', {{'weights', 'y_offset'}})));

fitSubstack([],10^-2);
STACK{end}{1}.fit_fields = {'weights'}; % No more y_offset