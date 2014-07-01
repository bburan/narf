function wc02n()
% function wc02n()
%
% Weight all input channels, producing 2 output channels
% Works on 'stim' by default. 
% n:  initialize all weights to 0.1 
%
global MODULES STACK XXX;

% normalize first - doesn't help!
%append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

signal = 'stim';
n_output_chans = 2;

% Compute number of input channels
fns = fieldnames(XXX{end}.dat);
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);

append_module(MODULES.weight_channels.mdl(...
       struct('weights', 0.1*ones(n_input_chans, n_output_chans), ...
              'y_offset', repmat(0.01, n_output_chans, 1), ...
              'fit_fields', {{'weights', 'y_offset'}})));

