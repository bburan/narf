function wc03a()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES STACK XXX;
signal = 'stim';
n_output_chans = 3;

% Fit an FIR filter first and use its principal components
firtemp();
coefs = STACK{end}{1}.coefs';
B = zscore(abs(coefs));
V = princomp(B);
pop_module();

append_module(MODULES.weight_channels.mdl(...
       struct('weights', V(:, 1:n_output_chans), ...
              'y_offset', zeros(n_output_chans, 1), ...
              'fit_fields', {{'weights', 'y_offset'}})));

end