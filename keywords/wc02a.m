function wc02a()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES STACK XXX;
signal = 'stim';
n_output_chans = 2;

% Fit an FIR filter first and use its principal components
firtemp();
coefs = STACK{end}{1}.coefs';
B = zscore(abs(coefs));
V = princomp(B);
pop_module();

append_module(MODULES.weight_channels.mdl(...
       struct('weights', V(:, 1:n_output_chans), ...
              'fit_fields', {{'weights'}})));

end