function linspec1()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 
global MODULES;

append_module(MODULES.lindeberg_spectral.mdl(...
       struct('order', 1, ...
              'fit_fields', {{'bf', 's', 'add_factor', 'norm_factor'}})));
