function linspecx()
% Weight all input channels, producing 1 output channel
% Works on 'stim' by default. 

% this keyword is initialized with order 0 - however, I will use it only in
% conjunction with "bestorder" fitting procedure, resulting in kernels
% whose order matches best the data

global MODULES;

append_module(MODULES.lindeberg_spectral.mdl(...
       struct('order', 0, ...
              'fit_fields', {{'bf', 's', 'add_factor', 'norm_factor'}})));
