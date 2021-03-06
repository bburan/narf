function wcg04()
% Weight channels, using a gaussian parameterization
% Works on 'stim' by default. 
global MODULES STACK XXX;

append_module(MODULES.weight_channels.mdl(...
       struct('phifn', @wc_gaussian, ...
              'phi', [2 5; 3 5; 4 5; 5 5], ...
              'y_offset', [0; 0; 0; 0], ...
              'fit_fields', {{'phi'}})));
