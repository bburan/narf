function wcg02()
% Weight channels, using a gaussian parameterization
% Works on 'stim' by default. 
global MODULES STACK XXX;

append_module(MODULES.weight_channels.mdl(...
       struct('phifn', @wc_gaussian, ...
              'phi', [7 3; 3 3], ...
              'y_offset', [0; 0], ...
              'fit_fields', {{'phi'}})));
