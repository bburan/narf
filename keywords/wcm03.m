function wcm03()
% Weight channels, using a morlet wavelet
% Works on 'stim' by default. 
global MODULES;

append_module(MODULES.weight_channels.mdl(...
       struct('phifn', @wc_morlet, ...
              'phi', [2 7 0.0; 3 7 0.0; 4 7 0], ...
              'y_offset', [0; 0; 0], ...
              'fit_fields', {{'phi'}})));
