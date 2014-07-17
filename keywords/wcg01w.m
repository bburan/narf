function wcg01w()
% Weight channels, using a gaussian parameterization with fixed width
% Works on 'stim' by default. 
global MODULES STACK XXX;
wc_gaussian_fixed_width = @(p,z) wc_gaussian([p 0.3],z); 
append_module(MODULES.weight_channels.mdl(...
       struct('phifn', wc_gaussian_fixed_width, ...
              'phi', [2], ...
              'fit_fields', {{'phi'}})));
