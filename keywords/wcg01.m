function wcg01()
% Weight all input channels, producing 1 output channel.
% Uses a gaussian parameterization.
% Works on 'stim' by default. 
global MODULES STACK XXX;

append_module(MODULES.weight_channels_with_gaussians.mdl(...
       struct('mu', 1, ...
              'sigma', 1, ...
              'fit_fields', {{'mu', 'sigma'}})));
