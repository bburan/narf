function pz4 ()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles_real', 'poles_img', 'zeros_real', 'gains', 'delays', 'y_offset'}}, ...
                       'order', 4)));
              
end