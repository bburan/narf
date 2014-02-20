function p4z4 ()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gains', 'delays', 'y_offset'}}, ...
                       'n_poles', 4, ...
                       'n_zeros', 4)));

end