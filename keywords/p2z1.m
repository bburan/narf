function p2z1 ()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gains', 'delays', 'y_offset'}}, ...
                       'n_poles', 2, ...
                       'n_zeros', 1)));

fitSubstack();

end