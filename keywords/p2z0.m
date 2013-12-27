function p2z0 ()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays'}}, ...
                       'n_poles', 2, ...
                       'n_zeros', 0)));

% The output is a normalized, weighted linear combination of the above. 
append_module(MODULES.concatenate_channels.mdl(struct('inputs', {{'stim'}})));
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
wc01();

fit04a(); pop_module();

end