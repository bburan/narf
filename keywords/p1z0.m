function p1z0 ()

global MODULES;

x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).stim);             

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays', 'gains', 'y_offset'}}, ...
                       'n_poles', 1, ...
                       'n_zeros', 0, ...
                       'n_inputs', num_chans)));


end