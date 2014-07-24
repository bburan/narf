function ap4z2 ()

global MODULES XXX STACK;

% Init the coefs to have the right dimensionality for 'stim'
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).stim);                  

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'delays'}}, ...
                       'n_poles', 4, ...
                       'n_zeros', 2, ...                       
                       'n_inputs', num_chans)));

append_module(MODULES.normalize_channels);


end