function ap1z0nd ()

global MODULES XXX STACK;

% Init the coefs to have the right dimensionality for 'stim'
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).stim);                  

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles'}}, ...
                       'n_poles', 1, ...
                       'n_zeros', 0, ...                       
                       'n_inputs', num_chans)));

STACK{end}{1}.delays = [15];
                   
append_module(MODULES.normalize_channels);


end