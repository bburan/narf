function ap1z0b ()

global MODULES XXX STACK;

% Init the coefs to have the right dimensionality for 'stim'
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).stim);                  

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays'}}, ...
                       'n_poles', 1, ...
                       'n_zeros', 0, ...                       
                       'n_inputs', num_chans)));

% Initialize poles and delays to be more broad
mdl = STACK{end}{1};
STACK{end}{1}.poles = repmat(-20*[1:mdl.n_poles], num_chans, 1); 
STACK{end}{1}.delays = repmat(5, num_chans, 1); 

append_module(MODULES.normalize_channels);

end