function p1z1x2 ()

global MODULES XXX STACK;

% Init the coefs to have the right dimensionality for 'stim'
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).stim);                  

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gains', 'delays', 'y_offset'}}, ...
                       'n_poles', 1, ...
                       'n_zeros', 1, ...
                       'n_inputs', num_chans)));

fitlastpz({'stim'});

STACK{end}{1}.fit_fields = {};
STACK{end}{1}.output = 'stim1';
calc_xxx(length(STACK)); % Recalc last element

% Add a second p1z1
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'delays'}}, ...
                       'output', 'stim2', ...
                       'n_poles', 1, ...
                       'n_zeros', 1)));
                
fitlastpz({'stim1', 'stim2'});

% Refit all PZ models 
[~, mod_idxs] = find_modules(STACK, 'pole_zeros');
for jj = 1:length(mod_idxs)
    idx = mod_idxs{jj};
    STACK{idx}{1}.fit_fields = {'poles', 'zeros', 'delays'};
end

% The output is a normalized, weighted linear combination of all PZ modules
append_module(MODULES.concatenate_channels.mdl(struct('inputs', {{'stim1', 'stim2'}})));
append_module(MODULES.normalize_channels);
wc01(); fitSubstack(length(STACK),10^-5); % Fit just the weights

fit04a(); pop_module(); % Now fit everything

end