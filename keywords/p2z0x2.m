function p2z0x2 ()

global MODULES STACK XXX;
                   
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays'}}, ...
                       'n_poles', 2, ...
                       'n_zeros', 0)));

fitlastpz({'stim'});

STACK{end}{1}.fit_fields = {};
STACK{end}{1}.output = 'stim1';
calc_xxx(length(STACK)); % Recalc last element

% Add a second p1z0
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays'}}, ...
                       'output', 'stim2', ...
                       'n_poles', 2, ...
                       'n_zeros', 0)));
                
fitlastpz({'stim1', 'stim2'});

% Refit all PZ models 
[~, mod_idxs] = find_modules(STACK, 'pole_zeros');
for jj = 1:length(mod_idxs)
    idx = mod_idxs{jj};
    STACK{idx}{1}.fit_fields = {'poles', 'delays'};
end

% The output is a normalized, weighted linear combination of all PZ modules
append_module(MODULES.concatenate_channels.mdl(struct('inputs', {{'stim1', 'stim2'}})));
append_module(MODULES.normalize_channels);
wc01(); fitSubstack(length(STACK),10^-5); % Fit just the weights

fit04a(); pop_module(); % Now fit everything

end