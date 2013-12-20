function wc02p1z0 ()

global MODULES STACK XXX;
                   

append_module(MODULES.normalize_channels);

% Compute number of input channels
signal = 'stim';
n_output_chans = 1;
fns = fieldnames(XXX{end}.dat);
n_input_chans = size(XXX{end}.dat.(fns{1}).(signal), 3);

% First linear transformation
append_module(MODULES.weight_channels.mdl(...
       struct('weights', ones(n_input_chans, n_output_chans), ...
              'y_offset', 0.01, ...
              'output', 'stim1', ...
              'fit_fields', {{'weights'}})));          

% First p1z0
append_module(MODULES.normalize_channels.mdl(struct('input', 'stim1', 'output', 'stim1')));
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays'}}, ...
                       'input', 'stim1', ...
                       'n_poles', 1, ...
                       'n_zeros', 0)));

fitlastpz({'stim'});
STACK{end}{1}.fit_fields = {}; 
STACK{end}{1}.output = 'stim1'; 
calc_xxx(length(STACK)); % Recalc last element because output changed

% Add a second linear transform and p1z0
append_module(MODULES.weight_channels.mdl(...
       struct('weights', ones(n_input_chans, n_output_chans), ...
              'y_offset', 0.01, ...
              'output', 'stim2', ...
              'fit_fields', {{'weights'}})));                    
append_module(MODULES.normalize_channels.mdl(struct('input', 'stim2', 'output', 'stim2')));
append_module(MODULES.pole_zeros.mdl(...
                struct('fit_fields', {{'poles', 'delays'}}, ...
                       'input', 'stim2', ...
                       'output', 'stim2', ...
                       'n_poles', 1, ...
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