function syn11 ()

global MODULES STACK XXX META;

% FAST DYNAMICS 
append_module(MODULES.pz_synapse.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gain', 'delayms', ...
                                       'prephi', 'postphi', 'y_offset'}}, ...
                       'poles', [-100 -60 -55 -5], ...
                       'zeros', [-470 -2 -1], ...
                       'gain', 140, ...
                       'prefn', @nl_log, ...
                       'prephi', [-2], ...
                       'postfn', @nl_softzero, ...
                       'postphi', [0.4 1.5 -2 2.6], ...
                       'input', 'stim', ...
                       'output', 'stim1')));                                      
                 
append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1'}})));
fit04a(); pop_module(); pop_module();
STACK{end}{1} = rmfield(STACK{end}{1}, 'fit_fields');

% SLOW DYNAMICS
append_module(MODULES.pz_synapse.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gain', 'delayms', ...
                                       'y_offset'}}, ...
                       'poles', [-2], ...
                       'zeros', [], ...
                       'gain', -10, ...
                       'prefn', @polyval, ... % Polyval [1 0] is identity
                       'prephi', [1 0], ...
                       'postfn', @polyval, ...
                       'postphi', [1 0], ...
                       'input', 'stim', ...
                       'output', 'stim2')));     

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2'}})));
fit04a(); pop_module(); pop_module();
STACK{end}{1} = rmfield(STACK{end}{1}, 'fit_fields');

% Another Dynamics Thing
append_module(MODULES.pz_synapse.mdl(...
                struct('fit_fields', {{'poles', 'zeros', 'gain', 'delayms', ...
                                       'prephi', 'postphi', 'y_offset'}}, ...
                       'poles', [-50], ...
                       'zeros', [], ...
                       'gain', 10, ...
                       'prefn', @nl_log, ...
                       'prephi', [-2], ...
                       'postfn', @nl_softzero, ...
                       'postphi', [-1 1 0 0], ...
                       'input', 'stim', ...
                       'output', 'stim3'))); 

append_module(MODULES.sum_fields.mdl(struct('inputs', {{'stim1', 'stim2', 'stim3'}})));
fit04a(); pop_module(); pop_module();
STACK{end}{1} = rmfield(STACK{end}{1}, 'fit_fields');

append_module(MODULES.concatenate_channels.mdl(struct('inputs', {{'stim1', 'stim2', 'stim3'}})));
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
append_module(MODULES.add_nth_order_terms);

% Normalize and combine 1st and 2nd order terms
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 1, ...
                                            'baseline', meanresp,...
                                            'fit_fields', {{'coefs','baseline'}})));

% Initialize weighting coefficients to be all 1                                        
STACK{end}{1}.coefs = ones(size(STACK{end}{1}.coefs));
fitSubstack([],10^-2);
siglog();
fit04a(); 

% Now go through and re-fit everything.
[~, mod_idxs] = find_modules(STACK, 'pz_synapse');

for jj = 1:length(mod_idxs)
    idx = mod_idxs{jj};
    STACK{idx}{1}.fit_fields = {'poles', 'zeros', 'gain', 'delayms', ...
                                'prephi', 'postphi'};
end

calc_xxx(2);
[a,b,c,d] = fit_scaat('InitStepSize', 100.0, ...
                      'StopAtAbsScoreDelta', 10^-4);


end