function syn12 ()

global MODULES STACK XXX META;

% FAST EXCITATORY DYNAMICS 
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

% SLOW INHIBITORY DYNAMICS
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

% Yet another 1st order filter
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
append_module(MODULES.weight_channels.mdl(struct('weights', [1;1;1;1;1;1], ...
                                                 'y_offset', 0, ...
                                                 'fit_fields', {{'y_offset', 'weights'}})));

% Add output nonlinearity                  
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:))-min(XXX{end}.dat.(ff).respavg(:));
curvature = 1 / resprange;
append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
           'phi', [0 meanresp*2 meanpred curvature curvature], ...
           'nlfn', @nl_sig_logistic)));
       
% Fit just the weights and output nonlinearity
fit12();

% Now go through and re-fit the filters
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