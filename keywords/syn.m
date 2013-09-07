function syn()

global MODULES XXX STACK;

% FIR FILTER
meanresp = nanmean(flatten_field(XXX{end}.dat,XXX{end}.training_set,'respavg'));
append_module(MODULES.fir_filter.mdl(struct('num_coefs', 12, ...
                                'baseline',meanresp,...
                                'fit_fields', {{'coefs','baseline'}})));
fitSubstack();

% SIGLOG
ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:)) - min(XXX{end}.dat.(ff).respavg(:));
curvature = 10 / resprange;
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                              'phi', [0 meanresp*2 meanpred curvature curvature], ...
                              'nlfn', @nl_sig_logistic)));

fitSubstack([],10^-2);

% DEPRESSION
append_module(MODULES.depression_filter_bank.mdl(...
                    struct('strength', [0.1], ...
                           'tau',      [5], ...
                           'tau_norm',100,...
                           'num_channels', 1, ...
                           'fit_fields', {{'strength', 'tau'}})));
append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
fitSubstack(length(STACK)-2,10^-2);

