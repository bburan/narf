function siglogIC3()

global MODULES XXX STACK;

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
stimrange = max(XXX{end}.dat.(ff).stim(:)) - min(XXX{end}.dat.(ff).stim(:));
curvature = stimrange / 10;
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                              'phi', [0 meanresp*2 meanpred curvature curvature], ...
                              'nlfn', @nl_sig_logistic)));
