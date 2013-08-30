function siglogIC()

global MODULES XXX STACK

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                              'phi', [0 meanresp*2 meanpred 0.5 0.5], ...
                              'nlfn', @nl_sig_logistic)));
fitSubStack();
