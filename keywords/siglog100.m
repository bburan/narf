function siglog100()

global MODULES XXX;

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
stimrange = max(XXX{end}.dat.(ff).stim(:)) - min(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:))-min(XXX{end}.dat.(ff).respavg(:));
%curvature = 1 / stimrange;
curvature = 1 / sqrt(stimrange*resprange).*100;

append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
            'fit_constraints', {{struct( ...
                                     'var', 'phi', ...
                                     'lower', [-100 0 -100 -100 -100], ...
                                     'upper', [100 200 200 100 100])}}, ...
           'phi', [0 meanresp*2 meanpred curvature curvature], ...
           'nlfn', @nl_sig_logistic100)));

fitSubstack();

