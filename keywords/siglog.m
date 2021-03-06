function siglog()

global MODULES XXX;

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
%stimrange = max(XXX{end}.dat.(ff).stim(:)) - min(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:))-min(XXX{end}.dat.(ff).respavg(:));
%curvature = 1 / stimrange;
curvature = 1 / resprange;

append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
           'phi', [0 meanresp*2 meanpred curvature curvature], ...
           'nlfn', @nl_sig_logistic)));

% If STIM has the right number of channels, we
% can fit it. If not, don't fit it. 
if size(XXX{end}.dat.(ff).stim, 3) == 1
    fitSubstack();
end

