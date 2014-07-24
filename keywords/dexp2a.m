function dexp2a()

global MODULES XXX;

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:))-min(XXX{end}.dat.(ff).respavg(:));
curvature = 1 / resprange;

nl_dexp3 = @(p,z) nl_dexp([0 20.0 p(1) p(2)], z);

append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
           'phi', [meanpred curvature], ...
           'nlfn', nl_dexp3)));

% If STIM has the right number of channels, we
% can fit it. If not, don't fit it. 
if size(XXX{end}.dat.(ff).stim, 3) == 1
    fitSubstack();
end

