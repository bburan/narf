function dexp3a()

global MODULES XXX;

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:))-min(XXX{end}.dat.(ff).respavg(:));
curvature = 1 / resprange;

nl_dexp3 = @(p,z) nl_dexp([3.0 p(1) p(2) p(3)], z);

append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
           'phi', [meanresp*2 meanpred curvature], ...
           'nlfn', nl_dexp3)));

% If STIM has the right number of channels, we
% can fit it. If not, don't fit it. 
if size(XXX{end}.dat.(ff).stim, 3) == 1
    fitSubstack();
end

