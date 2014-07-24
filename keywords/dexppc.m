function dexppc()

global MODULES XXX STACK;

ff=XXX{end}.training_set{1};
meanresp=nanmean(XXX{end}.dat.(ff).respavg(:));
meanpred=nanmean(XXX{end}.dat.(ff).stim(:));
resprange = max(XXX{end}.dat.(ff).respavg(:))-min(XXX{end}.dat.(ff).respavg(:));
curvature = 1 / resprange;

mdl = STACK{end}{1};
x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, C] = size(x.dat.(sf).(mdl.output)); 

append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
           'phi', repmat([0 meanresp*2 meanpred curvature], C, 1), ...
           'nlfn', @nl_dexp)));

% If STIM has the right number of channels, we
% can fit it. If not, don't fit it. 
if size(XXX{end}.dat.(ff).stim, 3) == 1
    fitSubstack();
end
