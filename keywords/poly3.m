function poly3()

global MODULES XXX;

ff=XXX{end}.training_set{1};
meanresp = nanmean(XXX{end}.dat.(ff).respavg(:));
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', [meanresp 1 0 0], ...
                                                  'nlfn', @polyval)));
