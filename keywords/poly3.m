function poly3()

global MODULES XXX;

ff=XXX{end}.training_set{1};
resp = XXX{end}.dat.(ff).respavg(:);
stim = XXX{end}.dat.(ff).stim(:);

phi_init = polyfit(stim(1:length(resp)), resp, 2);
append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                                  'phi', phi_init, ...
                                                  'nlfn', @polyval)));
