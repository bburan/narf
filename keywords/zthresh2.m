function zthresh2()

global MODULES;
disp('zthresh2...');

append_module(MODULES.nonlinearity.mdl(struct('fit_fields', {{'phi'}}, ...
                                              'phi', [0 100 0], ...
                                              'nlfn', @nl_zerothresh)));
fitSubstack()