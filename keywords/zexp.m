function zexp()

global MODULES XXX;

ff=XXX{end}.training_set{1};

append_module(MODULES.nonlinearity.mdl(struct('phi', [1 10 -2 0], ...
                                              'nlfn', @nl_linexp, ...
                                              'fit_fields', {{'phi'}})));
if size(XXX{end}.dat.(ff).stim, 3) == 1
    fitSubstack();
end