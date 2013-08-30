function cp2free()

global MODULES;
    
append_module(MODULES.nonlinearity.mdl(struct('phi', [0.5 0.1 1 0], ...
                                              'nlfn', @nl_sigmoid, ...
                                              'fit_fields', {{'phi'}})));