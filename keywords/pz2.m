function pz2 ()

global MODULES;

append_module(MODULES.pole_zeros.mdl(...
                  struct('fit_fields', {{'poles', 'B', 'delay'}}, ...
                         'order', 2)));
              
end