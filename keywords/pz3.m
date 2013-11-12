function pz3 ()

global MODULES;

append_module(MODULES.pole_zeros.mdl(...
                  struct('fit_fields', {{'poles', 'B', 'delay'}}, ...
                         'order', 3)));
              
end