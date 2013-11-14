function pz2off ()

global MODULES;

append_module(MODULES.pole_zeros.mdl(...
                  struct('fit_fields', {{'poles', 'B', 'delay', 'y_offset'}}, ...
                         'order', 2)));
              
end