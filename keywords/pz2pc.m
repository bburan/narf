function pz2pc ()

global MODULES;

append_module(MODULES.pole_zeros.mdl(...
                  struct('fit_fields', {{'poles', 'B', 'delay'}}, ...
                         'delay_per_chan', true, ...
                         'order', 2)));
              
end