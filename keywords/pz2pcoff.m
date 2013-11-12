function pz2pcoff ()

global MODULES;

append_module(MODULES.pole_zeros.mdl(...
                  struct('fit_fields', {{'poles', 'B', 'delay', 'y_offset'}}, ...
                         'delay_per_chan', true, ...
                         'order', 2)));
              
end