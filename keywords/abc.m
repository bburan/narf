function abc ()

global MODULES;

append_module(MODULES.state_space_diffeq.mdl(...
                  struct('fit_fields', {{'A', 'B', 'C', 'D', 'delay_B', 'delay_B_amount'}})));

end