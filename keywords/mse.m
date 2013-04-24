function mse()

global MODULES;

append_module(MODULES.mean_squared_error.mdl(struct('output', 'score')));
append_module(MODULES.correlation);