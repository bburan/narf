function booster()

global STACK;

phi_init = pack_fittables(STACK);
max_n_steps = max(50, length(phi_init(:)) * 5);
min_stepsize = 10^-99;
min_scoredelta = 10^-12;
relative_delta = false;
vary_stepsize = true;
starting_stepsize= 1;
fit_boost(max_n_steps, min_stepsize, min_scoredelta, relative_delta, vary_stepsize, starting_stepsize);