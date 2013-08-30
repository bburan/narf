function boosto3()

global STACK;
phi_init = pack_fittables(STACK);
max_n_steps = max(50, floor(length(phi_init(:)) * 2.5));
fit_boostomatic(max_n_steps, 0, 0, sqrt(2), 4, 5);
