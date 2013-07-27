function boosto2()

global STACK;
phi_init = pack_fittables(STACK);
max_n_steps = max(50, length(phi_init(:)) * 5);
fit_boostomatic(max_n_steps, 0, 0, sqrt(2), 2, 5);
