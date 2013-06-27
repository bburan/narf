function qboost()
global STACK;
phi_init = pack_fittables(STACK);
fit_boost(length(phi_init));