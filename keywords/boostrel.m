function boostrel()
global STACK;
phi_init = pack_fittables(STACK);
fit_boost(length(phi_init(:)) * 5, 10^-9, 10^-3, true);

