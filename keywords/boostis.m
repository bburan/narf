function boostis()
global STACK;
phi_init = pack_fittables(STACK);
fit_iteratively(@(~) fit_boost(max(5, length(phi_init)/5)), 100);

