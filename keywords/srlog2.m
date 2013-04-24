function srlog2()

global MODULES;

append_module(MODULES.smooth_respavg);
append_module(MODULES.nonlinearity.mdl(struct('phi', [-2 -log(0+10^-2)], 'nlfn', @nl_log)));