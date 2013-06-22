function npfnl0()

global MODULES;

append_module(MODULES.nonparm_filter_nonlinearity.mdl(struct(...
    'leftzero', true)));
