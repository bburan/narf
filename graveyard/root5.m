function root5()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [5], 'nlfn', @nl_root)));