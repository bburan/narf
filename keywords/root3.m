function root3()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [3], 'nlfn', @nl_root)));