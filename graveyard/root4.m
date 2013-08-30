function root4()

global MODULES;

append_module(MODULES.nonlinearity.mdl(struct('phi', [4], 'nlfn', @nl_root)));