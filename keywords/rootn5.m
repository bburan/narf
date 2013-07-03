function rootn5()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.nonlinearity.mdl(struct('phi', [5], 'nlfn', @nl_root)));