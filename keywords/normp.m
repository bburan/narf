function normp()

global MODULES;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));
