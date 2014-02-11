function voltp()

global MODULES XXX;

append_module(MODULES.normalize_channels.mdl(struct('force_positive', true)));

append_module(MODULES.add_nth_order_terms);
