function npfnlx()

global MODULES;

append_module(MODULES.nonparm_filter_nonlinearity.mdl(struct('splitter', @split_by_respfile,...
                                                  'unifier', @unify_respfiles,...
                                                  'gwinval',2)));
