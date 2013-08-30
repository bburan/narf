function npnlx()

global MODULES;

append_module(MODULES.nonparm_nonlinearity.mdl(struct('splitter', @split_by_respfile,...
                                                      'unifier', @unify_respfiles)));
