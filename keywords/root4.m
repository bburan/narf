function root4()

global MODULES;
global STACK;

append_module(MODULES.nonlinearity.mdl(struct('phi', [4], ...
                                              'nlfn',@nl_root)));

STACK{end}{1}.auto_plot=[];
