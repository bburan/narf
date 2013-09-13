function root2()

global MODULES;
global STACK;

append_module(MODULES.nonlinearity.mdl(struct('phi', [2], ...
                                              'nlfn',@nl_root)));

STACK{end}{1}.auto_plot=[];
