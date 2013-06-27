function npnl2()

global MODULES;

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.nonparm_nonlinearity_2d.mdl(struct('bincount', 25,...
                                                  'smoothwindow',2)));
