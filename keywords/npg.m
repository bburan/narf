function npg()

global MODULES;

% default inputs are already stim1, stim2
% default output already stim
append_module(MODULES.nonparm_gain.mdl(struct('bincount', 5,...
   'smoothwindow',0)));
