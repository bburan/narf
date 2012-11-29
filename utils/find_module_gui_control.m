function obj = find_module_gui_control(stack, gui_name)
% Return the PLOT gui object with name 'gui_name'. 
% Useful for finding a module's PLOT GUI control programmatically

for idx = 1:length(stack)
    if isfield(stack{idx}, 'plot_gui') && ...
       isfield(stack{idx}.plot_gui, gui_name)
        obj = stack{idx}.plot_gui.(gui_name);
        return
    end
end