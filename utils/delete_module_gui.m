function delete_module_gui(idx)
% delete_module_gui()
%
% Deletes any gui handles, widgets, etc from the module at NARFGUI index
% number idx. 
%
% No arguments or return values. Use purely for side effects. 

global NARFGUI;

try
    delete(NARFGUI{idx}.plot_axes);
    delete(NARFGUI{idx}.plot_popup);
    delete(NARFGUI{idx}.plot_panel);
    delete(NARFGUI{idx}.fn_apply);
    delete(NARFGUI{idx}.fn_table);
    delete(NARFGUI{idx}.fn_popup);
    delete(NARFGUI{idx}.fn_panel);
    delete(NARFGUI{idx}.fn_recalc);
    delete(NARFGUI{idx}.fn_replot);
catch
    % Do nothing if a delete failed
end

if isfield(NARFGUI{idx}, 'plot_gui')
    NARFGUI{idx} = rmfield(NARFGUI{idx}, 'plot_gui');
end
