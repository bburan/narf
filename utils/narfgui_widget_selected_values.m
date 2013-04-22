function sel = narfgui_widget_selected_values(narfgui)
% This function is like calc_xxx(), but for plotting. 
% It goes through the entire narfgui structure passed to it and returnse
% the selected stimfile, stim_idx, and channel at the end.
%
% TODO: If desired, this could be sped up by caching the values instead of
% recomputing them from scratch every time...although when debugging it
% might be easy to get into a strang state so this is a safer default.

sel = [];

for ii = 1:length(narfgui)
    gui = narfgui{ii};

    if isfield(gui, 'plot_gui')
        if isfield(gui.plot_gui, 'selected_stimfile_popup')
            sel.stimfile = popup2str(narfgui{1}.plot_gui.selected_stimfile_popup);
        end
        if isfield(gui.plot_gui, 'selected_stim_idx_popup')
            sel.stim_idx = popup2num(narfgui{1}.plot_gui.selected_stim_idx_popup);
        end
        if isfield(gui.plot_gui, 'selected_stim_chan_popup')
            sel.chan_idx = popup2num(narfgui{1}.plot_gui.selected_stim_chan_popup);    
        end
    end
end
