function sel = narfgui_widget_selected_values(narfgui)
    sel = [];
    %TODO: This is just a stub!
    
    sel.stimfile = popup2str(narfgui{1}.plot_gui.selected_stimfile_popup);
    sel.stim_idx = popup2num(narfgui{1}.plot_gui.selected_stim_idx_popup);
    sel.chan_idx = popup2num(narfgui{1}.plot_gui.selected_stim_chan_popup);    
end