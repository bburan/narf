function [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack)

    % I end up using this code so much, I just made it into a little macro
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    sf = popup2str(baphy_mod.plot_gui.selected_stimfile_popup);
    stim_idx = popup2num(baphy_mod.plot_gui.selected_stim_idx_popup);
    chan_idx = popup2num(baphy_mod.plot_gui.selected_stim_chan_popup);
    
end