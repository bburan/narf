function [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack)
    global XXX; % TODO: Make get_baphy_plot_Controls take XXX as an argument

    [baphy_mod, ~] = find_modules(STACK, 'load_stim_resps_from_baphy', true);
    
    % If we can't find the plot_gui, just return a default in case we want
    % to use one of STACK's plot function from a script (ie, when there is
    % no narf_modelpane gui connected and things are headless)
    if isfield(baphy_mod, 'plot_gui')
        sf = popup2str(baphy_mod.plot_gui.selected_stimfile_popup);
        stim_idx = popup2num(baphy_mod.plot_gui.selected_stim_idx_popup);
        chan_idx = popup2num(baphy_mod.plot_gui.selected_stim_chan_popup);
    else
        sf = XXX{end}.training_set{1};
        stim_idx = 1;
        chan_idx = 1;
    end
    
end