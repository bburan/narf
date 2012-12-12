function do_plot_time_series(stack, xxx, xseries, yseries)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');

    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    filt_idx = get(filt_pop, 'Value');
       
    dat = x.dat.(sf);
    
    plot(dat.(xseries), squeeze(dat.(yseries)(stim_idx,:,filt_idx)), ...
         pickcolor(filt_idx));
    axis tight;
    drawnow;
end