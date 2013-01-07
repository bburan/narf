function do_plot_channel_vs_time(stack, xxx, xfield, yfield)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, baphy_chan_idx] = get_baphy_plot_controls(stack);
    chan_idx = popup2num(mdl.plot_gui.selected_chan_popup);
    dat = x.dat.(sf);
    
    plot(dat.(xfield), dat.(yfield)(:, stim_idx, chan_idx), 'k-');
    axis tight;    
end