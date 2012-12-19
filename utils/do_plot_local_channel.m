function do_plot_local_channel(stack, xxx, xseries, yseries)
    mdl = stack{end};
    x = xxx{end};

    % Find the GUI controls
    [sf, stim_idx, baphy_chan_idx] = get_baphy_plot_controls(stack);
    chan_idx = popup2num(mdl.plot_gui.selected_chan_popup);
    dat = x.dat.(sf);
    
    plot(dat.(xseries), dat.(yseries)(:, stim_idx, chan_idx), 'k-');
    axis tight;
    drawnow;
    
    do_plot_time_series(stack, xxx, mdl.input_time, mdl.output);
    
end