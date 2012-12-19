function do_plot_output_vs_time(stack, xxx, xseries, yseries)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, baphy_chan_idx] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);
    
    plot(dat.(xseries), dat.(yseries)(:, stim_idx), 'k-');
    axis tight;    
end