function do_plot_channels_as_heatmap(stack, xxx, signal)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, baphy_chan_idx] = get_baphy_plot_controls(stack);
    
    imagesc(squeeze(x.dat.(sf).(signal)(:,stim_idx,:))');
    set(gca,'YDir','normal');
    axis tight;
    
end