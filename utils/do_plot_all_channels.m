function do_plot_all_channels(stack, xxx, xfield, yfield)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);
    
    [T, S, C] = size(x.dat.(sf).(yfield));
    
    hold on;
    for c = 1:C
        plot(dat.(xfield), dat.(yfield)(:, stim_idx, c), ...
             pickcolor(c));
    end
    axis tight;
    hold off;
end
