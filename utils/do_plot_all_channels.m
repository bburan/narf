function do_plot_all_channels(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);
    
    [T, S, C] = size(x.dat.(sf).(mdl.output));
    
    hold on;
    for c = 1:C
        plot(dat.(mdl.time), dat.(mdl.output)(:, stim_idx, c)), ...
             pickcolor(c));
    end
    axis tight;
    hold off;
end
