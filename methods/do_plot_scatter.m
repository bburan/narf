function do_plot_scatter(stack, xxx, field1, field2)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    plot(dat.(field1)(:, stim_idx), ...
         dat.(field2)(:, stim_idx), 'k.');
    axis tight;
end