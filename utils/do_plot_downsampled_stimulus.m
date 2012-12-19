function do_plot_downsampled_stimulus(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
        
    [sf, stim_idx, chan_idx ] = get_baphy_plot_controls(stack);
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    filt_idx = get(filt_pop, 'Value');
       
    dat = x.dat.(sf);
    
    plot(dat.(mdl.output_time), ...
         dat.(mdl.output)(:, stim_idx, chan_idx, filt_idx), ...
         pickcolor(filt_idx));
    axis tight;
    drawnow;
end