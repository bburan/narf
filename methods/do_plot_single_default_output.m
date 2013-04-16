function do_plot_single_default_output(sel, stack, xxx)
    [mdls, xins, xouts] = do_calc_paramsets(stack, xxx);    
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Output [-]');
end
