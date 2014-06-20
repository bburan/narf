function do_plot_all_default_outputs(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Output [-]');
end
