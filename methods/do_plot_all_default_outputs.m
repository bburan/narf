function do_plot_all_default_outputs(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    
    for jj=1:length(mdls),
        if ~isempty(xouts{jj}.dat.(sel.stimfile).(mdls{jj}.output)),
            outfield=mdls{jj}.output;
        end
    end
    
    do_plot(xouts, mdls{1}.time, outfield, ...
            sel, 'Time [s]', 'Output [-]');
end
