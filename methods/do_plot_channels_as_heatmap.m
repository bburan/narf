function do_plot_channels_as_heatmap(sel, stack, xxx)

    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));  
             
    ii=1;
    while ~isfield(xouts{ii}.dat,sel.stimfile) && ii<length(xouts),
       ii=ii+1;
    end
    h = imagesc(xouts{ii}.dat.(sel.stimfile).(mdls{1}.time)(:),...
                1:size(xouts{ii}.dat.(sel.stimfile).(mdls{1}.output),3),...
                squeeze(xouts{ii}.dat.(sel.stimfile).(mdls{1}.output)(:, sel.stim_idx, :))');
    do_xlabel('Time [s]');
    do_ylabel('Output [-]');
    
    set(gca,'YDir','normal');
    axis xy tight;
    
end