function do_plot_channels_as_heatmap(sel, stack, xxx)

    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end));  
    
    ii=1;
    while ~isfield(xouts{ii}.dat,sel.stimfile) && ii<length(xouts),
       ii=ii+1;
    end
    
    if isfield(mdls{1},'output'),
        outfield=mdls{1}.output;
    elseif isfield(mdls{1},'outputs'),
        for jj=1:length(mdls{1}.outputs),
            if ~isempty(xouts{ii}.dat.(sel.stimfile).(mdls{1}.outputs{jj})),
                outfield=mdls{1}.outputs{jj};
            end
        end
    end
    if size(xouts{ii}.dat.(sel.stimfile).(outfield),2)>=sel.stim_idx,
        h = imagesc(xouts{ii}.dat.(sel.stimfile).(mdls{1}.time)(:),...
                    1:size(xouts{ii}.dat.(sel.stimfile).(outfield),3),...
                    squeeze(xouts{ii}.dat.(sel.stimfile).(outfield)(:, sel.stim_idx, :))');
        do_xlabel('Time [s]');
        do_ylabel('Output [-]');
        set(gca,'YDir','normal');
        axis xy tight;
    end
    
    
end