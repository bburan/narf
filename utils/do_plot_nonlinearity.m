function do_plot_nonlinearity(stack, xxx, thesig, thefn, showhist)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    [bins, centers] = hist(dat.(thesig)(:, stim_idx), 50);
    xs = linspace(centers(1), centers(end), 200);
    
    if ~showhist
        plot(xs, thefn(xs), 'k-');
        axis tight;
    else
        [AX, H1, H2] =plotyy(centers, bins, xs, thefn(xs), 'bar', 'plot');
        set(AX(2),'XTick',[]); % Don't display xticks
        axis(AX(1), [xs(1) xs(end) 0 max(bins(2:end))]); % Ignore zero bin
        % axis(AX(2), tight);
    end
    
end