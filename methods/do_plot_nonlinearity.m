function do_plot_nonlinearity(stack, xxx, thesig, thefn, showhist)
    mdl = stack{end};
    x_pre = xxx{end-1};
    x_post = xxx{end};
    
    [sf, stim_idx, unused] = get_baphy_plot_controls(stack);  
    
    [bins, centers] = hist(x_pre.dat.(sf).(thesig)(:, stim_idx), 50);
    xs = linspace(centers(1), centers(end), 200);
    
    if ~showhist
        plot(xs, thefn(xs), 'k-');
        axis tight;
    else
        hold on;
        bar(centers, bins);
        max_bins = max(bins);
        max_fn = max(thefn(xs));
        plot(xs, (max_bins/max_fn)*thefn(xs)); 
        axis([xs(1) xs(end) 0 max_bins]); 
        hold off;
        % Plot YY would be better here, but it is buggy because it will
        % try to look for handles automatically and fail.
        %         [AX, H1, H2] = plotyy(centers, bins, xs, thefn(xs), 'bar', 'plot');
        %         set(AX(2),'XTick',[]); % Don't display xticks
        %         axis(AX(1), [xs(1) xs(end) 0 max(bins(2:end))]); % Ignore zero bin
        %         % axis(AX(2), tight);
    end
    
end