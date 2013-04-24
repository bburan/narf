function do_plot_scatter(sel, xxxs, field1, field2, n_plotpts)
% Plots two fields vs each other 
% n_plotpts determines how many points will be plotted. 
% If not provided, all points will be plotted. 
% If you want a smoother plot, try n_plotpts=100

hold on;

% These variables will hold the limits of ALL plots
xmin = nan;
xmax = nan;
ymin = nan;
ymax = nan;

for ii = 1:length(xxxs)
    x = xxxs{ii}{end};
    
    if ~isfield(x.dat, sel.stimfile)
        % Do nothing if there is no matching selected stimfile
        continue;
    end
    
    dat = x.dat.(sel.stimfile);
    
    if ~isequal(size(dat.(field1)), size(dat.(field2)))
        cla;
        text(0.35, 0.5, 'ERROR! Inputs of unequal size cannot be scatter plotted.');
        axis([0, 1, 0 1]);
        return
    end
    
    D = [dat.(field1)(:) dat.(field2)(:)];
    D = D(~isnan(D(:,1)) & ~isnan(D(:,2)),:);
    
    if exist('n_plotpts', 'var') && ~isnan(n_plotpts) && n_plotpts > 0
        % Sort and average them by groups of 100
        D = sortrows(D);
        D = conv_fn(D, 1, @nanmean, ceil(size(D, 1)/n_plotpts), 0);
    end
       
    c = pickcolor(ii);
    if length(xxxs) == 1
        c = [0 0 0];
    end
    
    plot(D(:,1), D(:,2), 'Color', c, 'Marker', '.', 'LineStyle', 'none');
    
    axis tight;
    v = axis;
    xmin = min(xmin, v(1));
    xmax = max(xmax, v(2));
    ymin = min(xmin, v(3));
    ymax = max(xmax, v(4));
    
end

axis([xmin, xmax, ymin, ymax]);
hold off;