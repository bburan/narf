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
    x = xxxs{ii};
    
    if ~isfield(x.dat, sel.stimfile)
        % Do nothing if there is no matching selected stimfile
        continue;
    end
    
    dat = x.dat.(sel.stimfile);
    [T1, S1, C1] = size(dat.(field1));
    [T2, S2, C2] = size(dat.(field2));
    
    if ~isequal([T1 S1], [T2 S2])
        cla;
        text(0.35, 0.5, 'ERROR! Inputs of unequal size cannot be scatter plotted.');
        axis([0, 1, 0 1]);
        return
    end
    
    if C2 ~= 1
        error('I can only accomodate a single channel in input2');
    end
    
    for cc = 1:C1
        tmp = dat.(field1)(:,:,cc);
        D = [tmp(:) dat.(field2)(:)];
        D = D(~isnan(D(:,1)) & ~isnan(D(:,2)),:);
    
        if exist('n_plotpts', 'var') && ~isnan(n_plotpts) && n_plotpts > 0        
            D = sortrows(D);
            D = conv_fn(D, 1, @nanmean, max(1,ceil(size(D, 1)/n_plotpts)), 0);
        end
    
        c = pickcolor(cc);
        if C1 == 1
            c = [0 0 0];
        end
    
        plot(D(:,1), D(:,2), 'Color', c, 'Marker', '.', 'LineStyle', 'none');
    end
    axis tight;
    v = axis;
    xmin = min(xmin, v(1));
    xmax = max(xmax, v(2));
    ymin = min(xmin, v(3));
    ymax = max(xmax, v(4));
    
end

axis([xmin, xmax, ymin, ymax]);
hold off;