function do_plot_scatter(sel, mdl, x, field1, field2, forcelinestyle, n_plotpts)
% Plots two fields vs each other 
% n_plotpts determines how many points will be plotted. 
% If not provided, all points will be plotted. 
% If you want a smoother plot, try n_plotpts=100

if ~isfield(x.dat, sel.stimfile)
    % Do nothing if there is no matching selected stimfile
    return;
end

dat = x.dat.(sel.stimfile);  
    
if ~isequal(size(dat.(field1)), size(dat.(field2)))
    text(0.35, 0.5, 'Cannot scatter plot a multi-channel input.');
    axis([0, 1, 0 1]);
    return;
end    

% Sort and average them by groups of 100
D = [dat.(field1)(:) dat.(field2)(:)]; 
D=D(~isnan(D(:,1))&~isnan(D(:,2)),:);
D = sortrows(D);
if exist('n_plotpts', 'var')
    D = conv_fn(D, 1, @nanmean, ceil(size(D, 1)/n_plotpts), 0);
end

if exist('forcelinestyle', 'var')
    forcelinestyle = 'k.';
end

plot(D(:,1), D(:,2), forcelinestyle);

axis tight;
