function plot_scatter(X, names)
% CORRPLOT  Plots a big correlation graph of correlations and residuals
%   Each column of X is a different data set that you want to plot
%   names is a cell array, with strings for each.

% 2012/10/03, Ivar Thorson.

% Set plot borders and sizes
tmargin = 0.1;
lmargin = 0.1;
rmargin = 0.1;
bmargin = 0.1;

% Set text position on each subplot, from 0 to 1
textx = 0.05;  % 0.1=Left
texty = 0.95;  % 0.9=Top

axmin = -0.2;
axmax = 0.9;
axmin = min(min(X(:)), 0);
axmax = max(X(:));

bins = linspace(axmin, axmax, 20);

[Npoints, Nsets] = size(X);

if not(exist('names')) || length(names) ~= Nsets
    error('I need as many names as the size of the y dimension!');
end

R = corrcoef(X);

w = (1.0 - lmargin - rmargin) / (Nsets - 1);    % Plot width
h = (1.0 - tmargin - bmargin) / (Nsets - 1);    % Plot height

for i = 1:Nsets
    for j = 1:Nsets;
        
        % Calculate draw location
        if i < j
            subplot('Position', [(lmargin+(i-1)*w) (bmargin + (Nsets-j)*h) w h]);
        else
            continue;
        end
        
        axis([axmin, axmax, axmin, axmax]);
        
        hold on;
        if j == Nsets
             qtys = hist(X(:,i), bins);
             bar(bins, (qtys ./ sum(qtys)), 'g');
        end
%          
        if i == 1
             qtys = hist(X(:,j), bins);
             bh = barh(bins, (qtys ./ sum(qtys)), 'y');
         end
         
        plot(X(:,i),X(:,j), 'k.', [axmin, axmax], [axmin, axmax], 'k--');
% 
%         xtmp = sort(X(:,i));
%         ytmp = sort(X(:,j));
%         
%         cs = 2* cumsum((xtmp > ytmp) - 0.5);
        %plot(linspace(0, 1, length(cs)), cs / length(cs), 'r-');
            %[N,X]
%                 sprintf('r=%f', R(i,j)));
        
        % Turn off tick labels unless in the bottom or leftmost rows
        if j == Nsets 
            hl = xlabel(sprintf('%s\nmean:%.3f med:%.3f', names{i}, ...
                                nanmean(X(:,i)),nanmedian(X(:,i))));
            set(hl,'interpreter','none');
        else
            set(gca, 'XTickLabel', '');
        end
        
        if i == 1 
            hl = ylabel(sprintf('%s\nmean:%.3f med:%.3f', names{i}, ...
                                nanmean(X(:,j)),nanmedian(X(:,j))));
            set(hl,'interpreter','none');
        else
            set(gca, 'YTickLabel', '');
        end
        hold off;
    end
end
