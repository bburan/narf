function fh = corrplot(X, names)
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
textx = 0.1;  % 0.1=Left
texty = 0.9;  % 0.9=Top

[Npoints, Nsets] = size(X);

if not(exist('names'))
    names = [1:Nsets];
end

R = corrcoef(X);

w = (1.0 - lmargin - rmargin) / Nsets;    % Plot width
h = (1.0 - tmargin - bmargin) / Nsets;    % Plot height

for i= 1:Nsets
    for j = 1:Nsets;
        
        % Calculate draw location
        subplot('Position', [(lmargin+(i-1)*w) (bmargin + (Nsets-j)*h) w h]);
        
        if i == j                   % Plot histograms on the diagonal
            hist(X(:,i), 20);             
            
        elseif i < j                % Plot correlations below the diagonal
            plot(X(:,i),X(:,j), 'k.');
            xl = xlim();
            yl = ylim();
            text(textx * (xl(2) - xl(1)) + xl(1), ...
                texty * (yl(2) - yl(1)) + yl(1), ...
                sprintf('r=%f', R(i,j)));
            
        else                        % Plot residuals above the diagonal
            % plot( 
            % TODO
        end       
        
        % Turn off tick labels unless in the bottom or leftmost rows, and
        % make the labels in those cases dynamic and active
        if j == Nsets 
            xlabel(names(i));
        else
            set(gca, 'XTickLabel', '');
        end
        
        if i == 1
            ylabel(names(j));
        else
            set(gca, 'YTickLabel', '');
        end
    end
end
