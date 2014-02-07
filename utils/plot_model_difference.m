function plot_model_difference(X, names)
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

% Shrink the names to their minimum
names = shorten_modelnames(names);

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
        
        
        hold on;
        dd = X(:,i) - X(:,j);
        % d = sign(dd) .* sqrt(0.5 * abs(dd).^2);
        d = dd;
        lim = max(abs(d));
        xlim([-lim lim]);
        histfit(d, ceil(Npoints/10));        
        
        % set(gca,'CameraUpVector',[-1,0,1]);

        % Test if it's a gaussian
        [hchi, pchi] = chi2gof(d);
        if ~hchi,
            isgauss = 'Probably';
            saturation = 0;
        else
            isgauss = 'false';
            saturation = 0.8;
        end
        
        % Test if the gaussian is zero mean
        [htt, ptt] = ttest(d);
        bias = nanmean(d);
        if htt, 
            istt = 'Probably';
            if bias >= 0
                winner = names{i};
                barcolor = [saturation 1 saturation];
            else
                winner = names{j};
                barcolor = [1 saturation saturation];
            end
        else
            istt = 'false';             
            winner = 'Uncertain';    
            barcolor = [0.8 0.8 0.8];
        end

        % Change the color of the histogram to indicate the winner
        hh = get(gca,'Children');
        set(hh(1),'Color', [0 0 0]);
        set(hh(2),'FaceColor', barcolor);
        v=axis;
        plot([0,0], [0, v(4)], 'k--');
        
        %textLoc(sprintf('Gaussian? %5s (Chi-Square p-value: %f)\nNonzero mean? %5s (T-test p-value: %f)\nMean: %f\nWinner: %s', ...
        %    isgauss, pchi, istt, ptt, bias, winner), ...
        %    'NorthWest', 'interpreter', 'none');
        
        textLoc(sprintf('Winner: %s', winner), ...
            'NorthWest', 'interpreter', 'none');
        
        % Turn off tick labels unless in the bottom or leftmost rows
        if j == Nsets 
            hl = xlabel(sprintf('%s\nmean:%f', names{i}, nanmean(X(:,i))));
            set(hl,'interpreter','none');
        else
            set(gca, 'XTickLabel', '');
        end
        
        if i == 1 
            hl = ylabel(sprintf('%s\nmean:%f', names{j}, nanmean(X(:,j))));
            set(hl,'interpreter','none');
        else
            set(gca, 'YTickLabel', '');
        end
        hold off;
    end
end
