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
        d = sign(dd) .* sqrt(0.5 * abs(dd).^2);
        lim = max(abs(d));
        xlim([-lim lim]);
        histfit(d, ceil(Npoints/5));
        
        hh = get(gca,'Children');
        set(hh(2),'FaceColor',[.8 .8 .9]);
        
        % set(gca,'CameraUpVector',[-1,0,1]);

        % Test if it's a gaussian
        [hchi, pchi] = chi2gof(d);
        if hchi, isgauss = 'true'; else isgauss = 'false'; end
        
        % Test if the gaussian is zero mean
        [htt, ptt] = ttest(d);
        bias = nanmean(d);
        if htt, 
            istt = 'true';
            if bias >= 0
                winner = names{i};
            else
                winner = names{j};
            end
        else
            istt = 'false';             
            winner = 'Uncertain';            
        end
        
        textLoc(sprintf('Gaussian? %5s (Chi-Square p-value: %f)\nMean not 0? %5s (T-test p-value: %f)\nMean: %f\nWinner: %s', ...
            isgauss, pchi, istt, ptt, bias, winner), ...
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
