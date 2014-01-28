function plot_fit_time_histogram(batch, cellids, modelnames)

% CHANGE THESE VALUES HERE TO CHANGE THE PLOT AS DESIRED
PLOT_SORTED = true;
PLOT_LOGTIME = true;

    function ret = getfitval(x)
        if ~isfield(x, 'fit_time')
            ret = [nan nan];
        else
            ret = x.fit_time;
        end
    end

stack_extractor = @(x) 0;
meta_extractor = @getfitval;

[~, times, ~] = load_model_batch(batch, cellids, modelnames, ...
                       stack_extractor, meta_extractor);

times = cell2mat(times)';

figure('Name', 'Fit Time Histogram', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
n = size(times,2);
for ii = 1:n
    subplot(n,1,ii);
    hist(times(:,ii), 50);   
    ll = xlabel(sprintf('Fit Time for %s', modelnames{ii}));
    set(ll, 'Interpreter', 'none');
end

figure('Name', 'Fit Time Bars', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
if PLOT_SORTED
    means = nanmean(times);
    stds = nanstd(times);
    D = [means' times'];
    [sD, idxs] = sortrows(D, -1);
    data = sD(:, 2:end)';
else
    data = times;
end

hold on; 
%bar(1:n, nanmean(data), 0.3, 'w'); 
%         plot(1:n, data, 'k.'); 
plot(1:n, data, 'k.');
errorbar(nanmean(data), nanstd(data), 'r.'); 
hold off;

set(gca,'XTick',1:n);        
set(gca,'XTickLabel', modelnames(idxs));
if PLOT_LOGTIME
    set(gca,'YScale','log');
end
set(gca,'CameraUpVector',[-1,0,0]);

title(sprintf('Average and Std Dev of Model Fit Time, Batch %d', batch));

end

