function compare_models(batch, cellids, modelnames)

    function ret = getfitval(x)
        if isfield(x, 'fit_time') && ~isempty(x.fit_time);
            ret = x.fit_time;
        else
            ret = nan;
        end
    end

stack_extractor = @(x) 0;
meta_extractor = @getfitval;

[~, times, ~] = load_model_batch(batch, cellids, modelnames, ...
                       stack_extractor, meta_extractor);

times(cellfun(@iscell, times)) = {nan}
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
B = repmat(1:n, size(data,1), 1);
if PLOT_JITTER
    B = B - 0.1 + 0.2*rand(size(B));
end    
plot(B, data, 'k.')
errorbar(nanmean(data), nanstd(data), 'r.'); 
hold off;

title(sprintf('Average and Std Dev of Model Fit Time, Batch %d', batch));
set(gca,'XTick',1:n);        
set(gca,'XTickLabel', modelnames(idxs));
if PLOT_LOGTIME
    set(gca,'YScale','log');
end
set(gca,'CameraUpVector',[-1,0,0]);

end

