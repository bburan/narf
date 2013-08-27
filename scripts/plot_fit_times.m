function plot_fit_times(batch, cellids, modelnames)

    function ret = getit(x, idx)
        if isempty(x)
            ret = nan;
        else
            ret = x(idx);
        end
    end

    function ret = getfitval(x)
        if ~isfield(x, 'fit_time') || ~isfield(x, 'perf_val_corr')
            ret = [nan nan];
        else
            ret = [x.fit_time x.perf_val_corr];
        end
    end

stack_extractor = @(x) 0;
meta_extractor = @getfitval;

[~, metas, ~] = load_model_batch(batch, cellids, modelnames, ...
                       stack_extractor, meta_extractor);

times = cellfun(@(x) getit(x,1), metas);
corrs = cellfun(@(x) getit(x,2), metas);

figure('Name', 'Fit Time vs Corr', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
plot(times', corrs', '.');
h = legend(modelnames);
set(h, 'Interpreter', 'none');
xlabel('Fit Time');
ylabel('Val. Set Correlation');
title('Val. Set Correlation vs Fit Time');
n = ceil(sqrt(size(times,1)));
figure;
for ii = 1:size(times, 1)
    subplot(n,n,ii);
	hist(times(ii,:), ceil(size(times,2)/5));
    title([modelnames{ii} '(' num2str(nanmean(times(ii))) ')'], 'Interpreter', 'none');
end

end