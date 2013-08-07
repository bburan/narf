function plot_selfdist_vs_corr(batch, cellids, modelnames)

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
            ret = [x.metric_est_selfdist x.metric_val_selfdist x.perf_est_corr x.perf_val_corr];
        end
    end

stack_extractor = @(x) 0;
meta_extractor = @getfitval;

[~, metas, ~] = load_model_batch(batch, cellids, modelnames, ...
                       stack_extractor, meta_extractor);

est_selfdist = cellfun(@(x) getit(x,1), metas);
val_selfdist = cellfun(@(x) getit(x,2), metas);
est_corr = cellfun(@(x) getit(x,3), metas);
val_corr = cellfun(@(x) getit(x,4), metas);

figure('Name', 'Correlation vs Self Dist', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
plot(est_selfdist, est_corr, 'x', val_selfdist, val_corr, 'o');
h = legend('Est. Set', 'Val. Set');
set(h, 'Interpreter', 'none');
xlabel('Cell''s RESP Self-Dist');
ylabel('Correlation');
title('Correlation vs Self-Dist');

end