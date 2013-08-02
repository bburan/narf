function plot_fit_times(batch, cellids, modelnames)

stack_extractor = @(x) 0;
meta_extractor = @(x) [x.fit_time x.perf_val_corr];

[~, metas, ~] = load_model_batch(batch, cellids, modelnames, ...
                       stack_extractor, meta_extractor);
                                    
times = cellfun(@(x) x(1), metas);
corrs = cellfun(@(x) x(2), metas);

figure('Name', 'Fit Time vs Corr', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
%hold on;
%for ii = 1:size(times, 1)
%    p = pickpoint(ii);
%	plot(times(ii,:), corrs(ii,:), p);
%end
%hold off;
plot(times', corrs', '.');
h = legend(modelnames);
set(h, 'Interpreter', 'none');
xlabel('Fit Time');
ylabel('Val. Set Correlation');
title('Val. Set Correlation vs Fit Time');
