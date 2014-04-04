function plot_npnls(batch, cellids, modelnames)
global STACK META XXX;

[stacks, metas, x0s] = load_model_batch(batch, cellids, modelnames);

figure('Name', 'NPNLs', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
hold on;

for ii = 1:length(stacks)
    STACK = stacks{ii};
    META = metas{ii};
    XXX = x0s(ii);
       
    % If there are no NPFNLs, skip
    [~, mod_idx] = find_modules(STACK, 'nonparm_filter_nonlinearity', true);
    if isempty(mod_idx)
        continue;
    end
    
    % Otherwise, plot the NPFNL
    calc_xxx(1);
    
    m = STACK{mod_idx}{1};
    sel = [];  
    sel.stim_idx = 1;
    sel.chan_idx = 1;
    sel.stimfile = XXX{1}.test_set{1};  
    pfn = m.auto_plot;   
    pfn(sel, STACK(1:(mod_idx)), XXX(1:(mod_idx+1)));
    %keyboard;  
    
end
% keyboard;
% end
% 
% figure('Name', 'Fit Time Bars', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
% if PLOT_SORTED
%     means = nanmean(times);
%     stds = nanstd(times);
%     D = [means' times'];
%     [sD, idxs] = sortrows(D, -1);
%     data = sD(:, 2:end)';
% else
%     data = times;
% end
% 
% hold on;
% B = repmat(1:n, size(data,1), 1);
% if PLOT_JITTER
%     B = B - 0.1 + 0.2*rand(size(B));
% end
% plot(B, data, 'k.')
% errorbar(nanmean(data), nanstd(data), 'r.');
% hold off;
% 
% title(sprintf('Average and Std Dev of Model Fit Time, Batch %d', batch));
% set(gca,'XTick',1:n);
% set(gca,'XTickLabel', modelnames(idxs));
% if PLOT_LOGTIME
%     set(gca,'YScale','log');
% end
% set(gca,'CameraUpVector',[-1,0,0]);

end

