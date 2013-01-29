function plot_summary_heatmaps(summary_files)
% A high level function which creates several PNG summary plots extracted
% from the cellid summary files. A good example of how to use
% select_summaries.m to create useful plots.
% 
% ARGUMENTS:
%    SUMMARY_FILES    A cell array of files to use

global NARF_SAVED_ANALYSIS_PATH;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M, xl, yl] = select_summaries(summary_files, ...
                @(c) 100*getfield(c, 'score_test_corr'));
figure; clf;
heatmap(M, xl, yl, '%2.0f', 'TickAngle', 90,...
        'ShowAllTicks', true, 'TickFontSize', 6);
set(gca,'Position',[.05 .2 .9 .75])
title('Test Set Correlation Squared (r^2)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M, xl, yl] = select_summaries(summary_files, ...
                @(c) 100*getfield(c, 'score_test_corr'), ...
                @(c) strfind(c.modelname, 'log'), ...
                @(cc) sprintf('%f', nanmean(cell2mat(extract_field(cc, 'score_test_corr')))));
            
figure; clf;
heatmap(M, xl, yl, '%2.0f', 'TickAngle', 90,...
        'ShowAllTicks', true, 'TickFontSize', 6);
set(gca,'Position',[.05 .2 .9 .75])
title('Test Set Correlation Squared (r^2)');

% [M, xl, yl] = select_summaries(summary_files, ...
%                 @(c) 100*getfield(c, 'score_train_corr'), ...
%                 @(c) strfind(c.modelname, 'log'), ...
%                 @(cc) nanmean(cell2mat(extract_field(cc, 'score_train_corr'))), ...
%                 @(cc) cc{1}.score_test_corr, ...
%                 'score_train_corr', ...
%                 'score_test_corr');
%             
% figure; clf;
% plot(M, 'k.')
% set(gca,'Position',[.05 .2 .9 .75])
% title('Test Set Correlation Square (r^2)');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M, xl, yl] = select_summaries(summary_files, ...
%                 @(c) 100*getfield(c, 'score_train_corr'));
% figure; clf;
% heatmap(M, xl, yl, '%2.0f', 'TickAngle', 90,...
%         'ShowAllTicks', true, 'TickFontSize', 6);
% set(gca,'Position',[.05 .2 .9 .75])
% title('Train Set Correlation Squared (r^2)');

% 
% 
% pngfile = [NARF_SAVED_ANALYSIS_PATH filesep cellid '.png'];
% if savetodisk
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf,'InvertHardcopy','off');
%     print(fh, pngfile, '-dpng');    
%     close(fh);
%     fh = nan;
% end
% 
% % Test performance heatmap
% % Training performance heatmap
% % Fitting time heatmap
% % Exit code heatmap
% % Test/Train scatter grouped by fitter 
% 
% 
% 
% 
% plot_saved_analyses('Test Corr', @(c) 100*getfield(c, 'score_test_corr'), ...
%                     summary_files);
% 
% plot_saved_analyses('Train Corr', @(c) 100*getfield(c, 'score_train_corr'), ...
%                     summary_files, modelnames);
% 
% plot_saved_analyses('Fit Time', @(c) getfield(c, 'fit_time'), ...
%                     summary_files, modelnames); 
% 
% plot_saved_analyses('Exit Code', @(c) getfield(c, 'exit_code'), ...
%                     summary_files, modelnames); 
% 
% 
% % Make another figure with models ranked by average performance
% M2 = nanmean(M, 1);
% Mc = num2cell(M2');
% 
% for ii = 1:length(Mc);
%     Mc{ii} = {Mc{ii}, xlabs{ii}, {M(:, ii)}};
% end
% 
% [vals,order] = sort(cellfun(@(v) v{1}, Mc));
% sortedM = Mc(order);
% sxl = cellfun(@(v) v{2}, sortedM, 'UniformOutput', false);
% sM = cellfun(@(v) v{3}, sortedM);
% sM = cat(2, sM{:});
% figure; clf;
% heatmap(sM, sxl, ylabs, '%2.0f', 'TickAngle', 90,...
%         'ShowAllTicks', true, 'TickFontSize', 6);
% set(gca,'Position',[.05 .2 .9 .75])
% title(sprintf('%s, sorted by column nanmean()', thetitle));
% 
% % Also, display a performance plot for each model token
% tokens = cellfun(@(l) regexp(l, '(.*?)(?:_|$)', 'tokens'), xlabs, 'UniformOutput', false);
% tokens = cat(1, tokens{:});
% tokens = unique(cat(2, tokens{:}));
% mu = [];
% top10 = [];
% figure;
% hold on;
% for ii = 1:length(tokens)
%     tok = tokens{ii};
%     mods = cellfun(@(x)~isempty(x), regexp(xlabs, tok));
%     tmp =  M(:, mods);
%     scores = tmp(:);
%     sorted_scores = sort(scores(:));
%     sorted_scores(any(isnan(sorted_scores),2),:)=[];
%     mu(ii) = nanmean(scores(:));
%     top10(ii) = nanmean(sorted_scores(max(1,floor(end-length(sorted_scores)/10)):end));
%     plot(ii*ones(1, length(scores))+0.01*randn(1, length(scores)), scores, 'k.');
% end
% plot(1:length(tokens), mu, 'b-');
% plot(1:length(tokens), top10, 'r-');
% set(gca,'ButtonDownFcn','selectmoveresize');
% set(gca,'Position',[.05 .05 .9 .9])
% hold off;
% title(sprintf('Mean (and Top10 Percent) Performance Across All Cellids: %s', thetitle));
% xticks(1:length(tokens), tokens);
