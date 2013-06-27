function plot_summaries(summaries, prefix, savetodisk)
% A high level function which creates several PNG summary plots extracted
% from multiple cellid summary files. A good example of how to use
% select_summaries.m and heat maps to create useful plots.
% 
% ARGUMENTS:
%    SUMMARY_FILES    A cell array of files to use

global NARF_SAVED_ANALYSIS_PATH;

if nargin < 3
    savetodisk = false;
end

function saveit (fh, filename)
    if savetodisk
        savethefig(fh, filename, prefix);
    end
end

function fh = phm(M, xl, yl, the_title, filename) 
    fh = figure('Position', [0 0 2000 1200]); clf;
    heatmap(M, xl, yl, '%2.0f', 'TickAngle', 90,...
                'ShowAllTicks', true, 'TickFontSize', 12);
    set(gca,'Position',[.05 .2 .9 .75])
    title(the_title);
end


% ------------------------------------------------------------------------
% Abscissa and Ordinate sorted by value

[M, xl, yl] = select_summaries(summaries, ...
                @(c) 100*getfield(c, 'score_test_corr'), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))));
fh = phm(M, xl, yl, 'Test Set r^2, Value-Sorted Abscissa and Ordinate');
saveit(fh, 'heat0.png');

% ------------------------------------------------------------------------
% Abscissa sorted alphabetically by 1st Token (Envelope)

[M, xl, yl] = select_summaries(summaries, ...
                @(c) 100*getfield(c, 'score_test_corr'));

fh = phm(M, xl, yl, 'Test Set r^2, Alphabetical Abscissa and Ordinate');
saveit(fh, 'heat1.png');

% ------------------------------------------------------------------------
% Abscissa sorted alphabetically by 2nd Token (Usually compressor)

[M, xl, yl] = select_summaries(summaries, ...
                @(c) 100*getfield(c, 'score_test_corr'), ...
                @(cc) rotate_tokens(cc{1}.modelname, 2), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))));
fh = phm(M, xl, yl, 'Test Set r^2, Abscissa Sorted by 2nd Token, Value-Sorted Ordinate');
saveit(fh, 'heat2.png');

% ------------------------------------------------------------------------
% Abscissa sorted alphabetically by 3rd Token (Usually FIR type)
[M, xl, yl] = select_summaries(summaries, ...
                @(c) 100*getfield(c, 'score_test_corr'), ...
                @(cc) rotate_tokens(cc{1}.modelname, 3), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))));
fh = phm(M, xl, yl, 'Test Set r^2, Abscissa Sorted by 3rd Token, Value-Sorted Ordinate');
saveit(fh, 'heat2.png');

% ------------------------------------------------------------------------
% Abscissa sorted alphabetically by 4th Token (Usually Nonlinearity)
[M, xl, yl] = select_summaries(summaries, ...
                @(c) 100*getfield(c, 'score_test_corr'), ...
                @(cc) rotate_tokens(cc{1}.modelname, 4), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))));
fh = phm(M, xl, yl, 'Test Set r^2, Abscissa Sorted by 4th Token, Value-Sorted Ordinate');
saveit(fh, 'heat3.png');

% ------------------------------------------------------------------------
% Abscissa sorted alphabetically by 5th Token (Usually Fitter)
[M, xl, yl] = select_summaries(summaries, ...
                @(c) 100*getfield(c, 'score_test_corr'), ...
                @(cc) rotate_tokens(cc{1}.modelname, 5), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))));
fh = phm(M, xl, yl, 'Test Set r^2, Abscissa Sorted by 5th Token, Value-Sorted Ordinate');
saveit(fh, 'heat4.png');

% ------------------------------------------------------------------------
% Sorted by Fitting time
[M, xl, yl] = select_summaries(summaries, ...
                @(c) getfield(c, 'fit_time'), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'fit_time')))));
fh = phm(M, xl, yl, 'Fitting Time, Value-Sorted Abscissa');
saveit(fh, 'heat_time.png');

% ------------------------------------------------------------------------
% Exit code heatmap
[M, xl, yl] = select_summaries(summaries, ...
                @(c) getfield(c, 'exit_code'), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))), ...
                @(cc) sprintf('%.9f', nanmean(cell2mat(getfieldforcellarray(cc, 'score_test_corr')))));            
fh = phm(M, xl, yl, 'Exit Code, Sorted Time-Sorted Abscissa');
saveit(fh, 'heat_exit.png');

% ------------------------------------------------------------------------
% Labeled scatter plot?

test_scores = cell2mat(getfieldforcellarray(summaries, 'score_test_corr'));
train_scores = cell2mat(getfieldforcellarray(summaries, 'score_train_corr'));
fh = figure;
plot(train_scores, test_scores, 'k.');
xlabel('Training Set r^2');
ylabel('Test Set r^2');
title('Training/Test Performance for all Models and Cellids');
saveit(fh, 'test_train.png');

% ------------------------------------------------------------------------
% Average Token Performance

% tokens = cellfun(tokenize_modelname(c.modelname), S, 'UniformOutput', false);
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
% title(sprintf('Mean Performance Across All Cellids: %s', thetitle));
% xticks(1:length(tokens), tokens);

% ------------------------------------------------------------------------

end

