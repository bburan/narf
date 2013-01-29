function [M, xlabs, ylabs] = plot_saved_analyses(thetitle, thefn, summary_files, maskfn)
% Plot a heat map of every saved analysis file given as an argument
% The value of the heat map is equal to the return value of THEFN, which is
% applied to each element of the 'summary' cell array. Assumes that the
% value of THEFN is a number between 0 and 100.
% 
% The heatmap has the following properties:
%    On the X axis are the models, sorted alphabetically
%    On the Y axis are the cellids, sorted alphabetically
%    If a model exists for one cellid but not others, the others will take on
% values of NaN in the heatmap.
% 
% ARGUMENTS:
%     THETITLE         Title of the heat map plot
%     THEFN            Function to apply to each saved model
%     SUMMARY_FILES    Which cached summary files to use
%     MASK             Optional. If present, ONLY modelnames present in
%                      MASK will be used to build the heatmap.
%
% RETURNS:
%     the matrix, x labels, and y labels used in the plot

H = []; % A big struct of structs to save data in

if nargin < 4
    mask = {};
end

% Put all analysis summaries into a big struct array
for ii = 1:length(summary_files);
    af = summary_files{ii};
    
    % Load the summary of the analysis
    if exist(af, 'file') ~= 2
        fprintf('Skipping nonexistant summary file: %s\n', af);
        continue;
    end
    summary = getfield(load(af, 'summary'), 'summary');

    % Rip out all the analysis summaries and store them in the big struct
    for ii = 1:length(summary)
        s = summary{ii};
        if isempty(s) 
            continue;
        end
        if isempty(mask) | any(strcmp(s.modelname, mask))
            cellid = regexprep(s.cellid, '-', '_');
            H.(cellid).(s.modelname) = thefn(s);
        end
    end
end

% Convert the struct of structs into a sorted matrix and sorted labels
ylabs = sort(fieldnames(H));
vals = cellfun(@(f) getfield(H, f), fieldnames(H), 'UniformOutput', false);
keys = cellfun(@fieldnames, vals, 'UniformOutput', false);
xlabs = sort(unique(cat(1, keys{:})));

M = NaN * zeros(length(ylabs), length(xlabs));
for ii = 1:length(ylabs),
    for jj = 1:length(xlabs),
        if isfield(H.(ylabs{ii}), xlabs{jj})
            if isnan(H.(ylabs{ii}).(xlabs{jj}))
                M(ii,jj) = nan;
            else
                M(ii,jj) = H.(ylabs{ii}).(xlabs{jj});
            end
        end
    end
end

figure; clf;
heatmap(M, xlabs, ylabs, '%2.0f', 'TickAngle', 90,...
        'ShowAllTicks', true, 'TickFontSize', 6);
set(gca,'Position',[.05 .2 .9 .75])
title(sprintf('%s', thetitle));


% Make another figure with models ranked by average performance
M2 = nanmean(M, 1);
Mc = num2cell(M2');

for ii = 1:length(Mc);
    Mc{ii} = {Mc{ii}, xlabs{ii}, {M(:, ii)}};
end

[vals,order] = sort(cellfun(@(v) v{1}, Mc));
sortedM = Mc(order);
sxl = cellfun(@(v) v{2}, sortedM, 'UniformOutput', false);
sM = cellfun(@(v) v{3}, sortedM);
sM = cat(2, sM{:});
figure; clf;
heatmap(sM, sxl, ylabs, '%2.0f', 'TickAngle', 90,...
        'ShowAllTicks', true, 'TickFontSize', 6);
set(gca,'Position',[.05 .2 .9 .75])
title(sprintf('%s, sorted by column nanmean()', thetitle));


% Also, display a performance plot for each model token
tokens = cellfun(@(l) regexp(l, '(.*?)(?:_|$)', 'tokens'), xlabs, 'UniformOutput', false);
tokens = cat(1, tokens{:});
tokens = unique(cat(2, tokens{:}));
mu = [];
top10 = [];
figure;
hold on;
for ii = 1:length(tokens)
    tok = tokens{ii};
    mods = cellfun(@(x)~isempty(x), regexp(xlabs, tok));
    tmp =  M(:, mods);
    scores = tmp(:);
    sorted_scores = sort(scores(:));
    sorted_scores(any(isnan(sorted_scores),2),:)=[];
    mu(ii) = nanmean(scores(:));
    top10(ii) = nanmean(sorted_scores(max(1,floor(end-length(sorted_scores)/10)):end));
    plot(ii*ones(1, length(scores))+0.01*randn(1, length(scores)), scores, 'k.');
end
plot(1:length(tokens), mu, 'b-');
plot(1:length(tokens), top10, 'r-');
set(gca,'ButtonDownFcn','selectmoveresize');
set(gca,'Position',[.05 .05 .9 .9])
hold off;
title(sprintf('Mean (and Top10 Percent) Performance Across All Cellids: %s', thetitle));
xticks(1:length(tokens), tokens);
