function [M, xlabs, ylabs] = plot_saved_analyses(thefn, analysis_files)
% Plot a heat map of every saved analysis file given as an argument
% The value of the heat map is equal to the return value of THEFN, which is
% applied to each element of the 'results' cell array. Assumes that the
% value of THEFN is a number between 0 and 100, such as 100/r^2, where r is
% the correlation coefficient.
% 
% The heatmap has the following properties:
%    On the X axis are the models, sorted alphabetically
%    On the Y axis are the cellids, sorted alphabetically
%    If a model exists for one cellid but not others, the others will take on
% values of NaN in the heatmap.
% 
% Returns the matrix, x labels, and y labels used in the plot

H = []; % A big struct of structs to save data in

% Put all analysis results into a big struct array
for ii = 1:length(analysis_files);
    af = analysis_files{ii};
    
    % Load the results of the analysis
    results = getfield(load(af, 'results'), 'results');
    
    % Rip out all the analysis results and store them in the big struct
    for ii = 1:length(results)
        res = results{ii};
        if isempty(res) 
            continue;
        end
        modelname = res.modelname;
        cellid = regexprep(res.cellid, '-', '_');
        H.(cellid).(modelname) = thefn(res);
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
            M(ii,jj) = H.(ylabs{ii}).(xlabs{jj});
        end
    end
end

figure; clf;
heatmap(M, xlabs, ylabs, '%2.0f', 'TickAngle', 90,...
        'ShowAllTicks', true, 'TickFontSize', 6);

% Also, display a mean/var plot for each model token between underscores
tokens = cellfun(@(l) regexp(l, '_(.*?)_', 'tokens'), xlabs, 'UniformOutput', false);
tokens = cat(1, tokens{:});
tokens = unique(cat(2, tokens{:}));
mu = [];
figure;
hold on;
for ii = 1:length(tokens)
    tok = tokens{ii};
    mods = cellfun(@(x)~isempty(x), regexp(xlabs, tok));
    tmp =  M(:, mods);
    scores = tmp(:);
    mu(ii) = nanmean(scores(:));
    plot(ii*ones(1, length(scores)), scores, 'k.');
%     s_min(ii) = nanmin(scores(:));
%     s_max(ii) = nanmax(scores(:));
end
plot(1:length(tokens), mu, 'b-');
hold off;
title('Mean Model Performance Across All Cellids');
xticks(1:length(tokens), tokens);
