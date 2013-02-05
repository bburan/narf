function fh = plot_token_performance(summaries, titlename, tokidx, tokens)
% Plots performance of model tokens across all the summary files given.
 
S = summaries;

colormap = [];
for ii = 1:length(tokens)
    colormap.(tokens{ii}) = ii;
end

function g = pickgroup(modelname)
    toks = tokenize_modelname(modelname); 
    k = toks{tokidx};
    g = colormap.(k{:});    
end

test_scores = cell2mat(extract_field(S, 'score_test_corr'));
train_scores = cell2mat(extract_field(S, 'score_train_corr'));
token_group  = cell2mat(cellfun(@pickgroup, ...
                                 extract_field(S, 'modelname'),...
                                 'UniformOutput', false));

fh = figure;
gscatter(train_scores, test_scores, token_group, 'bgrcmky', 'o');
xlabel('Training Set r^2');
ylabel('Test Set r^2');
title(['Training/Test Performance for ' titlename]);
legend(tokens, 'Location', 'SouthEast');
end