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

    function mat = get_as_mat(s, fieldname)
        mat = cell2mat(getfieldforcellarray(s, fieldname));
    end

test_scores = get_as_mat(S, 'scorr_test_corr');
train_scores = get_as_mat(S, 'score_train_corr'); 
token_group  = cell2mat(cellfun(@pickgroup, ...
                                 getfieldforcellarray(S, 'modelname'),...
                                 'UniformOutput', false));

fh = figure;
gscatter(train_scores, test_scores, token_group, 'bgrcmky', 'o');
xlabel('Training Set r');
ylabel('Test Set r');
title(['Training/Test Performance for ' titlename]);
legend(tokens, 'Location', 'SouthEast');
end