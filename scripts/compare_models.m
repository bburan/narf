function compare_models(batch, cellids, modelnames)
% Computes several interesting graphs for a combination of models

    function ret = eval_perfmetric(x)
        global META;
        calc_xxx(1);
        ret = META.perf_metric();        
    end

stack_extractor = @(x) 0;
meta_extractor = @eval_perfmetric;

[~, score, ~] = load_model_batch(batch, cellids, modelnames, ...
                       stack_extractor, meta_extractor);

%times(cellfun(@iscell, score)) = {nan};
score = cell2mat(score)';

score(score>2) = 2;

figure('Name', 'Model Comparison', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
n = size(score,2);
for ii = 1:n
    subplot(n,1,ii);
    hist(score(:,ii), 50);   
    ll = xlabel(sprintf('NMSE for %s', modelnames{ii}));
    set(ll, 'Interpreter', 'none');
end

figure('Name', 'Model Comparison', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
n = size(score,2);
for ii = 1:n
    subplot(n,1,ii);
    hist(log(score(:,ii)), 50);   
    ll = xlabel(sprintf('log(NMSE) for %s', modelnames{ii}));
    set(ll, 'Interpreter', 'none');
end

figure('Name', 'Model Comparison', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
n = size(score,2);
for ii = 1:n
    subplot(n,1,ii);
    hist(1./(score(:,ii)), 50);   
    ll = xlabel(sprintf('1/(NMSE) for %s', modelnames{ii}));
    set(ll, 'Interpreter', 'none');
end

end

