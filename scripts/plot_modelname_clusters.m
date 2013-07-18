function plot_modelname_clusters(batch, cellids, modelnames)
% For all modelnames
%  1. Calculate the avg loss to best-model performance if it were dropped    
%  2. Eliminate the least-impactful model by this metric
%  3. When a given number of categories has been reached, stop. 

N = 5; % Number of clusters to find

if N >= length(modelnames)
    error('There is no point in trying to cluster when comparing less than five modelnames.');
end

function relative_scores = scores_for_cellid(bat, cellid, models)

	sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(bat) ''];
	sql = [sql ' AND cellid="' cellid '"'];
    sql = [sql ' AND modelname in (' interleave_commas(models) ')'];
    results = mysql(sql);  
    max_score = max([results.r_test]);
    
    relative_scores = zeros(1, length(models));
    
    for ii = 1:length(models)
        m = models{ii};
        
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(bat) ''];
        sql = [sql ' AND cellid="' cellid '" AND modelname="' m '"'];
        results = mysql(sql);
        if isempty(results)
            fprintf('Skipping missing model...\n');
            continue;
        elseif length(results) > 1
            error('Too many results found!!? WTF!?\n');
        else
            relative_scores(ii) = results.r_test / max_score;
            if relative_scores(ii) < 0 || max_score < 0 
                relative_scores(ii) = NaN;
            elseif relative_scores(ii) > 1
                error('WTF? how did that get large?');
            end
        end       
    end   
end

% Create a scores cache
scores = zeros(length(cellids), length(modelnames));
for ci = 1:length(cellids)
    fprintf('Caching scores for cellid: %s [%d/%d]\n', cellids{ci}, ci, length(cellids)); 
    ret = scores_for_cellid(batch, cellids{ci}, modelnames);
    scores(ci, :) = ret(:);
end

% Create an objective function for optimizing modelname choices
function perf = obj_fn(indexes)    
    top_scores = max(scores(:,indexes), [], 2); % Find best score for each cellid
    perf = mean(top_scores); % TODO: Also try simple sum instead of mean here
end
    
% Start with all the first ten models, and try to 'boost'
chosen_models = 1:N;
N_max_steps = 1000;
for step = 1:N_max_steps;
    fprintf('Step %d of at most %d\n', step, N_max_steps); 
    
    best_score = obj_fn(chosen_models);
    best_pos = NaN;
    best_mdl_idx = NaN;
    for position = 1:N
        for mdl_idx = 1:length(modelnames)
            cache = chosen_models;
            cache(position) = mdl_idx;
            del = obj_fn(cache);
            if del > best_score
                best_score = del;
                best_pos = position;
                best_mdl_idx = mdl_idx;
            end
        end
    end
          
    % Stop iterating when the score can improve no more
    if isnan(best_pos)
        break;
    else % Otherwise take that one best step
        chosen_models(best_pos) = best_mdl_idx;
    end
end

thelabels = modelnames(chosen_models);

% Extract just the count of how many models fall into each
scrap = zeros(length(cellids), N);
chosen_scores = nan(size(scores));
for ci = 1:length(cellids)
    [~, idx] = max(scores(ci, chosen_models), [], 2);
    scrap(ci, idx) = 1;        
    chosen_scores(ci, idx) = scores(ci, chosen_models(idx));
    thelabels{idx} = [thelabels{idx} sprintf('\n%s (%5.3f)', cellids{ci}, chosen_scores(ci, idx))];    
end



for li = 1:N
    val = nanmean(chosen_scores(:, li));
    thelabels{li} = [thelabels{li} sprintf('\nMean Fraction of Maximum: %6.4f', val)];
end


thecount = sum(scrap, 1);
figure('Name', 'Model Clustering', 'NumberTitle', 'off');

h = pie(thecount, thelabels);
for ni = 1:N
    set(h(2*ni), 'Interpreter', 'none');
end


end