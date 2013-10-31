function cluster_population(batch, cellids, modelnames)
% For each cellid
%    1. Check all pairings of models (n^2 comparisons for n models)
%    2. Record which pairings are "significantly" better
%    3. Create a population PIE chart from those pairings.

name_statistic = 'r_test';
name_noise     = 'r_floor';

function [scores, noises] = scores_for_cellid(bat, cellid, models)
    scores = [];
    noises = [];
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
            scores(ii) = abs(results.(name_statistic));
            noises(ii) = abs(results.(name_noise));
        end       
    end
end

N_models = length(modelnames);
N_cells = length(cellids);
scores = zeros(N_cells, N_models);
noises = zeros(N_cells, N_models);

for ci = 1:N_cells
    fprintf('Computing scores for cellid: %s [%d/%d]\n', cellids{ci}, ci, N_cells); 
   	[ret1, ret2] = scores_for_cellid(batch, cellids{ci}, modelnames);
    scores(ci, :) = ret1(:);
    noises(ci, :) = ret2(:);
end

noises = mean(noises, 2); % Average noise floor across all models

modelnames = shorten_modelnames(modelnames);

% SCIENCE WARNING: This type of statistical testing is p-value fishing. 
%  By comparing many different pairings, we are much more likely to
%  find 'significantly better' things just by chance. We should really
%  account for the (N_models-1)^2 comparisons we are doing here, or use
%  a Bayesian method instead of significance testing. 

results_by_cellid = {};
for ci = 1:N_cells,
    results_by_cellid{ci} = {};
    for mi = 1:N_models,
        for oi = 1:N_models,
            if mi == oi
                continue;
            end
            % Gaussianized significance test
            if (scores(ci, mi) > scores(ci, oi) + noises(oi))
                if isempty(results_by_cellid{ci})
                    results_by_cellid{ci}{end+1} = [char(modelnames(mi)) '>' char(modelnames(oi))];
                else
                    results_by_cellid{ci}{end+1} = [', ' char(modelnames(mi)) '>' char(modelnames(oi))];
                end
            end
        end
    end
end

% Flatten the results into single strings
for ci = 1:N_cells,
    results_by_cellid{ci} = char([results_by_cellid{ci}{:}]);
end

% Make unique groupings
frequency = [];
thelabels = {};
uniqs = unique(results_by_cellid);
for li = 1:length(uniqs)
    boolshit = strcmp(uniqs{li}, results_by_cellid);
    count = sum(boolshit);
    frequency(li) = count / N_cells;
    if isempty(uniqs{li})
        thelabels{li} = '[NULL HYPOTHESIS]';
    else
        thelabels{li} = sprintf('[%s]', uniqs{li});
        thelabels{li} = [thelabels{li} sprintf('\n%s', cellids{boolshit})];
    end
end

% Make pie chart
figure('Name', 'Population Clusters', 'NumberTitle', 'off');
h = pie(frequency, thelabels);
for ni = 1:N
    set(h(2*ni), 'Interpreter', 'none');
end

end