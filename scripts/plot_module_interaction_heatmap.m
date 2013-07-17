function plot_module_interaction_heatmap(batch, cellids, modelnames)
% Plots a 2 by 2 heatmap using the mean validation set performance.
% Expects that all modelnames only differ by TWO tokens at most
% Plots the NANMEAN across all cellids

tokens = {};

% Reconstruct which model keywords are different across all models
for ii = 1:length(modelnames)
    mn = modelnames{ii};
    toks = tokenize_string(mn);
    % toks = unique([toks{:}]); 
    for jj = 1:length(toks)        
        tokens{jj}{ii} = toks{jj};        
    end
end

for ii = 1:length(tokens)
    tokens{ii} = unique([tokens{ii}{:}]);
end
    
% If there are more or less than two varying keywords, we can't do a heatmap. 
lens = cellfun(@length, tokens);
if sum(lens ~= 1) ~= 2
    disp(tokens{:});
    error('Cannot plot a 2x2 heatmap unless exactly 2 tokens vary simultaneously.');
end

idxs = 1:length(lens);
idxs = idxs(lens~=1);

abscissa = tokens{idxs(1)};
ordinate = tokens{idxs(2)};

data = size(length(ordinate), length(abscissa));

default_modelname = {};
for ii = 1:length(tokens)
    default_modelname{ii} = tokens{ii}{1};
end

dbopen; 
for ii = 1:length(ordinate)
    for jj = 1:length(abscissa)
        mn = default_modelname;
        mn{idxs(2)} = ordinate{ii};
        mn{idxs(1)} = abscissa{jj};

        tmp = cellfun(@(n) sprintf('%s_', n), mn, 'UniformOutput', false);
        modelname = strcat(tmp{:});
        modelname = modelname(1:end-1);
        
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(batch) ''];
        sql = [sql ' AND cellid in (' interleave_commas(cellids) ')'];
        sql = [sql ' AND modelname="' modelname '"'];              
        results = mysql(sql);  

        % NOTE: Ask Ivar if you want a different metric and he'll refactor
        % this so we can extract an arbitrary thing
        data(ii,jj) = nanmean([results.r_test]);
    end
end

figure;
heatmap(data, abscissa, ordinate, '',  'TickAngle', 90, 'ShowAllTicks', true);
title(sprintf('Mean Validation Set Correlation for Each Keyword Interactions (%d cells)', length(cellids)));
