function ax = plot_perf_vs_truncation(data, modelnames, metric_name)
% A plot of model performance vs the amount of data that has been truncated
        
if ~any(strcmp(metric_name, {'r_ceiling', 'r_test', 'r_fit', 'r_test - r_fit'}))
    error('I only understand correlation metrics right now!');
end
    
% Compute the means and standard error of each model
d_nparm  = nanmean(data(:,:,1));  % Means should result in no change
D = data(:,:);
D(D==0) = NaN;
d_metric = abs(D);

d_means = nanmean(d_metric);
d_count = sum(~isnan(d_metric));
d_stddev = sqrt(nanvar(d_metric));
d_stderr = d_stddev ./ sqrt(d_count);
len = length(modelnames);

%---------------
% Create a structure to hold model names differing only by "trunc"
% Start by tokenizing and removing any "trunc" keywords
names = {};
for ii = 1:length(modelnames)
    tmp = tokenize_string(modelnames{ii});
    tmp = [tmp{:}];
    tmp = tmp( cellfun(@isempty, (strfind(tmp, 'trunc'))));
    mm = cellfun(@(x) [x '_'], tmp, 'UniformOutput', false);
    names{ii} = [mm{:}];
    names{ii} = names{ii}(1:end-1); % Trim off last _, if any    
end

% Create list of unique names
unames = unique(names);

scores = {};
for ii = 1:length(unames)    
    % Collect indexes that match each unique name
    idxs = strcmp(unames{ii}, names);      
    
    scores{ii} = [0, 0];
    for jj = 1:length(idxs)       
        idx = idxs(jj);
        if ~idx  % Skip false indexes
            continue;
        end
        % Get the truncation amount
        trunc = regexp(modelnames{jj}, 'trunc(\d+)', 'tokens');
        if isempty(trunc)
            trunc = 100;
        else
            trunc = str2double(trunc{1});
        end
        
        % Store it in a structure
        if ~isnan(d_means(jj))
            scores{ii} = cat(1, scores{ii}, [trunc d_means(jj)]);
        end
    end
    
    scores{ii} = sortrows(scores{ii});
    hold on;
    plot(scores{ii}(:,1), scores{ii}(:,2), 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerSize', 2, 'LineWidth', 2, 'color', pickcolor(ii), 'LineStyle', pickline(1+floor(ii/12)))
    hold off;
end

legend(unames, 'Location', 'SouthEast', 'Interpreter', 'none');
title('Model Performance vs Data Fraction');
xlabel('Percent of Estimation Set Data Used');
ylabel(sprintf('%s', metric_name), 'interpreter', 'none');

end