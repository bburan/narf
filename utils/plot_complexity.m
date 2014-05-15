function ax = plot_complexity(data, modelnames, metric_name)
% A plot of model complexity (number of parameters) vs FOV remaining
        
if ~any(strcmp(metric_name, {'r_ceiling', 'r_test', 'r_fit'}))
    error('I only understand correlation metrics right now!');
end
    
% Compute the means and standard error of each model
d_nparm  = nanmean(data(:,:,1));  % Means should result in no change
D = data(:,:,2);
D(D==0) = NaN;
d_metric = abs(D);
%d_fittime = nanmean(data(:,:,3));

d_means = nanmean(d_metric);
d_count = sum(~isnan(d_metric));
d_stddev = sqrt(nanvar(d_metric));
d_stderr = d_stddev ./ sqrt(d_count);
len = length(modelnames);

figure('Name', 'Complexity Plot', 'NumberTitle', 'off', 'Position', [10 10 1000 1000]);
ax = axes();

% UGLY HACK TO ADD NUMBER OF PARMS FOR INITIAL FILTERS
for ii = 1:len
    if regexp(modelnames{ii}, '^but', 'once')
        d_nparm(ii) = d_nparm(ii) + 3;
    end
    if regexp(modelnames{ii}, '^ell', 'once')
        d_nparm(ii) = d_nparm(ii) + 3;
    end
    if regexp(modelnames{ii}, '^gam', 'once')
        d_nparm(ii) = d_nparm(ii) + 2;
    end
    if regexp(modelnames{ii}, '^apgt\d\dx\d', 'once')
        d_nparm(ii) = d_nparm(ii) + 2;
    end        
    if regexp(modelnames{ii}, '^ozgf\d\dx\d', 'once')
        d_nparm(ii) = d_nparm(ii) + 3;
    end
    if regexp(modelnames{ii}, '^szgf\d\dx\d', 'once')
        d_nparm(ii) = d_nparm(ii) + 3;
    end
end
names = shorten_modelnames(modelnames);

% Scatter plot with text labels
hold on;       

jitter = randn(size(d_means));

for pass = 1:2
    for ii = 1:len
        name = names{ii};
        yc = 1 - d_means(ii);
        yt = yc + d_stderr(ii); 
        yb = yc - d_stderr(ii);
        x = d_nparm(ii) + 0.1*jitter(ii);                
        
        if pass == 1
            line([x+0.1 x-0.1 x x x-0.1 x+0.1], [yt, yt, yt, yb, yb, yb], 'Linewidth', 2, 'Linestyle', '-', 'Color', pickcolor(ii));
        else
            %fittime = d_fittime(ii);
            %text(x, yc, sprintf('%s [%0.1f]', name, fittime/60));
            text(x, yc, name);    
        end
    end
end

hold off
xlabel('Number of Parameters');
ylabel(sprintf('1 - %s', metric_name));
xticks(1:ceil(max(d_nparm(:))));

end