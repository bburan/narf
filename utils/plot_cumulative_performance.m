function plot_cumulative_performance (batch, cellid, holdtokens, freetokens, ordinate)
% plot_cumulative_performance (batch, cellid, holdtokens, freetokens, ordinate)
%
% Plots cumulative performance curves for all 
%
% ARGUMENTS:
%    batch       Number of the batch
%    cellid      Cellid
%    holdtokens  
%    freetokens  
%    ordinate    
%
% RETURNS: Nothing

n_pts = 200;
token_count = length(freetokens);
x = zeros(n_pts, token_count);
y = zeros(n_pts, token_count);
mdlcount = 0;

for token_idx = 1:length(freetokens)
    
    models = db_get_models(batch, cellid, cat(2, holdtokens, freetokens{token_idx}));
    mdlcount = mdlcount + length(models);
    
    vals = [models(:).(ordinate)];
    
    if isempty(vals)
        fprintf('No matching models found. Unable to plot.\n');
        return
    end
    
    x(:,token_idx) = linspace(min(vals(:)), max(vals(:)), n_pts);
    
    n = length(models);
    
    for ii = 1:n_pts;
        y(ii, token_idx) = sum(vals < x(ii, token_idx) ) / n;
    end
end

figure;

set(gca, 'ColorOrder', [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0;], ...
         'LineStyleOrder',{'-','--',':'},'NextPlot','ReplaceChildren');

plot(x, y);
xlabel(ordinate, 'Interpreter', 'none');
ylabel(['Cumulative fraction of models with ' ordinate ' > value'], 'Interpreter', 'none')
legend(freetokens, 'Interpreter', 'none');
title(['Performance of Tokens (# Models: ' num2str(mdlcount) ')']);
