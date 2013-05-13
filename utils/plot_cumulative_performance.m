function fh = plot_cumulative_performance (batch, cellid, holdtokens, freetokens, ordinate)
% plot_cumulative_performance (batch, cellid, holdtokens, freetokens, ordinate)
%
% Plots cumulative performance curves for all models found whose batch and
% cellid matches the arguments to this function. If cellid is an empty
% string, all cellids will match it. 
%
% The abscissa (x) will be the range of values taken by the ordinate. The
% ordinate (y) will be the fraction of the models with values larger than
% the value on the abcissa at that point. 
%
% ARGUMENTS:
%    batch       Number of the batch
%    cellid      Cellid
%    holdtokens  Which tokens to use common across all model queries. 
%    freetokens  Which tokens to use one-at-a-time for model queries. 
%    ordinate    Which value to plot along the ordinate axis. 
%
% RETURNS: 
%    fh          A newly created figure handle showing the cumulative
%                performance of the models.

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
    
    % This is just so the number is listed in the legend:
    freetokens{token_idx} = [freetokens{token_idx} '(' num2str(length(models)) ')'];
end

fh = figure;

set(gca, 'ColorOrder', [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0;], ...
         'LineStyleOrder',{'-','--',':'},'NextPlot','ReplaceChildren');

plot(x, y);
xlabel(ordinate, 'Interpreter', 'none');
ylabel(['Cumulative fraction of models with ' ordinate ' > value'], 'Interpreter', 'none')
legend(freetokens, 'Interpreter', 'none');
title(['Performance of Tokens (# Models: ' num2str(mdlcount) ') \{' write_readably(holdtokens) '\}']);
