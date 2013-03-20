function plot_cumulative_performance (batch, cellid, tokens, ordinate)

n_pts = 200;
token_count = length(tokens);
x = zeros(n_pts, token_count);
y = zeros(n_pts, token_count);

for token_idx = 1:length(tokens)
    
    models = db_get_models(batch, cellid, {tokens{token_idx}});
    
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
plot(x, y);
xlabel(ordinate, 'Interpreter', 'none');
ylabel(['Cumulative fraction of models with ' ordinate ' > value'], 'Interpreter', 'none')
legend(tokens, 'Interpreter', 'none');
title(['Performance of Tokens (# Models: ' num2str(numel(x)) ')']);