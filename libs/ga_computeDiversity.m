function score_2 = ga_computeDiversity( population, varargin )
% This function computes the genotypic L2 norm of the population

% hack to exclude the last 5 parameters of the distance (because they span such a long space without much change)
selected_parameters = 1:(size(population,2)-5);

if numel(varargin) == 0
    n_sol = size(population,1);
    dists = squareform(pdist(population(:,selected_parameters)));
    score_2 = -min(dists+max(max(dists))*eye(n_sol));
else
    newpop = varargin{1};
    score_2 = [];
    for i=1:size(population,1)
        dists = pdist([population(i,selected_parameters); newpop(:,selected_parameters)]);
        if isequal(size(dists),[1 1])
            score_2 = [score_2 -dists];
        else
            score_2 = [score_2 -min(dists(1,1:size(newpop,1)))];
        end
    end
end

end

