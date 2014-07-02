function score_2 = ga_computeDiversity( population, varargin )
% This function computes the genotypic L2 norm of the population

if numel(varargin) == 0
    n_sol = size(population,1);
    dists = squareform(pdist(population));
    score_2 = -min(dists+max(max(dists))*eye(n_sol));
else
    newpop = varargin{1};
    score_2 = [];
    for i=1:size(population,1)
        dists = pdist([population(i,:); newpop]);
        if isequal(size(dists),[1 1])
            score_2 = [score_2 -dists];
        else
            score_2 = [score_2 -min(dists(1,1:size(newpop,1)))];
        end
    end
end

end

