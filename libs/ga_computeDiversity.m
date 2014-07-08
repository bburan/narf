function score_2 = ga_computeDiversity( population, varargin )
% This function computes the genotypic L2 norm of the population

% hack to exclude the last 5 parameters of the distance (because they span such a long space without much change)
selected_parameters = 1:(size(population,2)-5);

    function minvec = getminvec(pop)
        n_sol = size(pop,1);
        dists = squareform(pdist(pop));
        minvec = -min(dists+max(max(dists))*eye(n_sol));
    end

if numel(varargin) == 0
    score_2 = getminvec(population(:,selected_parameters));
else
    newpop = varargin{1};
    dists = dist([population(:,selected_parameters);  newpop(:,selected_parameters)]');
    score_2 = -min(dists(1:size(population,1),(end+1-size(newpop,1)):end),[],2)';
%     score_2 = zeros(1,size(population,1));
%     for i=1:size(population,1)
%         dists = pdist([population(i,selected_parameters); newpop(:,selected_parameters)]);
%         if isequal(size(dists),[1 1])
%             score_2(i) = -dists;
%         else
%             score_2(i) = -min(dists(1,1:size(newpop,1)));
%         end
%     end
end

end

