function mutationChildren = ga_mutation(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation)

global GA_LowerBound GA_UpperBound GA_Bounded;

% mutation procedure
% for unbounded variables, new variables are drawn from a Laplacian
%  distribution, with a hybrid scale parameter which is
%  (50% prob) scaled by the previous generation range
%  (50% prob) scaled by 1 - (current stall time) / (max stall time)
% for bounded variables, new variables are drawn with polynomial
% mutation, with eta=20 and hybrid disruption:
% (30% prob) highly disruptive
% (70% prob) non disruptive

[~,r_idx] = sort(thisScore(:,1));
r = range(thisPopulation(r_idx(1:10),:));
r(r==0) = 1;

% mutationChildren = zeros(length(parents),GenomeLength);
mutationChildren = thisPopulation(parents,:);
for i=1:length(parents)
    parent = thisPopulation(parents(i),:);
    mask = 0;
    while (sum(mask) == 0)
        mask = floor(randi(length(parent),1,length(parent))/length(parent));
    end
    % half of the gaussian mutations are created with a scale
    % depending of the current population (good-scoring) range, the other have 1
    if rand>0.5
        scale = r;
    else
        scale = 1 - options.genSinceLastChange / options.StallGenLimit;
    end
    
    % let's draw random values from a Laplace distribution for
    % unbounded variables
    if (sum(~GA_Bounded))
        u = rand(1,sum(~GA_Bounded)) - 0.5;
        mutationChildren(i,~GA_Bounded) = parent(~GA_Bounded) + scale .* (sign(u).*log(1-2*abs(u))) .* mask(~GA_Bounded);
    end
    
    % let's use the polynomial mutation for bounded variables
    for j = find(GA_Bounded)
        if mask(j)
            eta = 20+1;
            if rand > 0.3
                % 70% of non-disruptive polynomial mutation (following Hamdan 2010)
                delta = min(GA_UpperBound(j) - parent(j), parent(j) - GA_LowerBound(j)) / ...
                    (GA_UpperBound(j) - GA_LowerBound(j));
                r = rand;
                if r<0.5
                    delta = (2*r + (1-2*r) * (1-delta)^eta)^(1/eta) - 1;
                else
                    delta = 1 - (2*(1-r) + 2*(r-0.5) * (1-delta)^eta)^(1/eta);
                end
            else
                % 30 % of highly-disruptive polynomial mutation
                delta1 = (parent(j) - GA_LowerBound(j)) / ...
                    (GA_UpperBound(j) - GA_LowerBound(j));
                delta2 = (GA_UpperBound(j) - parent(j)) / ...
                    (GA_UpperBound(j) - GA_LowerBound(j));
                r = rand;
                if r<0.5
                    delta = (2*r + (1-2*r) * (1-delta1)^eta)^(1/eta) - 1;
                else
                    delta = 1 - (2*(1-r) + 2*(r-0.5) * (1-delta2)^eta)^(1/eta);
                end
            end
            newval = parent(j) + delta * (GA_UpperBound(j) - GA_LowerBound(j));
%             mutationChildren(i,:) = parents(i,:);
            mutationChildren(i,j) = real(newval);
        end
    end
    %             mutationChildren(i,:) = parent  + scale .* randn(1,length(parent)) .* mask ;
end

end