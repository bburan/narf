function xoverKids  = ga_crossover(parents,options,GenomeLength,FitnessFcn,unused,thisPopulation)

global GA_MaskFittable;

% crossover procedure
% the cross-over is performed on the MODULE subdivision with random weight
% average, meaning that two modules with parameters [x1 x2 ... xn] and 
% [y1 y2 ... yn] will produce with r a random number in [0,1]:
% [x1 x2 ... xn] * r + [y1 y2 ... yn] * (1-r)

% How many children to produce?
nKids = length(parents)/2;
% Extract information about linear constraints, if any
linCon = options.LinearConstr;
constr = ~isequal(linCon.type,'unconstrained');
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);

% To move through the parents twice as fast as thekids are
% being produced, a separate index for the parents is needed
index = 1;
% for each kid...
for i=1:nKids
    % get parents
    r1 = parents(index);
    index = index + 1;
    r2 = parents(index);
    index = index + 1;
    % Randomly select half of the genes from each parent
    % This loop may seem like brute force, but it is twice as fast as the
    % vectorized version, because it does no allocation.
    for j = 1:size(GA_MaskFittable,1)
        % this has 1/2 prob to take either one of the parents
        %                 if(rand > 0.5)
        %                     xoverKids(i,logical(GA_MaskFittable(j,:))) = thisPopulation(r1,logical(GA_MaskFittable(j,:)));
        %                 else
        %                     xoverKids(i,logical(GA_MaskFittable(j,:))) = thisPopulation(r2,logical(GA_MaskFittable(j,:)));
        %                 end
        % this does a randomly weighted average
        alpha = rand;
        xoverKids(i,logical(GA_MaskFittable(j,:))) = alpha * thisPopulation(r1,logical(GA_MaskFittable(j,:))) + ...
            (1 - alpha) * thisPopulation(r2,logical(GA_MaskFittable(j,:)));
    end
%     % Make sure that offspring are feasible w.r.t. linear constraints
%     if constr
%         feasible  = isTrialFeasible(xoverKids(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
%             linCon.beq,linCon.lb,linCon.ub,sqrt(options.TolCon));
%         if ~feasible % Kid is not feasible
%             % Children are arithmetic mean of two parents (feasible w.r.t
%             % linear constraints)
%             alpha = rand;
%             xoverKids(i,:) = alpha*thisPopulation(r1,:) + ...
%                 (1-alpha)*thisPopulation(r2,:);
%         end
%     end
end
end