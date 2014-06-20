function fitgen00()

semse();

PopSize = 100;
TolFun = 1e-6;
Gen = 150; % should replace by 1000



    function crowdingDistance = distancecrowdingModified(pop,score,options)
        % THIS VERSION IS ADAPTED TO THE DIVERSIFICATION OBJECTIVE BY
        % RETAINING ALWAYS ALL THE BEST FRONT MEMBERS TO BE THE BEST FITTING
        % SOLUTION FOR OBJECTIVE 1.
        
        %DISTANCECROWDING Assign local crowding distance to each individual
        %   CROWDINGDISTANCE = DISTANCECROWDING(POP,SCORE,OPTIONS,SPACE) Calculates
        %   crowding distance for each individuals on a non-dominated front. The
        %   fourth argument SPACE can be 'phenotype' or 'genotype' for distance to
        %   be in function space or decision variable space respectively.
        %
        %   Example:
        %   Create an options structure using DISTANCECROWDING as the distance
        %   function in decision variable space
        %     options = gaoptimset('DistanceMeasureFcn',{@distancecrowding,'genotype'});
        
        %   Reference: Kalyanmoy Deb, "Multi-Objective Optimization using
        %   Evolutionary Algorithms", John Wiley & Sons ISBN 047187339, pg: 245 -
        %   253
        
        %   Copyright 2007 The MathWorks, Inc.
        %   $Revision: 1.1.6.1 $  $Date: 2007/05/23 18:49:35 $
        
        
        [~, sorted_distance] = sort(score(:,1));
        crowdingDistance = 1/sorted_distance;
        crowdingDistance(sorted_distance(1)) = Inf;
        
%         popSize = size(y,1);
%         numData = size(y,2);
%         crowdingDistance = zeros(popSize,1);
%         
%         for m = 1:numData
%             data = y(:,m);
%             % Normalize obective before computing distance
%             data = data./(1 + max(abs(data(isfinite(data)))));
%             [sorteddata,index] = sort(data);
%             % The best and worst individuals are at the end of Pareto front and
%             % they are assigned Inf distance measure
%             crowdingDistance([index(1),index(end)]) = Inf;
%             % Distance measure of remaining individuals
%             i = 2;
%             while i < popSize
%                 crowdingDistance(index(i)) = crowdingDistance(index(i)) + ...
%                     min(Inf, (data(index(i+1)) - data(index(i-1))));
%                 i = i+1;
%             end
%         end
    end




    function [termcond, term_score, n_iters] = two_dimensional_fitter_loop(fittername, highlevel_fn)
        
        global STACK META;
        
        phi_init = pack_fittables(STACK);
        
        if isempty(phi_init)
            fprintf('Skipping because there are no parameters to fit.\n');
            term_cond = NaN;
            term_score = NaN;
            n_iters = 0;
            return
        end
        
        n_iters = 0;
        start_depth = find_fit_start_depth(STACK);
        
        function score_diversity = my_obj_fn(phi2)
            n_sol = size(phi2,1);
            score_diversity = Inf(n_sol,2);
            
            for phi_i=1:n_sol
                
                unpack_fittables(phi2(phi_i,:));
                calc_xxx(start_depth);
                
                [m, p] = META.perf_metric();
                score_diversity(phi_i,1) = m + p;
                
            end
            
            if (n_sol>1)
                dists = squareform(pdist(phi2));
                score_diversity(:,2) = -min(dists+max(max(dists))*eye(n_sol));
            end
            
            % Print a newline and important info every 100 iterations
            if isequal(mod(n_iters, 100), 1)
                fprintf('\n[Generation: %5d, score:%f...%f...%f]', round(n_iters/PopSize), min(score_diversity(:,1)), mean(score_diversity(:,1)), max(score_diversity(:,1)));
                % dbtickqueue(n_iters);
            end
            
            % Print 1 progress dot for every 5 iterations no matter what
            if isequal(mod(n_iters, 5), 0) %% && ~bequiet
                fprintf('.');
            end
            
            n_iters = n_iters + n_sol;
        end
        
        fprintf('Fitting %d variables with %s', length(phi_init), fittername);
        
        [X, FVAL, termcond] = highlevel_fn(@my_obj_fn, phi_init);
        
        [term_score, minvali] = min(FVAL(:,1));
        term_phi = X(minvali,:);
        
        unpack_fittables(term_phi);
        calc_xxx(start_depth);
        
        fprintf('Complete fit with %d objective function evaluations.\n', n_iters);
        fprintf('----------------------------------------------------------------------\n');
        
    end



[termcond, term_score, n_iters] = two_dimensional_fitter_loop('gamultiobj()', ...
    @(obj_fn, phi_init) gamultiobj(obj_fn, length(phi_init), [], [], [], [], [], [], ...
    gaoptimset('TolFun', TolFun, ...
    'PopulationSize', PopSize, ...
    'Generations', Gen, ...
    'Vectorized', 'on', ...
    'StallGenLimit', 10, ...
    'ParetoFraction', 0.5, ...
    'CrossoverFraction', 0.5, ...
    'DistanceMeasureFcn', {@distancecrowdingModified}, ...
    'InitialPopulation', repmat(phi_init,1,PopSize)' )));




function [a,b,c,d] = step_until_10neg6(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-6, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

% NOW DO THE LAST FIT05C STEP:

MaxStepsPerIteration=10;
StepGrowth=1.1;



fit_iteratively(@step_until_10neg6, ...
                create_term_fn());

end