function fitgen00()

global GA_XXXHistory GA_XXXPointer GA_MaskFittable;

semse();

PopSize = 100;
TolFun = 1e-6; % 1e-6 is not enough?
Gen = 10000;
StallGen = 200;

    function xoverKids  = crossoverscatteredModified(parents,options,GenomeLength,FitnessFcn,unused,thisPopulation)
        
        %CROSSOVERSCATTERED Position independent crossover function.
        %   XOVERKIDS = CROSSOVERSCATTERED(PARENTS,OPTIONS,GENOMELENGTH, ...
        %   FITNESSFCN,SCORES,THISPOPULATION) creates the crossover children XOVERKIDS
        %   of the given population THISPOPULATION using the available PARENTS.
        %   In single or double point crossover, genomes that are near each other tend
        %   to survive together, whereas genomes that are far apart tend to be
        %   separated. The technique used here eliminates that effect. Each gene has an
        %   equal chance of coming from either parent. This is sometimes called uniform
        %   or random crossover.
        %
        %   Example:
        %    Create an options structure using CROSSOVERSCATTERED as the crossover
        %    function
        %     options = gaoptimset('CrossoverFcn' ,@crossoverscattered);
        
        %   Copyright 2003-2007 The MathWorks, Inc.
        %   $Revision: 1.9.4.4 $  $Date: 2007/08/03 21:23:15 $
        
        
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
            % Make sure that offspring are feasible w.r.t. linear constraints
            if constr
                feasible  = isTrialFeasible(xoverKids(i,:)',linCon.Aineq,linCon.bineq,linCon.Aeq, ...
                    linCon.beq,linCon.lb,linCon.ub,sqrt(options.TolCon));
                if ~feasible % Kid is not feasible
                    % Children are arithmetic mean of two parents (feasible w.r.t
                    % linear constraints)
                    alpha = rand;
                    xoverKids(i,:) = alpha*thisPopulation(r1,:) + ...
                        (1-alpha)*thisPopulation(r2,:);
                end
            end
        end
    end

    function mutationChildren = mutationgaussianModified(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation)
        %MUTATIONGAUSSIAN Gaussian mutation.
        
        r = range(thisPopulation);
        r(r==0) = 1;
        
        mutationChildren = zeros(length(parents),GenomeLength);
        for i=1:length(parents)
            parent = thisPopulation(parents(i),:);
            mask = 0;
            while (sum(mask) == 0)
                mask = floor(randi(length(parent),1,length(parent))/length(parent));
            end
            % half of the gaussian mutations are created with a scale
            % depending of the current population range, the other have 1
            if rand>0.5
                scale = r;
            else
                scale = options.StallGenLimit / options.genSinceLastChange;
            end
            mutationChildren(i,:) = parent  + scale .* randn(1,length(parent)) .* mask ;
        end
    end




    function [termcond, term_score, n_iters, term_phi] = two_dimensional_fitter_loop(fittername, highlevel_fn)
        
        global STACK META XXX;
        
        
        phi_init = pack_fittables(STACK);
        
        if isempty(phi_init)
            fprintf('Skipping because there are no parameters to fit.\n');
            term_cond = NaN;
            term_score = NaN;
            n_iters = 0;
            return
        end
        
        % initialization of the MaskFittable
        mask_rows = 0;
        for ii = 1:length(STACK)
            if isfield(STACK{ii}{1}, 'fit_fields')
                mask_rows = mask_rows + 1;
            end
        end
        
        GA_MaskFittable = zeros(mask_rows, length(phi_init));
        mask_rowi = 1;
        variable_index = 1;
        for ii = 1:length(STACK)
            if isfield(STACK{ii}{1}, 'fit_fields')
                l = 0;
                for i =1:length(STACK{ii}{1}.fit_fields)
                    l = l + size(STACK{ii}{1}.(cell2mat(STACK{ii}{1}.fit_fields(i))),1) * ...
                        size(STACK{ii}{1}.(cell2mat(STACK{ii}{1}.fit_fields(i))),2);
                end
                GA_MaskFittable(mask_rowi,variable_index:(variable_index+l-1)) = ii;
                variable_index = variable_index+l;
                mask_rowi = mask_rowi + 1;
            end
        end

        
        % initialization of the StimHistory and StackHistory
        GA_PhiHistory = Inf((size(GA_MaskFittable,1)-1)*PopSize*2, length(phi_init));
        GA_XXXHistory = cell((size(GA_MaskFittable,1)-1)*PopSize*2,1);
        GA_XXXPointer = 1;
        
        n_iters = 0;
        start_depth = find_fit_start_depth(STACK);
        fprintf('The STACK computation start depth is %d\n', start_depth);
        
        function score_diversity = my_obj_fn(phi2)
            
            n_sol = size(phi2,1);
            n_params = size(phi2,2);
            score_diversity = Inf(n_sol,2);
            
            for phi_i=1:n_sol
                
                unpack_fittables(phi2(phi_i,:));
                
                % check for already computed stimulus
                s = start_depth-1;
                for c=(size(GA_MaskFittable,1)-1):-1:1
                    [~, last_param_position] = find(GA_MaskFittable(c+1,:),1);
                    cols = 1:(last_param_position - 1);
                    phi = Inf(1,n_params);
                    phi(cols) = phi2(phi_i,cols);
                    [is_stored, stored_index] = ismember(phi,GA_PhiHistory,'rows');
                    if is_stored
                        s = GA_MaskFittable(c,last_param_position - 1);
                        XXX{s}.dat = GA_XXXHistory{stored_index};
                        break
                    end
                end
                
                calc_xxx(s+1);
                
                % try to update the cached copies
                for c=(size(GA_MaskFittable,1)-1):-1:1
                    [~, last_param_position] = find(GA_MaskFittable(c+1,:),1);
                    cols = 1:(last_param_position - 1);
                    phi = Inf(1,n_params);
                    phi(cols) = phi2(phi_i,cols);
                    if ~ismember(phi,GA_PhiHistory,'rows');
                        s = GA_MaskFittable(c,last_param_position - 1);
                        GA_XXXHistory{GA_XXXPointer} = XXX{s}.dat;
                        GA_PhiHistory(GA_XXXPointer,:) = phi;
                        GA_XXXPointer = GA_XXXPointer + 1;
                        if GA_XXXPointer > length(GA_XXXHistory)
                            GA_XXXPointer = 1;
                        end
%                     else
%                         break
                    end
                end
                
                [m, p] = META.perf_metric();
                score_diversity(phi_i,1) = m + p;
                
            end
            
            if (n_sol>1)
                dists = squareform(pdist(phi2));
                score_diversity(:,2) = -min(dists+max(max(dists))*eye(n_sol));
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


[termcond, term_score, n_iters, term_phi] = two_dimensional_fitter_loop('gagamultiobj()', ...
    @(obj_fn, phi_init) ga_gamultiobj(obj_fn, length(phi_init), [], [], [], [], [], [], ...
    gaoptimset('TolFun', TolFun, ...
    'PopulationSize', PopSize, ...
    'Generations', Gen, ...
    'Vectorized', 'on', ...
    'StallGenLimit', StallGen, ...
    'ParetoFraction', 1, ...
    'MutationFcn', @mutationgaussianModified, ...
    'CrossoverFraction', 0.2, ...%     'SelectionFcn', @selectiontournament, ...
    'CrossoverFcn', @crossoverscatteredModified, ...
    'DistanceMeasureFcn', @ga_distancecrowding, ...
    'InitialPopulation', repmat(phi_init,1,PopSize)')));

% TODO: Write a custom termination function to stop when the improvement on
% the best item alone for objective 1 is not big enough

unpack_fittables(term_phi);

    function [a,b,c,d] = step_until_10neg3(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3, ...
                'StopAtStepSize', 10^-6,...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end

    function [a,b,c,d] = step_until_10neg35(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3.5, ...
                'StopAtStepSize', 10^-6,...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end

    function [a,b,c,d] = step_until_10neg4(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-4, ...
                'StopAtStepSize', 10^-7,...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end

    function [a,b,c,d] = step_until_10neg45(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-4.5, ...
                'StopAtStepSize', 10^-7,...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end

    function [a,b,c,d] = step_until_10neg5(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-5, ...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end

    function [a,b,c,d] = step_until_10neg55(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-5.5, ...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end


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


% fit_iteratively(@step_until_10neg3, ...
%     create_term_fn());
% 
% fit_iteratively(@step_until_10neg35, ...
%     create_term_fn());
% 
% fit_iteratively(@step_until_10neg4, ...
%     create_term_fn());
% 
% fit_iteratively(@step_until_10neg45, ...
%     create_term_fn());
% 
% fit_iteratively(@step_until_10neg5, ...
%     create_term_fn());
% 
% fit_iteratively(@step_until_10neg55, ...
%     create_term_fn());
% 
% [~, term_score, ~] = fit_iteratively(@step_until_10neg6, ...
%     create_term_fn());

end