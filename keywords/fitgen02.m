function fitgen02()

global GA_PhiHistory GA_XXXHistory GA_XXXPointer GA_MaskFittable GA_LowerBound GA_UpperBound XXX META STACK;

semse();

PopSize = 100;
TolFun = 1e-6; % 1e-6 is not enough?
Gen = 10000;
StallGen = 200;


    function [termcond, term_score, n_iters, term_phi] = two_dimensional_fitter_loop(fittername, highlevel_fn)
                
        phi_init = pack_fittables(STACK);
        
        if isempty(phi_init)
            fprintf('Skipping because there are no parameters to fit.\n');
            term_cond = NaN;
            term_score = NaN;
            n_iters = 0;
            return
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

phi_init = pack_fittables(STACK);
initialize_GA_GLOBALS(phi_init);

phi0 = pack_fittables(STACK);

[termcond, term_score, n_iters, term_phi] = two_dimensional_fitter_loop('gagamultiobj()', ...
    @(obj_fn, phi_init) ga_gamultiobj(obj_fn, length(phi_init), [], [], [], [], [], [], ...
    gaoptimset('TolFun', TolFun, ...
    'PopulationSize', PopSize, ...
    'Generations', Gen, ...
    'Vectorized', 'on', ...
    'StallGenLimit', StallGen, ...
    'ParetoFraction', 1, ...
    'MutationFcn', @ga_mutation, ...
    'CrossoverFraction', 0.2, ...%     'SelectionFcn', @selectiontournament, ...
    'CrossoverFcn', @ga_crossover, ...
    'DistanceMeasureFcn', @ga_distancecrowding, ...
    'InitialPopulation', [phi0'; rand(PopSize-1,length(phi_init)).*repmat(GA_UpperBound-GA_LowerBound, PopSize-1,1) + repmat(GA_LowerBound, PopSize-1,1)])));

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

% NOW DO ALL THE FIT05C STEPS:

MaxStepsPerIteration=10;
StepGrowth=1.1;


fit_iteratively(@step_until_10neg3, ...
    create_term_fn());

fit_iteratively(@step_until_10neg35, ...
    create_term_fn());

fit_iteratively(@step_until_10neg4, ...
    create_term_fn());

fit_iteratively(@step_until_10neg45, ...
    create_term_fn());

fit_iteratively(@step_until_10neg5, ...
    create_term_fn());

fit_iteratively(@step_until_10neg55, ...
    create_term_fn());

[~, term_score, ~] = fit_iteratively(@step_until_10neg6, ...
    create_term_fn());

end