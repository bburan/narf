function fitgen04()

global GA_PhiHistory GA_XXXHistory GA_XXXPointer GA_MaskFittable GA_LowerBound GA_UpperBound STACK META XXX;

semse();

% genetic algorithm parameters
PopSize = 100; % 100
TolFun = 1e-6; % 1e-6 is not enough?
Gen = 10000;
StallGen = 200;

% fit05c parameters
MaxStepsPerIteration=10;
StepGrowth=1.1;

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
                    end
                end
                
                [m, p] = META.perf_metric();
                score_diversity(phi_i,1) = m + p;
                
            end
            
            if (n_sol>1)
                dists = squareform(pdist(phi2));
                score_diversity(:,2) = -min(dists+max(max(dists))*eye(n_sol));
            end
            
            % Print 1 progress dot for every generation
            fprintf('.');
            
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


    function fn = make_subfitter(del)
        function [a,b,c,d] = subfitter(prev_opts)
            
            % Detect whether the fittables are in a FIR block or not
            module_being_fit = '';
            for kk = 1:length(STACK)
                if isfield(STACK{kk}{1}, 'fit_fields') && ~isempty(STACK{kk}{1}.fit_fields)
                    module_being_fit = STACK{kk}{1}.name;
                    break;
                end
            end
            
            if strcmp(module_being_fit, 'fir_filter') || ...
                    (strcmp(module_being_fit, 'weight_channels') && ...
                    isempty(STACK{kk}{1}.phifn))
                if exist('prev_opts', 'var')
                    [a,b,c,d] = fit_boo(prev_opts);
                else
                    [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', del, ...
                        'StopAtStepNumber', 1, ...
                        'StepGrowth', 1.3);
                end
            else
                if exist('prev_opts', 'var')
                    [a,b,c,d] = fit_scaat(prev_opts);
                else
                    [a,b,c,d] = fit_scaat('StopAtAbsScoreDelta', del, ...
                        'StopAtStepNumber', 10);
                end
            end
        end
        
        fn = @subfitter;
        
    end




phi0 = pack_fittables(STACK);


% initialization of a best solution as found with fit05c
fit05c();
phi_fit05c = pack_fittables(STACK);
calc_xxx(1);
[m, p] = META.perf_metric();
best_indiv = phi_fit05c;
best_indivs = best_indiv;
best_score = m+p;
best_origin = 'initial fit05c';
best_scores = m+p;
fprintf('\n\n ===== SOLUTION FROM FIT05C ===== SCORE is %f\n\n', best_score);




% Do a run of fit05g+fit05c like
unpack_fittables(phi0);calc_xxx(1); % reset
fit05g();
fit05c();
phi_fit05gfit05c = pack_fittables(STACK);
calc_xxx(1);

[m, p] = META.perf_metric();
if m+p < best_score
    best_indiv = phi_fit05gfit05c;
    best_score = m+p;
    best_origin = 'initial fit05g+fit05c -like';
    fprintf('\n\n ===== NEW BEST SOLUTION ===== SCORE is %f\n\n', best_score);
end
best_indivs = [best_indivs phi_fit05gfit05c];
best_scores = [best_scores m+p];
fprintf('\n\n ===== SOLUTION FROM FIT05g+FIT05c ===== SCORE is %f\n\n', m+p);





% Do a run of fit09-like
unpack_fittables(phi0);calc_xxx(1);
% Initialization: If FIR filters are all zero, initialize them randomly
[~, mod_idxs] = find_modules(STACK, 'fir_filter');
for ii = 1:length(mod_idxs)
    for jj = 1:length(STACK{mod_idxs{ii}})
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
           any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields)) && ...
           all(0 == STACK{mod_idxs{ii}}{jj}.coefs(:))
            STACK{mod_idxs{ii}}{jj}.coefs = normrnd(0, 10^-3, size(STACK{mod_idxs{ii}}{jj}.coefs));
            fprintf('=====> Initializing FIR coefs to random numbers!\n');
        end
    end
end
% Unpack the initial stack
% Now gradually shrink the stopping criterion
scale=10^-1;
stop_at=10^-2;%10^-6;
while(scale > stop_at)
    fit_iteratively(make_subfitter(scale), create_term_fn('StopAtAbsScoreDelta', scale));
    scale = scale * 0.666; % Very conservative: 0.8. Probably 0.5 or even 0.25 is fine.
end
phi_fit09 = pack_fittables(STACK);
calc_xxx(1);
[m, p] = META.perf_metric();
if m+p < best_score
    best_indiv = phi_fit09;
    best_score = m+p;
    best_origin = 'initial fit09-like';
    fprintf('\n\n ===== NEW BEST SOLUTION ===== SCORE is %f\n\n', best_score);
end
best_indivs = [best_indivs phi_fit09];
best_scores = [best_scores m+p];
fprintf('\n\n ===== SOLUTION FROM FIT09 ===== SCORE is %f\n\n', m+p);




% initialization of all the global variables used by genetic algorithms
initialize_GA_GLOBALS(phi0)

% now do the optimization 10 times
for run=1:10
    fprintf('\n\n ===== STARTING RUN %d\n\n', run);
    % Unpack the initial stack
    unpack_fittables(phi0);calc_xxx(1);
    
    [termcond, term_score, n_iters, term_phi] = two_dimensional_fitter_loop('gagamultiobj()', ...
        @(obj_fn, phi_init) ga_gamultiobj(obj_fn, length(phi_init), [], [], [], [], [], [], ...
        gaoptimset('TolFun', TolFun, ...
        'PopulationSize', PopSize, ...
        'Generations', Gen, ...
        'Vectorized', 'on', ...
        'StallGenLimit', StallGen, ...
        'ParetoFraction', 0.1, ...
        'MutationFcn', @ga_mutation, ...
        'CrossoverFraction', 0.2, ...%     'SelectionFcn', @selectiontournament, ...
        'CrossoverFcn', @ga_crossover, ...
        'DistanceMeasureFcn', @ga_distancecrowding, ...
        'InitialPopulation', [phi0'; rand(PopSize-1,length(phi_init)).*repmat(GA_UpperBound-GA_LowerBound, PopSize-1,1) + repmat(GA_LowerBound, PopSize-1,1)])));
    %     'InitialPopulation', repmat(phi_init(1:(length(phi_init))),1,PopSize)')));
    
    unpack_fittables(term_phi);
    
    % NOW DO THE LAST FIT05C STEPS:
    
    fit_iteratively(@step_until_10neg3, ...
        create_term_fn());
    %     fit_iteratively(@step_until_10neg5, ...
    %         create_term_fn());
    %     fit_iteratively(@step_until_10neg55, ...
    %         create_term_fn());
    %     [~, term_score, ~] = fit_iteratively(@step_until_10neg6, ...
    %         create_term_fn());
    
    % compute the last score - update if necessary
    phi_candidate = pack_fittables(STACK);
    calc_xxx(1);
    [m, p] = META.perf_metric();
    if best_score > m+p
        % we have a winner
        best_indiv = phi_candidate;
        best_score = m+p;
        best_origin = ['genetic fit number ' num2str(run)];
        fprintf('\n\n ===== NEW BEST SOLUTION ===== SCORE is %f\n\n', best_score);
    end
    
    best_indivs = [best_indivs phi_candidate];
    best_scores = [best_scores m+p];
end

best_indivs = best_indivs';

fprintf('\n\n ===== STARTING LAST RUN %d\n\n', run);

% now the last round: using all the previously found optimal solutions
% Unpack the initial stack
unpack_fittables(phi0);calc_xxx(1);
[termcond, term_score, n_iters, term_phi] = two_dimensional_fitter_loop('gagamultiobj()', ...
    @(obj_fn, phi_init) ga_gamultiobj(obj_fn, length(phi_init), [], [], [], [], [], [], ...
    gaoptimset('TolFun', TolFun, ...
    'PopulationSize', PopSize, ...
    'Generations', Gen, ...
    'Vectorized', 'on', ...
    'StallGenLimit', StallGen, ...
    'ParetoFraction', 0.5, ... % higher elitism to keep as long as possible the previous optimal solutions
    'MutationFcn', @ga_mutation, ...
    'CrossoverFraction', 0.2, ...
    'CrossoverFcn', @ga_crossover, ...
    'DistanceMeasureFcn', @ga_distancecrowding, ...
    'InitialPopulation', [best_indivs; rand(PopSize - size(best_indivs,1),length(phi_init)).*repmat(GA_UpperBound-GA_LowerBound, PopSize - size(best_indivs,1),1) + repmat(GA_LowerBound, PopSize - size(best_indivs,1),1)])));

unpack_fittables(term_phi);

% NOW DO THE LAST FIT05C STEPS:

fit_iteratively(@step_until_10neg3, ...
    create_term_fn());
% fit_iteratively(@step_until_10neg5, ...
%     create_term_fn());
% fit_iteratively(@step_until_10neg55, ...
%     create_term_fn());
% [~, term_score, ~] = fit_iteratively(@step_until_10neg6, ...
%     create_term_fn());

% compute the last score - update if necessary
phi_candidate = pack_fittables(STACK);
calc_xxx(1);
[m, p] = META.perf_metric();
if best_score > m+p
    % we have a winner
    best_indiv = phi_candidate;
    best_score = m+p;
    best_origin = ['aggregated last genetic fit'];
    fprintf('\n\n ===== NEW BEST SOLUTION ===== SCORE is %f\n\n', best_score);
end

fprintf('\nBest Aggregated Score = %f\n%s',m+p,num2str(phi_candidate'));
for i=1:size(best_indivs,1)
    fprintf('\n\nScore = %f\n%s',best_scores,num2str(best_indivs(i,:)));
end




unpack_fittables(best_indiv);calc_xxx(1);

fprintf('\n\n\n\n\n!!!!!\nOVERALL BEST SOLUTION COMES FROM %s\nOVERALL BEST SCORE is %f\n!!!!!\n\n\n\n\n', best_origin, best_score);





end