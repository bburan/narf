function fit15rnd()


% RATIONALE:
% with low-parametric kernels, the default starting position may lead to a
% local minimum quickly (with bad luck and wrong initial changes of step
% size, the starting position may even be the local minimum!

% STUDIED SOLUTION:
% Employing 10 different random starting positions to perform a preliminary
% search. Keeping the best position for a real search (exactly as in fit05)

% RESULT:
% huge improvement for Lindeberg kernels that were stuck at the initial
% position

global STACK

nmse();

MaxStepsPerIteration=10;
StepGrowth=1.1;

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

    function [a,b,c,d] = step_until_10neg6(prev_opts)
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_boo(prev_opts);
        else
            [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-6, ...
                'StopAtStepNumber', MaxStepsPerIteration, ...
                'StepGrowth', StepGrowth);
        end
    end

    function [term_cond, term_score,  n_iters, options] = quick_search()
        
        % [term_cond, term_score,  n_iters, options] = ...
        %     fit_boo('StopAtAbsScoreDelta', 10^-2.0, ...
        %     'StepSize', 10, ...
        %     'StepGrowth', 1.1, ...
        %     'StepShrink', 0.5);
        
        %         [term_cond, term_score,  n_iters, options] = ...
        %             fit_boo('StepSize', 10, ...
        %             'StepGrowth', 1.1, ...
        %             'StepShrink', 0.5, ...
        %             'StopAtAbsScoreSize', 10^-1.0);
        
        [term_cond, term_score,  n_iters, options] = ...
            fit_boo('StepSize', 1, ...
            'StepGrowth', 1.1, ...
            'StepShrink', 0.5, ...
            'StopAtStepNumber', 50);
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
                    strcmp(module_being_fit, 'weight_channels')
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
                        'StopAtStepNumber', 1);
                end
            end
        end
        
        fn = @subfitter;
        
    end




% initialization: a first rough search
initial_stack = STACK;

[term_cond, term_score,  n_iters, options] = quick_search();

for i=1:10,
    % Store the previous state of the STACK and revert to default
    cached_stack = STACK;
    STACK = initial_stack;
    
    % Randomize variables in the fit fields
    for kk = 1:length(STACK)
        if ~isfield(cached_stack{kk}{1}, 'fit_fields')
            continue;
        else
            for ll = 1:length(STACK{kk}{:})
                if (strcmp(STACK{kk}{ll}.name, 'lindeberg_filter')),
                    mdl = STACK{kk}{ll};
                    params = zeros(7,1);
                    params(1) = (0 + (mdl.num_dims-0)).*rand(1); %xshift
                    params(2) = (0 + (mdl.num_coefs/4-0)).*rand(1); %tshift
                    params(3) = (0.1 + (mdl.num_dims/2-0.1)).*rand(1); %s
                    params(4) = (0.5 + (mdl.num_coefs/3-0.5)).*rand(1); % tau
                    params(5) = (-0.75 + 1.5).*rand(1); % v
                    params(6) = 10.*rand(1); % norm_factor
                    % add_factor is mdl.lincoefs(7) and is set to 0
                    
                    cc = 'lincoefs';
                    STACK{kk}{ll}.(char(cc)) = params;
                    cc = 'baseline';
                    STACK{kk}{ll}.(char(cc)) = 0;
                    
                    %                     for cc=STACK{kk}{ll}.fit_fields
                    %                         STACK{kk}{ll}.(char(cc)) = ...
                    %                             rand(length(STACK{kk}{ll}.(char(cc))),1)*30.0-10.0;
                    %                     end
                end
            end
        end
    end
    
    [a,b,c,d] = quick_search();
    
    fprintf('======================================================================\n');
    fprintf('======================================================================\n');
    fprintf('=================       Finished prelim search num %d       ==========\n', i);
    
    if term_score<b,
        % no improvement => revert
        STACK = cached_stack;
        fprintf('=================                                           ==========\n');
        fprintf('=================                                           ==========\n');
        fprintf('=================        Reverting to old solution          ==========\n');
        fprintf('======================================================================\n');
        fprintf('======================================================================\n');
        
    else
        % improvement => keep new
        term_score = b;
        fprintf('=================                                           ==========\n');
        fprintf('=================                                           ==========\n');
        fprintf('=================           Keeping new solution            ==========\n');
        fprintf('======================================================================\n');
        fprintf('======================================================================\n');
        
    end
end

% Now gradually shrink the stopping criterion - as in fit10
scale=10^1;
stop_at=10^-6;

while(scale > stop_at)
    fit_iteratively(make_subfitter(scale), create_term_fn('StopAtAbsScoreDelta', scale));
    scale = scale * 0.7; % Very conservative: 0.8. Probably 0.5 or even 0.25 is fine.
end


end
