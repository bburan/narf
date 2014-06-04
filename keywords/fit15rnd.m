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


% initialization: a first rough search
initial_stack = STACK;

[term_cond, term_score,  n_iters, options] = quick_search();

for i=1:20,
    % Store the previous state of the STACK and revert to default
    cached_stack = STACK;
    STACK = initial_stack;
    
    % Randomize variables in the fit fields
    for kk = 1:length(STACK)
        if ~isfield(cached_stack{kk}{1}, 'fit_fields')
            continue;
        else
            for ll = 1:length(STACK{kk}{:})
                if (strcmp(STACK{kk}{ll}.name, 'lindeberg_filter'))
                    mdl = STACK{kk}{ll};
                    params = zeros(7,1);
                    % start with a reasonable guess
                    params(1) = mdl.num_dims.*rand(1); %xshift
                    params(2) = (mdl.num_coefs/4).*rand(1); %tshift
                    params(3) = (0.1 + (mdl.num_dims/2-0.1).*rand(1)); %s
                    params(4) = (0.5 + (mdl.num_coefs/3-0.5).*rand(1)); % tau
                    params(5) = 2.*(rand(1)-0.5); % v in [-1;1]
                    params(6) = 10.*(rand(1)-0.5); % norm_factor in [-5;5]
                    % add_factor is mdl.lincoefs(7) and is set to 0
                    
                    cc = 'lincoefs';
                    STACK{kk}{ll}.(char(cc)) = params;
                    cc = 'baseline';
                    STACK{kk}{ll}.(char(cc)) = 0;
                    
                    %                     for cc=STACK{kk}{ll}.fit_fields
                    %                         STACK{kk}{ll}.(char(cc)) = ...
                    %                             rand(length(STACK{kk}{ll}.(char(cc))),1)*30.0-10.0;
                    %                     end
                elseif (strcmp(STACK{kk}{ll}.name, 'lindeberg_spectral'))
                    % start with a reasonable guess
                    STACK{kk}{ll}.bf = STACK{kk}{ll}.num_channels.*rand(1); %xshift
                    STACK{kk}{ll}.s = (0.1 + (STACK{kk}{ll}.num_channels/2-0.1).*rand(1)); %s
                    STACK{kk}{ll}.norm_factor = 10.*(rand(1)-0.5); % norm_factor in [-5;5]
                    STACK{kk}{ll}.add_factor = 0; % add_factor is set to 0
                elseif (strcmp(STACK{kk}{ll}.name, 'weight_channels'))
                    % start with a reasonable guess
                    STACK{kk}{ll}.weights = rand(size(STACK{kk}{ll}.weights)); %(order of magnitude 1)
                    STACK{kk}{ll}.y_offset = rand(size(STACK{kk}{ll}.y_offset))/100; %(order of magnitude 0.01)
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

% fit_iteratively(@step_until_10neg6, ...
%                create_term_fn());
end
