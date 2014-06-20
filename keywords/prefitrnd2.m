function prefitrnd()
% JL - June 2015

% RATIONALE:
% with low-parametric kernels, the default starting position may lead to a
% local minimum quickly (with bad luck and wrong initial changes of step
% size, the starting position may even be the local minimum!

% STUDIED SOLUTION:
% Trying 50 different random starting positions to perform a preliminary
% search. Keeping the best position for a real search, to happen later on.


global STACK

nmse();

    function [term_cond, term_score,  n_iters, options] = quick_search()
        [term_cond, term_score,  n_iters, options] = ...
            fit_boo('StepSize', 1, ...
            'StepGrowth', 1.1, ...
            'StepShrink', 0.5, ...
            'StopAtStepNumber', 10);
    end

% initialization: a first rough search with the default starting position
initial_stack = STACK;
[~, term_score,  ~, ~] = quick_search();

% random search: try 50 different positions
for i=1:50
    % Store the previous state of the STACK and revert to default
    cached_stack = STACK;
    STACK = initial_stack;
    
    % Randomize variables in the fit fields
    for kk = 1:length(STACK)
        if ~isfield(cached_stack{kk}{1}, 'fit_fields')
            continue;
        else
            for ll = 1:length(STACK{kk}{:})
                % different random initialization depending on the filter
                if (strcmp(STACK{kk}{ll}.name, 'lindeberg_filter'))
                    % LINDEBERG FILTER
                    mdl = STACK{kk}{ll};
                    params = zeros(7,1);
                    params(1) = mdl.num_dims.*rand(1); %xshift
                    params(2) = (mdl.num_coefs/4).*rand(1); %tshift
                    params(3) = (0.1 + (mdl.num_dims/2-0.1).*rand(1)); %s
                    params(4) = (0.5 + (mdl.num_coefs/3-0.5).*rand(1)); % tau
                    params(5) = 2.*(rand(1)-0.5); % v in [-1;1]
                    params(6) = 10.*(rand(1)-0.5); % norm_factor in [-5;5]
                    STACK{kk}{ll}.('lincoefs') = params;
                    STACK{kk}{ll}.('baseline') = 0;
                elseif (strcmp(STACK{kk}{ll}.name, 'lindeberg_spectral'))
                    % LINDEBERG SPECTRAL-ONLY
                    STACK{kk}{ll}.bf = STACK{kk}{ll}.num_channels.*rand(1); %xshift
                    STACK{kk}{ll}.s = (0.1 + (STACK{kk}{ll}.num_channels/2-0.1).*rand(1)); %s
                    STACK{kk}{ll}.norm_factor = 10.*(rand(1)-0.5); % norm_factor in [-5;5]
                    STACK{kk}{ll}.add_factor = 0; % add_factor is set to 0
                elseif (strcmp(STACK{kk}{ll}.name, 'weight_channels'))
                    % WEIGHT CHANNELS
                    STACK{kk}{ll}.weights = rand(size(STACK{kk}{ll}.weights)); %(order of magnitude 1)
                    STACK{kk}{ll}.y_offset = rand(size(STACK{kk}{ll}.y_offset))/100; %(order of magnitude 0.01)
                elseif (strcmp(STACK{kk}{ll}.name, 'fir_filter'))
                    % FIR
                    STACK{kk}{ll}.coefs = normrnd(0, 5, size(STACK{kk}{ll}.coefs)); %(0,5) gaussian random initialization
                    STACK{kk}{ll}.('baseline') = 0;
                end
            end
        end
    end
    
    [~,new_term_score,~,~] = quick_search();
    
    fprintf('======================================================================\n');
    fprintf('======================================================================\n');
    fprintf('=================       Finished prelim search num %d       ==========\n', i);
    
    if term_score < new_term_score,
        % no improvement => revert
        STACK = cached_stack;
        fprintf('=================                                           ==========\n');
        fprintf('=================                                           ==========\n');
        fprintf('=================        Reverting to old solution          ==========\n');
        fprintf('======================================================================\n');
        fprintf('======================================================================\n');
        
    else
        % improvement => keep new (do nothing but update the best score)
        term_score = new_term_score;
        fprintf('=================                                           ==========\n');
        fprintf('=================                                           ==========\n');
        fprintf('=================           Keeping new solution            ==========\n');
        fprintf('======================================================================\n');
        fprintf('======================================================================\n');
        
    end
end

end
