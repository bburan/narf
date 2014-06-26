function [termcond, s, n_iters] = fit_iteratively(fitter, outer_loop_term_fn)

global STACK;

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;  
    return 
end

if ~exist('fitter', 'var')
    error('You should provide a fitter to fit_iteratively();'); 
end

if ~exist('outer_loop_term_fn', 'var')
    outer_loop_term_fn = create_term_fn('StopAtStepNumber', 10);
end

cached_stack = STACK;
prev_opts = {}; 

n = 1;
s = nan;
d = nan;
n_iters = 0;

fprintf('*********************BEGINNING ITERATIVE FIT**************************\n');

% Terminate if the previous position is exactly the same as this one
termcond = outer_loop_term_fn(n, s, d, n_iters);
while ~termcond                 
    fprintf('Beginning Iterative Loop Iteration #%d\n', n);
    for jj = 1:length(STACK),
        if ~isfield(cached_stack{jj}{1}, 'fit_fields')
            prev_opts{jj} = nan; % Just a stupid placeholder
            continue;
        end
        for mm=1:length(STACK{jj}),
            fprintf('Running fitter on only module %d:%d (%s)\n', ...
                    jj, mm, STACK{jj}{mm}.name);
        
            % Erase all fit fields except the normal one
            for kk = 1:length(STACK)
                for ll=1:length(STACK{kk}),
                    if ~isfield(cached_stack{kk}{ll}, 'fit_fields') || ...
                            (jj==kk && mm==ll),
                        continue;
                    else
                        STACK{kk}{ll}.fit_fields = {};
                    end
                end
            end
            
            % Run the fitter once, passing it the previous opts if they exist.       
            if length(prev_opts) < jj
                [~, s_new, iters, prev_opts{jj}] = fitter();
            else
                [~, s_new, iters, prev_opts{jj}] = fitter(prev_opts{jj});
            end        
            
            if isnan(s)
                s = s_new;   % Just the first time, initialize it
            end
            
            % Restore the fit fields
            for kk = 1:length(STACK)
                if ~isfield(cached_stack{kk}{1}, 'fit_fields')
                    continue;
                else
                    for ll = 1:length(STACK{kk}(:))
                        STACK{kk}{ll}.fit_fields = cached_stack{kk}{ll}.fit_fields;
                    end
                end
            end
        end
    end
    
    n_iters = n_iters + iters;
    if isnan(s)
        d = s_new;
    else
        d = s - s_new;
    end
    s = s_new;
    
    fprintf('Ending Iterative loop %d. Score: %e, S_delta: %e\n', n, s, d);
    fprintf('======================================================================\n');  
    n = n + 1;
    termcond = outer_loop_term_fn(n, s, d, n_iters);
    
    % tell the queue daemon that job is still alive
    dbtickqueue;
end
