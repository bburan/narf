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

% Terminate if the previous position is exactly the same as this one
termcond = outer_loop_term_fn(n, s, d, n_iters);
while ~termcond
    fprintf('Outer Loop Iteration #%d. Score: %e\n', n, s);
    for jj = 1:length(STACK),
        if ~isfield(cached_stack{jj}{1}, 'fit_fields')
            prev_opts{jj} = nan; % Just a stupid placeholder
            continue;
        end
        
        fprintf('Running fitter on only module %d (%s)\n', jj, STACK{jj}{1}.name);
        
        % Erase all fit fields except the normal one
        for kk = 1:length(STACK)
            if ~isfield(cached_stack{kk}{1}, 'fit_fields') || jj==kk
                continue;
            else
                for ll = 1:length(STACK{kk}{:})
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
        
        n_iters = n_iters + iters; 
        if isnan(s)
            d = s_new;
        else
            d = s - s_new;
        end
        s = s_new;
        
        % Restore the fit fields
        for kk = 1:length(STACK)
            if ~isfield(cached_stack{kk}{1}, 'fit_fields')
                continue;
            else
                for ll = 1:length(STACK{kk}{:})
                    STACK{kk}{ll}.fit_fields = cached_stack{kk}{ll}.fit_fields;
                end
            end
        end
    end    
    fprintf('Iterative loop %d ended. Score: %e, S_delta: %e\n', n, s, d);
    n = n + 1;
    termcond = outer_loop_term_fn(n, s, d, n_iters);
end
