function [term_cond, best_score, n_iters, options] = fit_scaat(varargin)

global STACK;
phi = pack_fittables(STACK);
n_params = length(phi);

default_options = struct('StepsPerParam', 1, ...
                         'InitStepSize', 1.0, ...
                         'StopAtStepsize', 10^-6, ...
                         'StepGrowth', 2.0, ...
                         'StepShrink', 0.5);

if (length(varargin) == 1)
    options = merge_structs(default_options, varargin{1});
else
    options = merge_structs(default_options, struct(varargin{:}));
end

if ~isfield(options, 'TermFn') 
    options.TermFn = create_term_fn(options);
end

[term_cond, best_score, n_iters, term_stepsizes] = ...
    default_fitter_loop('fit_scaat()', ...
        @(obj_fn, phi_init) one_at_a_time(obj_fn, phi_init, options), true);

% If we are running this in an iterative fitter, next time start from here
options.InitStepSizeVector = term_stepsizes;

% ---------------

function [x, s, term_cond, stepsizes] = one_at_a_time(objfn, x_0, opts)        
    
    if ~isfield(opts, 'InitStepSizeVector') 
        stepsizes = opts.InitStepSize * ones(size(x_0));
    else
        stepsizes = opts.InitStepSizeVector;
    end
    x = x_0; 
    [s, o] = objfn(x); 
    n_params = length(x_0);        
    n = 1;
    d = inf;    
    s_outer = s;
    fprintf('\n');
    while ~(opts.TermFn(n, s_outer, d, o))
        fprintf('SCAAT Step #%d (Score: %e)\n', n, s_outer);
        for jj = 1:n_params
            stepsize = stepsizes(jj);
            
            % Comment out if desired
            fprintf('eval: %d, size: %.3e, coef#%3d, s_delta: %.3e, score:%e', o, stepsize, jj, d, s);
            
            ii = 1;
            while ii <= opts.StepsPerParam && stepsize > opts.StopAtStepsize 
                
                stepdir = zeros(n_params, 1);
                stepdir(jj) = stepsize;
                
                x_next = x;
                s_next = s;
                
                x_fwd = x + stepdir;
                x_bck = x - stepdir;            
            
                [s_fwd, o] = objfn(x_fwd);
                [s_bck, o] = objfn(x_bck);         
                
                if s_fwd < s
                    s_next = s_fwd;
                    x_next = x_fwd;
                end
                
                if s_bck < s
                    s_next = s_bck;
                    x_next = x_bck;
                end
                
                if all(x == x_next)
                    stepsize = stepsize * opts.StepShrink;
                    continue;
                else
                    stepsize = stepsize * opts.StepGrowth;
                    
                    x = x_next;
                    s = s_next;
                    ii = ii + 1;
                end
            end
            
            fprintf('\n');
            stepsizes(jj) = stepsize;
        end 
        n = n+1;
        d = s_outer - s;
        s_outer = s;
    end 
    term_cond = 0;
end

end