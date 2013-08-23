function [termcond, n_iters] = fit_scaat(options)

global STACK;
phi = pack_fittables(STACK);
n_params = length(phi);

default_options = struct('InitStepsize', 1.0, ...
                         'StepsPerParam', 1, ...
                         'OuterLoops', 100, ...
                         'StepGrowth', 2.0, ...
                         'StepShrink', 0.5, ...
                         'StopStepsize', 10^-7);

if ~exist('options', 'var')
    optz = default_options;
else
    optz = merge_structs(default_options,options);
end

function [x, s, termcond] = one_at_a_time(objfn, x_0, opts)        
    
    stepsizes = opts.InitStepsize * ones(size(x_0));
    x = x_0; 
    s = objfn(x); 
    n_params = length(x_0);
    
	fprintf('Initial score: %d \tphi:', s);
    p = pack_fittables(STACK);
    fprintf(' %.3f', p);
    fprintf('\n');
    
    for ii = 1:opts.OuterLoops
        fprintf('Outer loop %d/%d\n', ii, opts.OuterLoops);
        for jj = 1:n_params
            stepsize = stepsizes(jj);
            
            fprintf('Param# %03d/%03d, Score: %d, Stepsize: %d', jj, n_params, s, stepsize);
            p = pack_fittables(STACK);
            fprintf(' %.3f', p);
            fprintf(']\n');
            
            n = 1;
            while n <= opts.StepsPerParam && stepsize > opts.StopStepsize 
                
                stepdir = zeros(n_params, 1);
                stepdir(jj) = stepsize;
                
                x_next = x;
                s_next = s;
                
                x_fwd = x + stepdir;
                x_bck = x - stepdir;            
            
                s_fwd = objfn(x_fwd);
                s_bck = objfn(x_bck);         
                
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
                    n = n+1;
                end
            end            
            stepsizes(jj) = stepsize;
        end 
        fprintf('Score: %d \tphi:', s);
        p = pack_fittables(STACK);
        fprintf(' %.3f', p);
        fprintf('\n');
    end        
    termcond = 0;
end

[termcond, n_iters] = default_fitter_loop('fit_scaat()', ...
    @(obj_fn, phi_init) one_at_a_time(obj_fn, phi_init, optz), true);

end