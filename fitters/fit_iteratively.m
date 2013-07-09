function [termcond, n_iters] = fit_iteratively(fitter, n_outer_loops)

global STACK;

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;  
    return 
end

if ~exist('fitter', 'var')
    fitter = @(~) fit_boost(length(phi_init));
end

if ~exist('n_outer_loops', 'var')
    n_outer_loops = 10;
end

cached_stack = STACK;

for ii = 1:n_outer_loops,
    fprintf('Outer Loop Iteration %d/%d\n', ii, n_outer_loops);
    for jj = 1:length(STACK),
        if ~isfield(cached_stack{jj}{1}, 'fit_fields')
            continue;
        end
        
        fprintf('Running fitter on only module %d\n', jj);
        
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
        
        % Run the fitter once
        fitter();
        
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
end

n_iters = NaN;