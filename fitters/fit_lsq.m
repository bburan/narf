function termcond = fit_lsq(field1, field2, options)
% termcond = fit_lsq(field1, field2, options)
%
% Fits all parameters in STACK marked with 'fit_fields' such that signals
% FIELD1 and FIELD2 have a least-squared error that is minimized. Uses only
% data in XXX.training_set.
%
% ARGUMENTS:
%    field1   The name of the first signal
%             Defaults to 'stim' if no argument is passed
%    field2   The name of the second signal
%             Defaults to 'respavg' if no argument is passed
%    options  Options  passed to lsqcurvefit() to help fit.
%
% RETURNS:
%    termcond    Termination condition of optimization.
%
% TODO: There should be a way to add upper and lower bounds on search.

global XXX STACK;
cnt = 1;  % For printing progress dots

if nargin < 2
    field1 = 'stim';
    field2 = 'respavg';
end
if nargin < 3
    options = optimset('MaxIter', 3000, ...
                       'MaxFunEvals', 3000, ...
                       'TolFun', 1e-12, 'TolX', 1e-9);  
end

% Check if mse exists and has a sparseness penaltiy
% If so, append a bogus element to the objective target vector and 
% prediction vector in order to mimic the effect of a sparseness penalty.

[~, idxes] = find_modules(STACK, 'mean_squared_error');
if ~isempty(idxes) && STACK{idxes(1)}.sparseness_weight ~= 0
    append_penalty = true;
    sparseness_weight = STACK{idxes(1)}.sparseness_weight;
else
    append_penalty = false;
end

function error = my_obj_fn(phi, start_depth)
    % Compute the stack values from start_depth to the end of the XXX data 
    % structure and computing the L1 error between FIELD1 and FIELD2
        
    % Print 1 progress dot for every 20 iterations, and newline every 1000
    if isequal(mod(cnt, 1000), 1)
        fprintf('\n[Iteration %d]', cnt);
    end
    if isequal(mod(cnt, 20), 1)
        fprintf('.');
    end
    cnt = cnt + 1;

    unpack_fittables(phi);
    calc_xxx(start_depth);
    
    % Concatenate training set prediction and reality into a long vector
    pred    = flatten_field(XXX{end}.dat, XXX{end}.training_set, field1);
    respavg = flatten_field(XXX{end}.dat, XXX{end}.training_set, field2);       
    
    % Set error to zero where there is a NaN in respavg 
    % I would have preferred to just excise the NaNs, but lsqcurvefit
    % works on a constant length vector for X and Y, so that's not
    % possible.    
    % Make the error zero for any predictions with NaNs
    respavg(isnan(pred)) = 0; 
    pred(isnan(pred)) = 0;
    % Make the error zero for any respavg's with NaNs
    respavg(isnan(respavg)) = pred(isnan(respavg));
    
    % This will be compared with zero and squared to calc MSE
    error = pred - respavg; 
    
    if append_penalty
        % TODO: Square root or not?
        penalty = sparseness_weight * sparsity_metric(phi); 
%        fprintf('lsq metric: %f\n', sparsity_metric(phi));
        N = length(error);
        esum = nansum(error.^2);
        z = sqrt( (N+1)*((1/N)*esum + penalty) - esum);
        error(end+1) = z;
        %fprintf('lsq: %f\n', nanmean(error.^2)); 
    end
end

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

LB = [];
UB = [];

fprintf('Fitting %d variables with lsqcurvefit()\n', length(phi_init));

len = length(flatten_field(XXX{end}.dat, XXX{end}.training_set, 'respavg'));
start_depth = find_fit_start_depth(STACK);
if append_penalty
    fprintf('Appending penalty!');
    [phi_best, ~, ~, termcond] = lsqcurvefit(...
                                    @my_obj_fn, phi_init, start_depth,  ...
                                    zeros(len+1, 1), LB, UB, options);  
else
    [phi_best, ~, ~, termcond] = lsqcurvefit(...
                                    @my_obj_fn, phi_init, start_depth,  ...
                                    zeros(len, 1), LB, UB, options);
end
unpack_fittables(phi_best);
end
