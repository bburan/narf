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
%    options  Options  passed to lsqcurvefit() to help fit. Defaults are:
%                MaxIter 1000
%                MaxFunEvals 5000
%                TolFun 1e-12
%                TolX 1e-9
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
    options = optimset('MaxIter', 1000, ...
                       'MaxFunEvals', 1000, ...
                       'TolFun', 1e-12, 'TolX', 1e-9);  
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
    recalc_xxx(start_depth);
    
    % Concatenate training set prediction and reality into a long vector
    pred    = flatten_field(XXX{end}.dat, XXX{1}.training_set, field1);
    respavg = flatten_field(XXX{end}.dat, XXX{1}.training_set, field2);       
    
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

len = length(flatten_field(XXX{end}.dat, XXX{1}.training_set, 'respavg'));
start_depth = find_fit_start_depth(STACK);
[phi_best, ~, ~, termcond] = lsqcurvefit(...
                                    @my_obj_fn, phi_init, start_depth,  ...
                                    zeros(len, 1), LB, UB, options);
unpack_fittables(phi_best);
end
