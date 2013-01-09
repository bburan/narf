function fit_with_lsqcurvefit()
% Fits any fittable fields using the training set data only. 
% REQUIRES that there be a field 'prediction' defined at the end of the XXX
% data structure, which will be compared with the 'respavg' signal to
% compute the least squares error used to do the fitting.
% 
% TODO: Pass this function options for lsqcurvefit, and UB and LB

global XXX STACK;
cnt = 1;  % For printing progress dots

recalc_xxx(1); 

function error = my_fitter(phi, start_depth)
    % Perform the nonlinear fitting routine, starting at start_depth and
    % recalculating until the end of the XXX datastructure. 
        
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
    
    % Concatenate all the training set data into a large vector
    pred    = flatten_field(XXX{end}.dat, XXX{1}.training_set, 'prediction');
    respavg = flatten_field(XXX{end}.dat, XXX{1}.training_set, 'respavg');
    
    % Eliminate any MSE where there is a NaN in respavg 
    % I would have preferred to just excise the NaNs, but lsqcurvefit
    % works on a constant length vector for X and Y, so that's not
    % possible.
    % TODO: Consider nlinfit instead. 
    respavg(isnan(respavg)) = pred(isnan(respavg));
    
    error = pred - respavg;
end

% The optimization
phi_init = pack_fittables(STACK);
LB = [];
UB = [];
fprintf('Fitting %d variables with lsqcurvefit()\n', length(phi_init));

options = optimset('MaxIter', 1000, ...
                   'MaxFunEvals', 5000, ... % 100*numel(phi_init), ...
                   'TolFun', 1e-12, 'TolX', 1e-9);

len = length(flatten_field(XXX{end}.dat, XXX{1}.training_set, 'respavg'));
start_depth = 2; % TODO: Maybe not all fits should start here by default?
phi_best = lsqcurvefit(@my_fitter, phi_init, start_depth,  ...
               zeros(len, 1), LB, UB, options);
unpack_fittables(phi_best);
end