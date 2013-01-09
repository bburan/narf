function fit_with_lsqcurvefit()
% Fits any fittable fields using the training set data only. 
% 
% TODO: Pass this function options for lsqcurvefit, and UB and LB

global XXX STACK;
cnt = 0;  % For printing progress dots

recalc_xxx(1); 

function error = my_fitter(phi, start_depth)
    % Perform the nonlinear fitting routine, starting at start_depth and
    % recalculating until the end of the XXX datastructure. 
    
    cnt = cnt + 1;
    % Print 1 progress dot for every 20 iterations, and newline every 1000
    if isequal(mod(cnt, 1000), 1)
        fprintf('\n[Iteration %d]', cnt);
    end
    if isequal(mod(cnt, 20), 1)
        fprintf('.');
    end
    
    unpack_fittables(phi);
    recalc_xxx(start_depth);
    prediction = XXX{end}.dat.(XXX{1}.training_set{1}).stim(:);
    respavg = XXX{end}.dat.(XXX{1}.training_set{1}).respavg(:);
       
    % Eliminate any MSE where there is a NaN in respavg 
    % I would have preferred to just shorten the length of it, but
    % lsqcurvefit needs a constant-length vector passed to it. 
    % TODO: Consider nlinfit instead. 
    respavg(isnan(respavg)) = prediction(isnan(respavg));
    
    error = prediction - respavg;
end

% The optimization
phi_init = pack_fittables(STACK);
LB = [];
UB = [];
fprintf('Fitting %d variables with lsqcurvefit()\n', length(phi_init));

options = optimset('MaxIter', 1000, ...
                   'MaxFunEvals', 5000, ... % 100*numel(phi_init), ...
                   'TolFun', 1e-12, 'TolX', 1e-9);
phi_best = lsqcurvefit(@my_fitter, phi_init, 2,  ...
               zeros(size(XXX{end}.dat.(XXX{1}.training_set{1}).respavg(:))), ...
               LB, UB, options);
unpack_fittables(phi_best);
end