function termcond = fit_with_jacklsq(field1, field2, n_jacks, options)
% termcond = fit_with_jacklsq(field1, field2, n_jacks, options)
%
% Fits all parameters in STACK marked with 'fit_fields' such that signals
% FIELD1 and FIELD2 have a least-squared error that is minimized. Uses only
% data in XXX.training_set. Divides the data into N_JACKS groups of equally 
% sized data subsets, fits the parameters to each of those subsets, and
% then takes the mean of the parameters found with those subsets.
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

if nargin < 3
    field1 = 'stim';
    field2 = 'respavg';
    n_jacks = 10;
end
if nargin < 4
    options = optimset('MaxIter', 1000, ...
                       'MaxFunEvals', 1000, ...
                       'Display', 'off',... 
                       'TolFun', 1e-9, ...
                       'TolX', 1e-9);  
end

function error = my_obj_fn(phi, jackjackjack)
    % Compute the stack values from start_depth to the end of the XXX data 
    % structure and computing the L1 error between FIELD1 and FIELD2
        
    depth = jackjackjack{1};
    n = jackjackjack{2};
    idxs = jackjackjack{3};
    
    % Print 1 progress dot for every 20 iterations, and newline every 1000
    if isequal(mod(cnt, 1000), 1)
        fprintf('\n[Iteration %d]', cnt);
    end
    if isequal(mod(cnt, 20), 1)
        fprintf('.');
    end
    cnt = cnt + 1;

    unpack_fittables(phi);
    recalc_xxx(depth);
    
    % Concatenate training set prediction and reality into a long vector
    pred    = flatten_field(XXX{end}.dat, XXX{1}.training_set, field1);
    respavg = flatten_field(XXX{end}.dat, XXX{1}.training_set, field2);       
    
    % Set error to zero where there is a NaN in respavg 
    % I would have preferred to just excise the NaNs, but lsqcurvefit
    % works on a constant length vector for X and Y, so that's not
    % possible.
    respavg(isnan(respavg)) = pred(isnan(respavg));
    
    % This will be compared with zero and squared to calc MSE
    error = pred - respavg; 
    
    % Zero out the part which we don't want to use
    error(idxs(n):idxs(n+1)) = 0;
    
end

phi_init = pack_fittables(STACK);

if isempty(phi_init)
    fprintf('Skipping because there are no parameters to fit.\n');
    termcond = NaN;
    return 
end

LB = [];
UB = [];
recalc_xxx(1); 

len = length(flatten_field(XXX{end}.dat, XXX{1}.training_set, 'respavg'));
jack_idxs = floor(linspace(1, len, n_jacks+1));
start_depth = find_fit_start_depth(STACK);

for jj = 1:n_jacks
    fprintf('\nFitting %d variables with lsq jackknife [%d/%d]\n', ...
            length(phi_init), jj, n_jacks);
    unpack_fittables(phi_init); % Reset where you start from every time
    phi_jack(:, jj) = lsqcurvefit(@my_obj_fn, phi_init, {start_depth jj jack_idxs},  ...
                                    zeros(len, 1), LB, UB, options);
end                            
phi_best = mean(phi_jack, 2);
termcond = nan;
unpack_fittables(phi_best);
end