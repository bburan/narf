function termfn = create_term_fn(varargin)
% termfn = create_term_fn(varargin)
%
% Returns a termination function based on the arguments that you provide. 
% Intended for use with various fitters to avoid code bloat. See fit_boo,
% fit_scaat, etc. for examples. 
%
% The following parameters are expected:
%
%   'StopAtSeconds'        Stop after this many seconds. [Default 60*60*24]
%   'StopAtStepNumber'     Stop after this many steps. [10^6]
%   'StopAtObjFnNumber'    Stop after this many obj fn evaluations. [10^8]
%   'StopAtStepSize'       Stop when stepsize is this small. [10^-9]
%   'StopAtRelScoreSize'   Stop when the improvement to the score relative
%                          to the initial few improvements is this small.
%   'StopAtAbsScoreSize'   Stop when the absolute improvement to the score 
%                          is smaller than this level.
%   'RelDeltaSamples'      This is the number of samples used to compute
%                          the initial score delta for StopAtRelScoreSize.
%
% Keep in mind that StopAtRelScoreSize uses closed-over variables. If
% want to 'reset' the closed-over variables, you will need to re-create
% the termination function using create_term_fn(). 
%
% The returned TERMFN accepts 3 arguments:
%                     n    Step number of this iteration.
%                     s    Stepsize taken in this past step
%                     d    Score improvement (delta) vs previous point
%
% Boosting will terminate when TermFn returns a true.

default_options = struct('StopAtSeconds', 60*60*24, ...
                         'StopAtStepNumber', 10^6, ...
                         'StopAtObjFnNumber', 10^8, ...
                         'StopAtStepSize', 10^-9, ...
                         'StopAtRelScoreDelta', 10^-12, ...
                         'StopAtAbsScoreDelta', 10^-12, ...
                         'RelStopSamples', 5);
if (length(varargin) == 1)
    opts = merge_structs(default_options, varargin{1});
else
    opts = merge_structs(default_options, struct(varargin{:}));
end

t_start = clock();
first_deltas = [];

function stopcode = default_termfn(n,s,d,o)
    if (n > opts.StopAtStepNumber),   
        stopcode = 1; 
        fprintf('Took %d steps in total, stopping.\n', ...
            opts.StopAtStepNumber);
        return;
    end
    if (s < opts.StopAtStepSize),      
        stopcode = 2;  
        fprintf('Stepsize less than absolute limit %e, stopping.\n', ...
            opts.StopAtStepSize);
        return; 
    end
    if (d < opts.StopAtAbsScoreDelta), 
        stopcode = 3;
        fprintf('Absolute score improvement less than %e, stopping.\n', ...
            opts.StopAtAbsScoreDelta);
        return; 
    end
    
    if length(first_deltas) < 1+opts.RelStopSamples && ~isnan(d)
        first_deltas(end+1) = d;
    elseif isnan(d)
        % Do nothing
    else
        % Test if we need to stop
        if (d / mean(first_deltas(2:end)) < opts.StopAtRelScoreDelta)
            stopcode = 4;
            fprintf('Relative score improvement less than %e, stopping.\n', ...
                opts.StopAtRelScoreDelta);
            return
        end
    end
       
    if (etime(clock, t_start) > opts.StopAtSeconds)
        stopcode = 5; 
        fprintf('Fit time exceeded %e, stopping.\n', ...
                opts.StopAtSeconds);
        return
    end
    
    if (o > opts.StopAtObjFnNumber), 
        stopcode = 6;
        fprintf('Total number of iterations more than %d, stopping.\n', ...
            opts.StopAtObjFnNumber);
        return; 
    end

    stopcode = false; 
end

termfn = @default_termfn;

end

