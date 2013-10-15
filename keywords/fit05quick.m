function fit05quick()
% helper for perfile fits. just like fit05 but stop at 10^-4 error.
    
nmse();

MaxStepsPerIteration=100;

function [a,b,c,d] = step_until_10neg3(prev_opts)              
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', 1.3);
    end
end

function [a,b,c,d] = step_until_10neg35(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3.5, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', 1.3);
    end
end

function [a,b,c,d] = step_until_10neg4(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-4, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', 1.3);
    end
end

fit_boo('StopAtAbsScoreDelta', 10^-2, ...
        'StepGrowth', 1.3);

fit_iteratively(@step_until_10neg3, ...
                create_term_fn());

fit_iteratively(@step_until_10neg35, ...
                create_term_fn());

fit_iteratively(@step_until_10neg4, ...
                create_term_fn());

end
