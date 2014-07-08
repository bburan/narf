% function fit05cbounded()
%
% just like fit05 but normalize MSE by SE and stop at 1e-6
% also, this version turns on the 'ScaleStepSize' option
function fit05cbounded()

semse();

MaxStepsPerIteration=10;
StepGrowth=1.1;
DoScaleStepSize = true;

function [a,b,c,d] = step_until_10neg3(prev_opts)              
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepSize', 10^-6,...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg35(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3.5, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepSize', 10^-6,...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg4(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-4, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepSize', 10^-7,...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg45(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-4.5, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepSize', 10^-7,...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg5(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-5, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg55(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-5.5, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg6(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-6, ...
                            'ScaleStepSize', DoScaleStepSize, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

fit_boo('StopAtAbsScoreDelta', 10^-2.0, ...
        'ScaleStepSize', DoScaleStepSize, ...
        'StepGrowth', StepGrowth);

fit_iteratively(@step_until_10neg3, ...
                create_term_fn());

fit_iteratively(@step_until_10neg35, ...
                create_term_fn());

fit_iteratively(@step_until_10neg4, ...
                create_term_fn());

fit_iteratively(@step_until_10neg45, ...
                create_term_fn());

fit_iteratively(@step_until_10neg5, ...
                create_term_fn());

fit_iteratively(@step_until_10neg55, ...
                create_term_fn());

fit_iteratively(@step_until_10neg6, ...
                create_term_fn());
end

