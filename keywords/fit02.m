function fit02()

nmse();


function [a,b,c,d] = step_until_10neg3(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3, ...
                            'StepGrowth', 1.3);
    end
end

function [a,b,c,d] = step_until_10neg4(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-4, ...
                            'StepGrowth', 1.3);
    end
end
    
function [a,b,c,d] = step_until_10neg5(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-5, ...
                            'StepGrowth', 1.3);
    end
end

fit_boo('StopAtAbsScoreDelta', 10^-2, ...
        'StepGrowth', 1.3);

fit_iteratively(@step_until_10neg3, ...
                create_term_fn());

fit_iteratively(@step_until_10neg4, ...
                create_term_fn());

fit_iteratively(@step_until_10neg5, ...
                create_term_fn());
end