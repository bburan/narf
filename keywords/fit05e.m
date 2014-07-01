% function fit05a()
%
% just like fit05 semse instead of nmse, still stopping at 1e-5
%
function fit05e()

    global STACK
    
semse();

MaxStepsPerIteration=10;
StepGrowth=1.1;

function [a,b,c,d] = step_until_10neg3(prev_opts)              
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-3, ...
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
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg55(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-5.5, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

function [a,b,c,d] = step_until_10neg6(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', 10^-6, ...
                            'StopAtStepNumber', MaxStepsPerIteration, ...
                            'StepGrowth', StepGrowth);
    end
end

disp('Removing static spike NL for quick fit');
[nl_mods,nl_idx] = find_modules(STACK, 'nonlinearity', false);
nl_save=nl_mods{end};
nl_idx=nl_idx{end};

for jj=1:length(STACK{nl_idx}),
    STACK{nl_idx}{jj}.nlfn=@nl_dummy;
    STACK{nl_idx}{jj}=rmfield(STACK{nl_idx}{jj},'fit_fields');
end

% now run sequence of progressively tighter fits on model without
% spike NL, and don't run quite so deep.
fit_boo('StopAtAbsScoreDelta', 10^-2.0, ...
        'StepGrowth', StepGrowth);

fit_iteratively(@step_until_10neg3, ...
                create_term_fn());

fit_iteratively(@step_until_10neg35, ...
                create_term_fn());

fit_iteratively(@step_until_10neg4, ...
                create_term_fn());

fit_iteratively(@step_until_10neg45, ...
                create_term_fn());

disp('Restoring static spike NL for final fit');
%pop_module();  % remove MSE module
%update_xxx(2);
STACK{nl_idx}=nl_save;
%fitSubstack([],10^-3.0);
%semse();       % restore MSE module
update_xxx(2);

%
% now run through the series of fits again with static NL included.
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

%fit_iteratively(@step_until_10neg55, ...
%                create_term_fn());

%fit_iteratively(@step_until_10neg6, ...
%                create_term_fn());
end

