function boo03()

global STACK;

nmse();

function [a,b,c,d] = take_as_many_steps_as_parameters(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_boo(prev_opts);
    else
        [a,b,c,d] = fit_boo('StopAtStepSize', 10^-7, ...
                'InitStepSize', 1.0, ...
                'StopAtStepNumber', length(pack_fittables(STACK)), ...
                'StepRel', false, ...
                'Elitism', false, ...
                'StepGrowth', 1.3);
    end
end


fit_iteratively(@take_as_many_steps_as_parameters, ...
                create_term_fn('StopAtAbsScoreDelta', 10^-1));

end