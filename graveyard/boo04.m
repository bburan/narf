function boo04()

global STACK;

nmse();

function [a,b,c,d] = scaat_as_many_steps_as_params(prev_opts)
    if exist('prev_opts', 'var')
        [a,b,c,d] = fit_scaat(prev_opts);
    else
        [a,b,c,d] = fit_scaat('StopAtStepSize', 10^-7, ...
                              'StopAtStepNumber', length(pack_fittables(STACK)));
    end
end

fit_iteratively(@scaat_as_many_steps_as_params, ...
                create_term_fn('StopAtAbsScoreDelta', 10^-2));          
            
end