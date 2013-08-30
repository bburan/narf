function boo05()

global STACK;

nmse();

function [a,b,c,d] = boost_on_linear_scaat_on_nonlinear(prev_opts)    
    phi = pack_fittables(STACK);
    
    % Detect whether the fittables are in a FIR block or not    
    for kk = 1:length(STACK)
        if isfield(STACK{kk}{1}, 'fit_fields') && ~isempty(STACK{kk}{1}.fit_fields)
            module_being_fit = STACK{kk}{1}.name;
            break;
        end
    end
    
    if strcmp(module_being_fit, 'fir_filter')   
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
    else
        if exist('prev_opts', 'var')
            [a,b,c,d] = fit_scaat(prev_opts);
        else
            [a,b,c,d] = fit_scaat('StopAtStepSize', 10^-7, ...
                                  'StopAtStepNumber', 5*length(phi));
        end
    end
end

fit_iteratively(@boost_on_linear_scaat_on_nonlinear, ...
                create_term_fn('StopAtAbsScoreDelta', 10^-3));

end