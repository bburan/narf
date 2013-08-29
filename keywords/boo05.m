function boo05()

global STACK;

nmse();

function [o, stepsizes] = boost_on_linear_scaat_on_nonlinear(prev_opts)    
    phi = pack_fittables(STACK);
    
    % Detect whether the fittables are in a FIR block or not    
    for kk = 1:length(STACK)
        if isfield(STACK{kk}{1}, 'fit_fields') && ~isempty(STACK{kk}{1}.fit_fields)
            module_being_fit = STACK{kk}{1}.name;
            break;
        end
    end
    
    if strcmp(module_being_fit, 'fir_filter')
        if ~exist('init_stepsizes', 'var')
            init_stepsizes = 1;
        end        
        [~,o,stepsizes] = fit_boo('StopAtStepSize', 10^-7, ...
                                  'StepSize', init_stepsizes, ...
                                  'StopAtStepNumber', length(phi), ...
                                  'StepRel', false, ...
                                  'Elitism', false, ...
                                  'StepGrowth', 1.3); 
    else
        if ~exist('init_stepsizes', 'var')
            init_stepsizes = ones(1, length(phi));
        end        
    
        [~,o,stepsizes] = fit_scaat('StopAtStepSize', 10^-7, ...
                                    'StopAtStepNumber', 5*length(phi), ...
                                    'init_stepsizes', init_stepsizes);
    end
end

outer_loop_term_fn = create_term_fn('StopAtAbsScoreDelta', 10^-3);

fit_iteratively(@boost_on_linear_scaat_on_nonlinear, ...
                outer_loop_term_fn);

end