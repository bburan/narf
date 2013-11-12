function fit04()

global STACK;

nmse();

function fn = make_subfitter(del)
    function [a,b,c,d] = subfitter(prev_opts)    
    
        % Detect whether the fittables are in a FIR block or not    
        module_being_fit = '';
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
                [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', del, ...
                                    'StopAtStepNumber', 1, ...
                                    'StepGrowth', 1.3);
            end
        else
            if exist('prev_opts', 'var')
                [a,b,c,d] = fit_scaat(prev_opts);
            else
                [a,b,c,d] = fit_scaat('StopAtAbsScoreDelta', del);
            end
        end
    end
    
    fn = @subfitter;
    
end

fit_boo('StopAtAbsScoreDelta', 10^-2, 'StepGrowth', 1.3);
fit_iteratively(make_subfitter(10^1), create_term_fn('StopAtAbsScoreDelta', 10^1));
fit_iteratively(make_subfitter(10^0), create_term_fn('StopAtAbsScoreDelta', 10^0));
fit_iteratively(make_subfitter(10^-1), create_term_fn('StopAtAbsScoreDelta', 10^-1));
fit_iteratively(make_subfitter(10^-2), create_term_fn('StopAtAbsScoreDelta', 10^-2));
fit_iteratively(make_subfitter(10^-3), create_term_fn('StopAtAbsScoreDelta', 10^-3));
fit_iteratively(make_subfitter(10^-4), create_term_fn('StopAtAbsScoreDelta', 10^-4));
fit_iteratively(make_subfitter(10^-5), create_term_fn('StopAtAbsScoreDelta', 10^-5));           
%fit_iteratively(make_subfitter(10^-6), create_term_fn('StopAtAbsScoreDelta', 10^-6));
%fit_iteratively(make_subfitter(10^-7), create_term_fn('StopAtAbsScoreDelta', 10^-7));

end