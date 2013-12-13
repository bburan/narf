function fit13()

global STACK;

STACK{2}{1}.fit_fields = {'center_freq_khz', 'Q_factor', ...
                          'N_order', 'delayms'};

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
    
        if strcmp(module_being_fit, 'pz_synapse') || ...
                strcmp(module_being_fit, 'weight_channels')
            if exist('prev_opts', 'var')
                [a,b,c,d] = fit_scaat(prev_opts);
            else
                [a,b,c,d] = fit_scaat('StopAtAbsScoreDelta', del, ...
                                      'InitStepSize', 1, ...
                                      'StopAtStepNumber', 1, ...
                                      'StepGrowth', 1.3);
            end
        else
            if exist('prev_opts', 'var')
                [a,b,c,d] = fit_scaat(prev_opts);
            else
                [a,b,c,d] = fit_scaat('StopAtAbsScoreDelta', del, ...
                                      'InitStepSize', 100, ...
                                      'StopAtStepNumber', 1);
            end
        end
    end
    
    fn = @subfitter;
    
end

fit_boo('StopAtAbsScoreDelta', 10^-2, 'StepGrowth', 1.3);
fit_iteratively(make_subfitter(10^1), create_term_fn('StopAtAbsScoreDelta', 10^1));
fit_iteratively(make_subfitter(10^0.5), create_term_fn('StopAtAbsScoreDelta', 10^0.5));
fit_iteratively(make_subfitter(10^0), create_term_fn('StopAtAbsScoreDelta', 10^0));
fit_iteratively(make_subfitter(10^-0.5), create_term_fn('StopAtAbsScoreDelta', 10^-0.5));
fit_iteratively(make_subfitter(10^-1), create_term_fn('StopAtAbsScoreDelta', 10^-1));
fit_iteratively(make_subfitter(10^-1.5), create_term_fn('StopAtAbsScoreDelta', 10^-1.5));
fit_iteratively(make_subfitter(10^-2), create_term_fn('StopAtAbsScoreDelta', 10^-2));
fit_iteratively(make_subfitter(10^-2.5), create_term_fn('StopAtAbsScoreDelta', 10^-2.5));
fit_iteratively(make_subfitter(10^-3), create_term_fn('StopAtAbsScoreDelta', 10^-3));
fit_iteratively(make_subfitter(10^-3.5), create_term_fn('StopAtAbsScoreDelta', 10^-3.5));
fit_iteratively(make_subfitter(10^-4), create_term_fn('StopAtAbsScoreDelta', 10^-4));
fit_iteratively(make_subfitter(10^-4.5), create_term_fn('StopAtAbsScoreDelta', 10^-4.5));
fit_iteratively(make_subfitter(10^-5), create_term_fn('StopAtAbsScoreDelta', 10^-5));        
fit_iteratively(make_subfitter(10^-5.5), create_term_fn('StopAtAbsScoreDelta', 10^-5.5));            
fit_iteratively(make_subfitter(10^-6), create_term_fn('StopAtAbsScoreDelta', 10^-6));

end