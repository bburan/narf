function fit16()

global STACK META;

likelihood();

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
    
        if strcmp(module_being_fit, 'fir_filter') || ...
                strcmp(module_being_fit, 'weight_channels')   
            if exist('prev_opts', 'var')
                [a,b,c,d] = fit_boo(prev_opts);
            else
                [a,b,c,d] = fit_boo('StopAtAbsScoreDelta', del, ...
                                    'StopAtStepNumber', 10, ...
                                    'StepGrowth', 1.3);
            end
        else
            if exist('prev_opts', 'var')
                [a,b,c,d] = fit_scaat(prev_opts);
            else
                [a,b,c,d] = fit_scaat('StopAtAbsScoreDelta', del, ...
                                      'StopAtStepNumber', 10);
            end
        end
        phi = pack_fittables(STACK)';
        fprintf('phi=[');
        fprintf('%0.3f ', phi(:));
        fprintf(']\n');
    end
    
    fn = @subfitter;
end

% Initialization: If FIR filters are all zero, initialize them randomly
[~, mod_idxs] = find_modules(STACK, 'fir_filter');
for ii = 1:length(mod_idxs)
    for jj = 1:length(STACK{mod_idxs{ii}})
        if isfield(STACK{mod_idxs{ii}}{jj}, 'fit_fields') && ...
           any(strcmp('coefs', STACK{mod_idxs{ii}}{jj}.fit_fields)) && ...
           all(0 == STACK{mod_idxs{ii}}{jj}.coefs(:))
            STACK{mod_idxs{ii}}{jj}.coefs = normrnd(0, 10^-3, size(STACK{mod_idxs{ii}}{jj}.coefs));
            fprintf('=====> Initializing FIR coefs to random numbers!\n');
        end
    end
end

%fit_boo('StopAtAbsScoreDelta', 10^-2, 'StepGrowth', 1.3);

% Now gradually shrink the stopping criterion
scale = META.perf_metric()/100; % Start at the current performance value
stop_at=10^-12;

while(scale > stop_at)
    fit_iteratively(make_subfitter(scale), create_term_fn('StopAtAbsScoreDelta', scale));
    scale = scale * 0.6; % Very conservative: 0.8. Probably 0.5 or even 0.25 is fine.
end

nmse();

end
