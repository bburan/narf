function termcond = fit_nlsep(nl_index)
% Fits the FIR filter to all files in 
% but only fits the 

    global STACK XXX;
    % Rewire the model so the FIR filter connects directly to the output
    [tempfirmod, firidx] = find_modules(STACK, 'fir_filter', true);
    [tempcorrmod, corridx] = find_modules(STACK, 'correlation', true);
    STACK{firidx}.output = 'fittertempstim';
    STACK{corridx}.input1 = 'fittertempstim';
    % Cache the phi values, and remove any fits after the fir filter
    tempff = cell(1, length(STACK));
    for ii = 1:length(STACK)
        if isfield(STACK{ii}, 'fit_fields')
            tempff{ii} = STACK{ii}.fit_fields;
        else
            tempff{ii} = {};
        end
        if ii > firidx
            STACK{ii}.fit_fields = {};
        end
    end
    
    % Now we will just fit the FIR
    calc_xxx(2);
    fit_fminlsq('score', 'fittertempstim', 'respavg');
    
    % Wire the NL up again and fit the NL independently    
    STACK{firidx}.output = tempfirmod.output;
    STACK{corridx}.input1 = tempcorrmod.input1;
    for ii = 1:length(STACK)
        if ii <= firidx
            STACK{ii}.fit_fields = {};
        else
            STACK{ii}.fit_fields = tempff{ii};
        end
    end
    STACK{firidx}.output = 'stim';
    STACK{corridx}.input1 = 'stim';
    
    % Adjust the training matrix
    temptrain = XXX{1}.training_set;
    XXX{2}.training_set = {XXX{2}.training_set{nl_index}};
    
    % Now we will just fit the NL
    calc_xxx(2);
    fit_fminlsq();
    
    % Restore the stack phi values
	for ii = 1:length(STACK)
        STACK{ii}.fit_fields = tempff{ii};
    end
        
    XXX{2}.training_set = temptrain;
    termcond = NaN;
end
