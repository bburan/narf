function termcond = fit_twostep()
    global STACK XXX;
    % Rewire the model so the FIR filter connects directly to the output
    [tempfirmod, firidx] = find_module(STACK, 'fir_filter');
    STACK{firidx}.output = 'fittertempstim';
    [tempcorrmod, corridx] = find_module(STACK, 'correlation');
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
    fit_fminlsq();
    
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
    
    % Now we will just fit the NL
    fit_fminlsq();
    
    % Restore the stack phi values
	for ii = 1:length(STACK)
        STACK{ii}.fit_fields = tempff{ii};
    end
    
    termcond = NaN;
end