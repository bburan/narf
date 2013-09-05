function fitSubstack(startidx,stopdelta,output)
    
    global STACK META MODULES;
    
    if ~exist('startidx','var') || isempty(startidx),
        startidx=length(STACK); % default, only fit last module
    end
    if ~exist('stopdelta','var') || isempty(stopdelta),
        stopdelta=10^-2;
    end
    if ~exist('output','var') || isempty(output),
        output='stim';
    end
    
    savefitparms=cell(startidx-1,1);
    for ii=1:(startidx-1),
        if isfield(STACK{ii}{1},'fit_fields'),
            savefitparms{ii}=STACK{ii}{1}.fit_fields;
            STACK{ii}{1}.fit_fields={};
        end
    end
    
    % We must add the MSE module temporarily in a very specific way       
    append_module(MODULES.mean_squared_error.mdl(struct('input1', output)));   
    append_module(MODULES.correlation.mdl(struct('input1', output)));    
    META.perf_metric = @pm_nmse;
    STACK{end}{1}.input1=output;
    
    phi_init = pack_fittables(STACK);
    %fit_boost(length(phi_init)*5,min_stepsize,min_scoredelta,false,true);  
    
    % Force one step to occur no matter what
    fit_boo('StopAtStepNumber', 1, ...
            'StopAtStepSize', 10^-12);
        
    fit_boo('StopAtAbsScoreDelta', stopdelta, ...
            'StopAtStepNumber', length(phi_init)*5, ...
            'StopAtStepSize', 10^-7, ...
            'StepGrowth', 1.3, ...
            'RelStep', true, ...
            'Elitism', false);
        
    for ii=1:(startidx-1),
        if isfield(STACK{ii}{1},'fit_fields'),
            STACK{ii}{1}.fit_fields=savefitparms{ii};
        end
    end
    pop_module();
    pop_module();
    