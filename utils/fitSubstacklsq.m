function fitSubstacklsq(startidx,stopdelta,output)
    
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
    mods = find_modules(STACK, 'mean_squared_error', true);
    append_module(MODULES.mean_squared_error.mdl(struct('input1', output)));   
    % don't append correlation.  we don't need it, and it just
    % takes cycles.
    
    META.perf_metric = @pm_nmse;
    STACK{end}{1}.input1=output;
    
    phi_init = pack_fittables(STACK);   
          
    fit_scaat('InitStepSize', 1.0, ...
              'StopAtStepsize', 10^-1, ... 
              'StopAtStepNumber', 3, ...
              'StepGrowth', 2.0, ...
              'StepShrink', 0.5);

    %options = optimset('MaxIter', 500, ...
    %                   'MaxFunEvals', 500, ...
    %                   'TolFun', 1e-12, 'TolX', 1e-9);            
    %fit_lsq('stim', 'respavg', options);
    
    for ii=1:(startidx-1),
        if isfield(STACK{ii}{1},'fit_fields'),
            STACK{ii}{1}.fit_fields=savefitparms{ii};
        end
    end
    
    pop_module();  % trim the mean_squared_error module from the stack
    