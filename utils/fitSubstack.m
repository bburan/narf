function fitSubstack(startidx,stopdelta,output)
    
    global STACK
    
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
    
    % temporarily add mse module
    nmse();
    STACK{end}{1}.input1=output;
    
    min_stepsize = 10^-7; % stop boosting if abs step size smaller than this
    min_scoredelta = stopdelta; % stop boosting if mse improvement smaller than this
    
    phi_init = pack_fittables(STACK);
    fit_boost(length(phi_init)*5,min_stepsize,min_scoredelta,false,true);
    
    for ii=1:(startidx-1),
        if isfield(STACK{ii}{1},'fit_fields'),
            STACK{ii}{1}.fit_fields=savefitparms{ii};
        end
    end
    pop_module();
    pop_module();
    