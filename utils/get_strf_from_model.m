function strfs=get_strf_from_model(batch,cellids,modelnames);
   
global META;

    function strf = extract_strf(x)
        global STACK;
        w = [];
        h = [];
        % Take the first channel weights and FIR filter that exist
        for ii = 1:length(STACK)
            m = STACK{ii}{1};
            if strcmp(m.name, 'weight_channels') && isempty(w),
                w = m.weights;
            end
            if strcmp(m.name, 'fir_filter') && isempty(h),
                h = m.coefs;
            end
        end
                
        if isempty(w)
            w = [1];
        else
            % TODO: Proper normalization uses the stimulus power
            % w = w./repmat(sum(w), size(w,1), 1);  % Normalize wrongly
            w = w; % Don't normalize
        end
        strf = w * h;
        
    end

nullfn = @(x) 0;
[strfs, ~, ~] = load_model_batch(batch, cellids, modelnames, ...
                       @extract_strf, nullfn);

end

