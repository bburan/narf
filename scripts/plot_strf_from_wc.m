function plot_strf_from_wc(batch, cellids, modelnames)
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
        
        % I'm not sure I'm doing this normalization right, but it's a good
        % first try.
        if isempty(w)
            w = [1];
        else
            % w = w./repmat(sum(w), size(w,1), 1); 
            % w 
        end
        strf = w * h;
        
    end

nullfn = @(x) 0;

[strfs, ~, ~] = load_model_batch(batch, cellids, modelnames, ...
                       @extract_strf, nullfn);


for ii = 1:length(strfs)
    strf = strfs{ii};
    figure('Name', sprintf('Reconstituted STRF for %s', cellids{1}), ...
            'NumberTitle', 'off', 'Position', [20 50 900 900]);
    imagesc(strf);
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    xlabel('Time Latency'); 
    ylabel('Gammatone channel');
    hl = title(modelnames{ii});
    set(hl,'interpreter','none');
end
end

