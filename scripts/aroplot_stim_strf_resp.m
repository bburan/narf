function aroplot_stim_strf_resp(batch, cellids, modelnames)
global STACK XXX META MODULES;

for kk = 1:length(cellids)
    cellid = cellids{kk};
for ii = 1:length(modelnames)
    model = modelnames{ii};
    
    % Load the model
    sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(batch) ''];
    sql = [sql ' AND cellid="' cellid '"'];
    sql = [sql ' AND modelname="' model '"'];
    results = mysql(sql);
    if isempty(results)
        stacks{mi, ci} = {};
        metas{mi, ci} = {};
        x0s{mi, ci} = {};
        continue;
    elseif (length(results) > 1)
        error('Multiple DB Results found: for %s/%s\n', cellid, model);
    end    
    load_model(char(results(1).modelpath));
    calc_xxx(1);  
    
    w = [];
    h = [];
    % Take the first channel weights and FIR filter that exist
    for jj = 1:length(STACK)
        m = STACK{jj}{1};
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
        % w = w./repmat(sum(w), size(w,1), 1);  % Normalize wrongly
        
        % Normalize using stimulus power        
        rms = compute_channel_normalizers(XXX{2}, 'stim', false);               
        w = w ./ repmat(rms', 1, size(w,2));
    end
    %strf = w * h;    
    
    % Open a figure and begin plotting
    figure('Name', sprintf('Diagram for %s', modelnames{ii}), ...
            'NumberTitle', 'off', 'Position', [20 50 1500 300]);
    
    sf = XXX{1}.test_set{1};    
    
    % ----------------------------
    % Plot the stimulus
    sp = subplot(1,6,1); 
    imagesc(log(squeeze(XXX{2}.dat.(sf).stim(:,1,:)+10^-2))');    
    set(gca,'YDir','normal');
    xlabel('Time'); 
    ylabel('Stimulus Channel');
        
    % ----------------------------
    % Plot the weights
    subplot(1,6,2);
    % For plotting, renormalize the weights
    mm = sum(abs(w));
    w = w ./ repmat(mm', 1, size(w,1))';
    
    imagesc(w);
    set(gca,'YDir','normal');
    axis image;
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);    
    xlabel('Output Channel'); 
    ylabel('Stimulus Channel');
    
    % ----------------------------
    % Plot the FIR filter
    subplot(1,6,3);
    imagesc(h);
    set(gca,'YDir','normal');
    axis image;
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);    
    xlabel('Time Latency'); 
    ylabel('Input Channel');
        
    % ----------------------------
    % Plot the reconstructed STRF
    subplot(1,6,4);
    if size(w,2) == size(h,1)
        imagesc(w*h);
        axis image;
        set(gca,'YDir','normal');
        ca = caxis;
        lim = max(abs(ca));
        caxis([-lim, +lim]);    
        xlabel('Time Latency'); 
        ylabel('Input Channel');
    end
        
    % ----------------------------
    % Plot the nonlinearity
    subplot(1,6,5);
    sel = [];  
    sel.stim_idx = 1;
    sel.chan_idx = 1;
    sel.stimfile = XXX{1}.test_set{1};  
    STACK{end-2}{1}.auto_plot(sel, STACK(1:end-2), XXX(1:end-2));
	xlabel('Stimulus'); 
    ylabel('Prediction');
    
    % ----------------------------
    % Plot the Final prediction and response
    subplot(1,6,6);  
    %append_module(MODULES.smooth_respavg.mdl(struct('window', [.2 .2 .2 .2 .2]))); % SMOOTH RESPAVG
    %calc_xxx(length(STACK) -1);
    
    rp = XXX{end}.dat.(sf).respavg(:,1,1);
    st = XXX{end}.dat.(sf).stim(:,1,1);
    plot(XXX{end}.dat.(sf).stim_time, rp, 'k-', ...
         XXX{end}.dat.(sf).stim_time, st, 'r-');
    xlabel('Time'); 
    ylabel('Spike Rate');
        
end
end
end

