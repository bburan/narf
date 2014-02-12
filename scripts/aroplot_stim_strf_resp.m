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
    strf = w * h;    
    
    % Open a figure and begin plotting
    figure('Name', sprintf('Diagram for %s', modelnames{ii}), ...
            'NumberTitle', 'off', 'Position', [20 50 1500 300]);
    
    sf = XXX{1}.test_set{1};    
    
    % ----------------------------
    % Plot the stimulus
    sp = subplot(1,5,1);
    imagesc(squeeze(XXX{2}.dat.(sf).stim(:,1,:))');    
    set(gca,'YDir','normal');
    xlabel('Time'); 
    ylabel('Stimulus Channel');
        
    % ----------------------------
    % Plot the weights
    subplot(1,5,2);
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
    subplot(1,5,3);
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
    subplot(1,5,4);
    imagesc(strf);
    axis image;
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);    
    xlabel('Time Latency'); 
    ylabel('Input Channel');
    
    % ----------------------------
    % Plot the Final prediction and response
    subplot(1,5,5);  
    append_module(MODULES.smooth_respavg.mdl(struct('window', [1 2 4 2 1]))); % SMOOTH RESPAVG

    rp = XXX{end}.dat.(sf).respavg(:,1,1);
    st = XXX{end}.dat.(sf).stim(:,1,1);
    plot(XXX{end}.dat.(sf).stim_time, rp, 'k-', ...
         XXX{end}.dat.(sf).stim_time, st, 'r-');
    xlabel('Time'); 
    ylabel('Spike Rate');
        
end
end
end

