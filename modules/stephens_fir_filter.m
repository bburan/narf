function m = stephens_fir_filter(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @stephens_fir_filter;
m.name = 'stephens_fir_filter';
m.fn = @do_stephens_fir_filtering;
m.pretty_name = 'Stephen''s FIR Filter';
m.editable_fields = {'num_dims', 'altcore', 'maxlag', 'resampcount', 'sfscount', 'sfsstep'};
m.isready_pred = @fir_filter_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stimulus;
m.plot_fns{1}.pretty_name = 'FIR Filters vs Time';
m.plot_fns{2}.fn = @do_plot_fir_coefs;
m.plot_fns{2}.pretty_name = 'FIR Coefficients (Stem)';
m.plot_fns{3}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{3}.pretty_name = 'FIR Coefficients (Heatmap)';
m.plot_fns{4}.fn = @do_plot_summed_prediction;
m.plot_fns{4}.pretty_name = 'FIR Prediction';

% Module fields that are specific to THIS MODULE
m.num_dims = 2;
m.altcore = 'cdcore';
m.maxlag = 12;
m.resampcount = 12;
m.sfscount = 10;
m.sfsstep = 3;
m.rasterfs = 200;  % TODO: REMOVE ME
m.coefs = zeros(m.num_dims, 10);

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function x = do_stephens_fir_filtering(stack, xxx)
    global STACK;
    mdl = stack{end};
    x = xxx{end};
    
    % Apply the FIR filters to corresponding downsampled stimuli to get the model prediction
    % Since it is linear, the prediction is just the sum of the filters
    % We assume that there are no second order terms combining elements of both filters
     for sf = fieldnames(x.dat)', sf=sf{1};
         [S, T, P] = size(x.dat.(sf).ds_stim);       
         
         x.dat.(sf).lf_preds = zeros(S, T, P);        
         for s = 1:S
             for fir_dim = 1:mdl.num_dims,
                 x.dat.(sf).lf_preds(s, :, fir_dim) = ...
                     filter(squeeze(mdl.coefs(fir_dim, :)), [1], ...
                            squeeze(x.dat.(sf).ds_stim(s, :, fir_dim)))';
             end
         end       
         x.dat.(sf).lf_stim = sum(x.dat.(sf).lf_preds, 3);
     end
    
end

function do_plot_filtered_stimulus(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    filt_idx = get(filt_pop, 'Value');
       
    dat = x.dat.(sf);
    
    plot(dat.ds_stim_time, ...
         squeeze(dat.lf_stim(stim_idx,filt_idx,:)), ...
         pickcolor(filt_idx));
    axis tight;
    drawnow;
end

function do_plot_summed_prediction(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');

    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    filt_idx = get(filt_pop, 'Value');

    dat = x.dat.(sf);

    plot(dat.ds_stim_time, ...
         squeeze(dat.lf_stim(stim_idx, filt_idx, :)), ...
         pickcolor(filt_idx));

%        % Scale the response and prediction in case they have wildly
%        % different scales (a common problem when using a correlation
%         % coefficient-type performance metric is used to fit the model
%         respavg = squeeze(dat.ds_respavg(GS.selected_stim_idx,:));
%         rs = mean(respavg);
%         stim = squeeze(dat.ds_pred(GS.selected_stim_idx,:));
%         ss = mean(stim);
%         % Plot
%         plot(dat.ds_time, (1/rs)*respavg, 'k-', ...
%              dat.ds_time, (1/ss)*stim, 'r-');

    axis tight;
    drawnow;
end

function do_plot_fir_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
        
    hold on;
    for dim_idx = 1:(mdl.num_dims)
        stem([1:mdl.num_coefs], mdl.coefs(dim_idx,:), pickcolor(dim_idx));
    end
    hold off;
    axis tight;
end

function do_plot_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
       
    dat = x.dat.(sf);
    
    %TODO: This way of scaling the image intensity is specific to what is
    %selected and therefore probably wrong in general.
    tmp = mdl.coefs;
    [M, N] = size(tmp);
    for ii = 1:M
        tmp(ii,:) = tmp(ii,:) * abs(mean(squeeze(dat.lf_preds(stim_idx, :, ii))));
    end
    
    imagesc(tmp);
    set(gca,'YDir','normal');
    % axis tight;
end

end