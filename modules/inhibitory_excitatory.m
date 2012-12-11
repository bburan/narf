function m = inhibitory_excitatory(args)
% A pair of FIR filters with only excitory and inhibitory behavior.

% num_filts should always equal the number of preprocessed channels

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @inhibitory_excitatory;
m.name = 'inhibitory_excitatory';
m.fn = @do_inhib_excitatory_filtering;
m.pretty_name = 'Inhibitory Excitatory Filter';
m.editable_fields = {'num_coefs', 'num_dims', 'inhib_coefs', 'excit_coefs'};
m.isready_pred = @fir_filter_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_excitation;
m.plot_fns{1}.pretty_name = 'Excitation vs Time';
m.plot_fns{2}.fn = @do_plot_inhibition;
m.plot_fns{2}.pretty_name = 'Inhibition vs Time';
m.plot_fns{3}.fn = @do_plot_inhib_coefs;
m.plot_fns{3}.pretty_name = 'Inhibitory Coefficients (Stem)';
m.plot_fns{4}.fn = @do_plot_inhib_coefs_as_heatmap;
m.plot_fns{4}.pretty_name = 'Inhibitory Coefficients (Heat map)';
m.plot_fns{5}.fn = @do_plot_excit_coefs;
m.plot_fns{5}.pretty_name = 'Excitatory Coefficients (Stem)';
m.plot_fns{6}.fn = @do_plot_excit_coefs_as_heatmap;
m.plot_fns{6}.pretty_name = 'Excitatory Coefficients (Heat map)';
m.plot_fns{7}.fn = @do_plot_summed_prediction;
m.plot_fns{7}.pretty_name = 'Summed Prediction';

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_dims = 2;
m.inhib_coefs = zeros(m.num_dims, m.num_coefs);
m.excit_coefs = zeros(m.num_dims, m.num_coefs);

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_dims m.num_coefs], size(m.coefs))
    m.inhib_coefs = zeros(m.num_dims, m.num_coefs);
    m.excit_coefs = zeros(m.num_dims, m.num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function x = do_inhibitory_excitatory_filtering(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, T, P] = size(x.dat.(sf).ds_stim);
        
        if ~isequal(P, mdl.num_dims)
           error(['Dimensions of ds_stim don''t match filters. Check ' ...
                  'that preprocessing channel count matches num_dims.']);
        end
        
        
        % Excitatory filter
        x.dat.(sf).excit_preds = zeros(S, T, P);                
        for s = 1:S
            for fir_dim = 1:mdl.num_dims,
                x.dat.(sf).excit_preds(s, :, fir_dim) = ...
                    filter(squeeze(mdl.excit_coefs(fir_dim, :)), [1], ...
                           squeeze(x.dat.(sf).ds_stim(s, :, fir_dim)))';
            end
        end
        
        % Inhibitory filter
        x.dat.(sf).inhib_preds = zeros(S, T, P);                
        for s = 1:S
            for fir_dim = 1:mdl.num_dims,
                x.dat.(sf).inhib_preds(s, :, fir_dim) = ...
                    filter(squeeze(mdl.inhib_coefs(fir_dim, :)), [1], ...
                           squeeze(x.dat.(sf).ds_stim(s, :, fir_dim)))';
            end
        end        
        
        x.dat.(sf).inhib = sum(x.dat.(sf).inhib_preds, 3);
        x.dat.(sf).excit = sum(x.dat.(sf).excit_preds, 3);
        x.dat.(sf).lf_stim = x.dat.(sf).excit - x.dat.(sf).inhib;
    end
end

function do_plot_inhibition(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    [S, T, P] = size(x.dat.(sf).inhib_preds);
    hold on;
    for p = 1:P
        plot(dat.ds_stim_time, squeeze(dat.inhib_preds(stim_idx,:,p)), pickcolor(p));
    end
    axis tight;
    hold off;
    drawnow;
end

function do_plot_excitation(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    [S, T, P] = size(x.dat.(sf).excit_preds);
    hold on;
    for p = 1:P
        plot(dat.ds_stim_time, squeeze(dat.excit_preds(stim_idx,:,p)), pickcolor(p));
    end
    axis tight;
    hold off;
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
       
    dat = x.dat.(sf);
    
    % Scale the response and prediction in case they have wildly
    % different scales (a common problem when using a correlation
    % coefficient-type performance metric is used to fit the model
    rs = mean(squeeze(dat.raw_respavg(stim_idx, :)));
    ss = mean(squeeze(dat.lf_stim(stim_idx, :)));
    
    hold on;
    plot(dat.raw_resp_time, (1/rs)*dat.raw_respavg(stim_idx, :), 'k-');
    plot(dat.ds_stim_time, (1/ss)*squeeze(dat.lf_stim(stim_idx, :)), 'r-');
    hold off;
   
    axis tight;
    drawnow;
end

function do_plot_inhib_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
        
    hold on;
    for dim_idx = 1:(mdl.num_dims)
        stem([1:mdl.num_coefs], mdl.inhib_coefs(dim_idx,:), pickcolor(dim_idx));
    end
    hold off;
    axis tight;
end

function do_plot_excit_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
        
    hold on;
    for dim_idx = 1:(mdl.num_dims)
        stem([1:mdl.num_coefs], mdl.excit_coefs(dim_idx,:), pickcolor(dim_idx));
    end
    hold off;
    axis tight;
end

function do_plot_inhib_coefs_as_heatmap(stack, xxx)
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
    % tmp = mdl.incoefs;
    % [M, N] = size(tmp);
    % for ii = 1:M
    %     tmp(ii,:) = tmp(ii,:) * abs(mean(squeeze(dat.lf_preds(stim_idx, :, ii))));
    % end
    
    imagesc(mdl.inhib_coefs);
    set(gca,'YDir','normal');

end

function do_plot_excit_coefs_as_heatmap(stack, xxx)
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
    
    imagesc(mdl.excit_coefs);
    set(gca,'YDir','normal');

end

end
