function m = fir_filter(args)
% An N-dimensional FIR Filter is created for each input filter dimension. 
%
% The total number of filter coefficients = num_coefs * num_filts^2
% 
% num_filts should always equal the number of preprocessed channels

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @fir_filter;
m.name = 'fir_filter';
m.fn = @do_fir_filtering;
m.pretty_name = 'FIR Filter';
m.editable_fields = {'num_coefs', 'num_filts'};
m.isready_pred = @fir_filter_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stimulus;
m.plot_fns{1}.pretty_name = 'FIR Filters vs Time';
m.plot_fns{2}.fn = @do_plot_fir_coefs;
m.plot_fns{2}.pretty_name = 'FIR Coefficients (Stem)';
m.plot_fns{3}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{3}.pretty_name = 'FIR Coefficients (Heat map)';
m.plot_fns{4}.fn = @do_plot_summed_prediction;
m.plot_fns{4}.pretty_name = 'FIR Prediction';

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_filts = 2;
m.coefs = zeros(m.num_filts, m.num_filts, m.num_coefs);

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_filts m.num_coefs], size(m.coefs))
    m.coefs = zeros(m.num_filts, m.num_filts, m.num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function x = do_fir_filtering(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Apply the FIR filters to corresponding downsampled stimuli to get the model prediction
    % Since it is linear, the prediction is just the sum of the filters
    % We assume that there are no second order terms combining elements of both filters
    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, T, P] = size(x.dat.(sf).ds_stim);
        
        if ~isequal(P, mdl.num_filts)
           error('Dimensions of ds_stim don''t match filter.');
        end
       
        for s = 1:S
            for fir_idx = 1:mdl.num_filts
                preds = zeros(T, P);
                for fir_dim = 1:mdl.num_filts,
                    preds(:, fir_idx) = preds(:, fir_idx) + ...
                        filter(squeeze(mdl.coefs(fir_idx, fir_dim, :)), [1], ...
                               squeeze(x.dat.(sf).ds_stim(s, :, fir_dim)))'; 
                end
                x.dat.(sf).lf_stim(s, fir_idx, :) = sum(preds, 2);
            end
        end
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
         squeeze(dat.lf_stim(stim_idx,filt_idx, :)), ...
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
    
    tmp = [];
    for ii = 1:mdl.num_filts
        for jj = 1:mdl.num_filts
            tmp(ii*mdl.num_filts + jj, :) = mdl.coefs(ii,jj,:);
        end
    end
    
    for filt_idx = 1:(mdl.num_filts)^2
        stem([1:mdl.num_coefs], tmp(filt_idx,:), pickcolor(filt_idx));
    end
    axis tight;
end

function do_plot_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    tmp = [];
    for ii = 1:mdl.num_filts
        for jj = 1:mdl.num_filts
            tmp(ii*mdl.num_filts + jj, :) = mdl.coefs(ii,jj,:);
        end
    end
    imagesc(tmp);
    
    axis tight;
end

function isready = fir_filter_isready(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    if all(isfield(x, {'dat'}))
        isready = true;
        for sf = fieldnames(x.dat)', sf=sf{1};
            isready = isready && ~isfield(x.dat.(sf), 'lf_stim') && ...
                      all(isfield(x.dat.(sf), {'ds_stim', ...
                                               'ds_stim_time',...
                                               'ds_stim_fs'}));
        end
    else
        isready =false;
    end
end

end