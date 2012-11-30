function m = stephens_fir_filter(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @stephens_fir_filter;
m.name = 'stephens_fir_filter';
m.fn = @do_stephens_fir_filtering;
m.pretty_name = 'Stephen''s FIR Filter';
m.editable_fields = {'maxlag', 'resampcount', 'sfscount', 'sfsstep'};
m.isready_pred = @fir_filter_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stimulus;
m.plot_fns{1}.pretty_name = 'FIR Filters vs Time';
m.plot_fns{2}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{2}.pretty_name = 'FIR Coefficients (Raw STRF)';
m.plot_fns{3}.fn = @do_plot_fir_coefs_as_interpolated_heatmap;
m.plot_fns{3}.pretty_name = 'FIR Coefficients (Interpolated STRF)';
% m.plot_fns{4}.fn = @do_plot_summed_prediction;
% m.plot_fns{4}.pretty_name = 'FIR Prediction';

% Module fields that are specific to THIS MODULE
m.maxlag = 12;
m.resampcount = 12;
m.sfscount = 10;
m.sfsstep = 3;
m.rasterfs = 200;  % TODO: REMOVE ME

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function x = do_stephens_fir_filtering(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % choose fit algorithm and set various parameters
    params=[];
    params.altcore='cdcore';     % boosting, aka coordinate descent
    %params.altcore='xccorefet';  % Theunissen et al 2001 algorithm
    
    % min and max time lags in bins.  Min is typically zero
    params.maxlag=[0 mdl.maxlag];
    
    % some other hyperparameters
    params.resampcount = mdl.resampcount;
    params.sfscount = mdl.sfscount;
    params.sfsstep = mdl.sfsstep;

    % Prepare the stimuli and responses by concatenating everything?!
    % TODO: 
    %     for sf = fieldnames(x.dat)', sf=sf{1};      
    %         for s = 1:S
    %             x.dat.(sf).lf_stim(s, 1, :) = 0;
    %         end
    %     end
    
    % convert to PSTH
    %resp=nanmean(resp,2);
    %resp=resp(:);
    fns = fieldnames(x.dat);
    sf = fns{1};

    % TODO: Make an STRF for each input channel instead of just the 1st
    x.strf = cellxcdataloaded(x.dat.(sf).ds_stim(:,:,1), x.dat.(sf).raw_respavg, params);

    % TODO: Using the STRF, make predictions across every stimulus
    
    x.dat.(sf).lf_stim(s, 1, :) = 0; 
    
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

% function do_plot_summed_prediction(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     
%     % Find the GUI controls
%     baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
%     filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
%     
%     c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
%     sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
%     stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
%     filt_idx = get(filt_pop, 'Value');
% 
%     dat = x.dat.(sf);
%     
%     plot(dat.ds_stim_time, ...
%          squeeze(dat.lf_stim(stim_idx, filt_idx, :)), ...
%          pickcolor(filt_idx));
%      
% %        % Scale the response and prediction in case they have wildly
% %        % different scales (a common problem when using a correlation
% %         % coefficient-type performance metric is used to fit the model
% %         respavg = squeeze(dat.ds_respavg(GS.selected_stim_idx,:));
% %         rs = mean(respavg);
% %         stim = squeeze(dat.ds_pred(GS.selected_stim_idx,:));
% %         ss = mean(stim);
% %         % Plot 
% %         plot(dat.ds_time, (1/rs)*respavg, 'k-', ...
% %              dat.ds_time, (1/ss)*stim, 'r-');
% 
%     axis tight;
%     drawnow;
% end

function do_plot_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    plotastrf(x.strf(1).h, 0, x.stimparam.ff, mdl.rasterfs);
    
    axis tight;
end

function do_plot_fir_coefs_as_interpolated_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    plotastrf(mdl.strf(1).h, 2, mdl.stimparam.ff, mdl.rasterfs);
    
    axis tight;
end


end