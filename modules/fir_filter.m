function m = fir_filter(args)
% A Single N-dimensional FIR that spans the input space. 
%
% DOCUMENTATION TODO
%
% The total number of filter coefficients = num_coefs * num_dims * the
% number of groupings of respfiles
% 
% NUM_DIMS should always equal the size of the 'channel' input dimension.
%

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @fir_filter;
m.name = 'fir_filter';
m.fn = @do_fir_filtering;
m.pretty_name = 'FIR Filter';
m.editable_fields = {'coefs', 'num_coefs', 'num_dims', 'baseline', ...
                     'input', 'filtered_input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_dims = 2;
m.baseline = 0;
m.coefs = zeros(m.num_dims, m.num_coefs);
m.input =  'stim';
m.filtered_input = 'stim_filtered';
m.time =   'stim_time';
m.output = 'stim';
m.init_fit_sig = 'respavg';

% Optional fields
m.is_splittable = true;
m.plot_gui_create_fn = @create_chan_selector_gui;
m.auto_plot = @do_plot_fir_coefs_as_heatmap;
m.auto_init = @auto_init_fir_filter;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{1}.pretty_name = 'FIR Coefs (Heat map)';
m.plot_fns{2}.fn = @do_plot_fir_coefs;
m.plot_fns{2}.pretty_name = 'FIR Coefficients (Stem)';
m.plot_fns{3}.fn = @do_plot_all_filtered_channels;
m.plot_fns{3}.pretty_name = 'Filtered Channels (All)';
m.plot_fns{4}.fn = @do_plot_single_filtered_channel;
m.plot_fns{4}.pretty_name = 'Filtered Channels (Single)';
m.plot_fns{5}.fn = @do_plot_filter_output;
m.plot_fns{5}.pretty_name = 'FIR Output';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_dims m.num_coefs], size(m.coefs))
    m.coefs = zeros(m.num_dims, m.num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function mdl = auto_init_fir_filter(stack, xxx)
    % NOTE: Unlike most plot functions, auto_init functions get a 
    % STACK and XXX which do not yet have this module or module's data
    % added to them.    
    if ~isfield(m, 'fit_fields') 
        return
    end
    
    % Init the coefs to have the right dimensionality
    mdl = m;
    x = xxx{end};
    fns = fieldnames(x.dat);
    sf = fns{1};
    [T, S, C] = size(x.dat.(sf).(mdl.input));               
    
    mdl.num_dims = C; 
    mdl.coefs = zeros(mdl.num_dims, mdl.num_coefs);
end

function x = do_fir_filtering(mdl, x, stack, xxx)   
    % Apply the FIR filter across every stimfile
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
         sf=fns{ii};
         
         coefs = mdl.coefs;
    
         % Have a space allocated for computing initial filter conditions
         init_data_space = ones(size(coefs, 2) * 2, 1);
    
         % Compute the size of the filter
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         
         if ~isequal(C, mdl.num_dims)
            error('Dimensions of (mdl.input) don''t match channel count.');
         end

         tmp = zeros(T, S, C);        
         
         % Filter!
         for s = 1:S
             for c = 1:C,
                 % Find proper initial conditions for the filter
                 [~, Zf] = filter(coefs(c,:)', [1], ...
                     init_data_space .* x.dat.(sf).(mdl.input)(1, s, c));
                 
                 tmp(:, s, c) = filter(coefs(c,:)', [1], ...
                     x.dat.(sf).(mdl.input)(:, s, c), ...
                     Zf);
             end
         end
         
         % Now the input has been filtered.
         x.dat.(sf).(mdl.filtered_input) = tmp;
         
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = squeeze(sum(tmp, 3)) + mdl.baseline; 
    end
end

% ------------------------------------------------------------------------
% Plot methods

function do_plot_fir_coefs_as_heatmap(sel, stack, xxx)
    mdls = stack{end};

    % Find the min and max values so colors are scaled appropriately
    c_max = 0;
    for ii = 1:length(mdls)
        c_max = max(c_max, max(abs(mdls{ii}.coefs(:))));
    end
    
    % Plot all parameter sets' coefficients. Separate them by white pixels.
    xpos = 1;
    hold on;
    for ii = 1:length(mdls)
        coefs = mdls{ii}.coefs;
        [w, h] = size(coefs');
        imagesc([xpos xpos+w-1], [1, h], coefs, [-c_max-eps c_max+eps]);
        text(xpos, 1, sprintf('Sparsity: %f\nSmoothness: %f', ...
             sparsity_metric(coefs), ...
             smoothness_metric(coefs)));
        xpos = xpos + 1 + size(mdls{ii}.coefs, 2);
    end
    hold off;
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    axis tight;    
    do_xlabel('Coef Time Index');
    do_ylabel('Coef Channel Index');
end

function do_plot_all_filtered_channels(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.filtered_input, ...
            sel, 'Time [s]', 'Filtered Channel [-]');
end

function do_plot_single_filtered_channel(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Filtered Channel [-]');
end

function do_plot_filter_output(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'FIR Output [-]');
end

function do_plot_fir_coefs(sel, stack, xxx)
    mdls = stack{end};
    
    % Plot all parameter sets' coefficients. Separate them by white pixels    
    hold on;
    handles = [];
    names = {};
    n_mdls = length(mdls)
    for ii = 1:n_mdls
        mdl = mdls{ii};
        coefs = mdl.coefs;
        [w, h] = size(coefs);
        for c = 1:w
            handles(end+1) = stem((0.3*(ii/n_mdls))+(1:mdl.num_coefs), coefs(c, :), 'Color', pickcolor(c), 'LineStyle', pickline(ii));
            names{end+1} = ['PS' num2str(ii) '/CH' num2str(c)];
        end 
    end
    hold off;
    axis tight;
    do_xlabel('Coef Time Index');
    do_ylabel('Coef Magnitude');
    legend(handles, names{:}, 'Location', 'NorthWest');
end

end