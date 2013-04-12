function m = fir_filter(args)
% A Single N-dimensional FIR that spans the input space. 
%
% Updates for version 2:
%     - You are expected to initialize the FIR coefficients with a fitter.
%     - Don't access the coefs structure yourself, please!
%     - Instead, use these methods to read or write them
%          .get_coefs(respfile)
%          .set_coefs(respfile)
%     - By defining the above two methods as desired, you can implement more
%       sophisticated partitioning. For example, you could use one FIR for
%       all the respfiles, one FIR for each type of respfile, or one FIR
%       for each different respfile.
%     - The default setter and getter use one FIR for all respfiles.
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
m.editable_fields = {'num_coefs', 'num_dims', 'baseline', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_dims = 2;
m.baseline = 0;
m.coefs = zeros(m.num_dims, m.num_coefs);
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.auto_plot = @do_plot_fir_coefs_as_heatmap;
m.auto_init = @auto_init_fir_filter;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{1}.pretty_name = 'FIR Coefficients (Heat map)';
m.plot_fns{1}.fn = @do_plot_all_fir_coefs_as_heatmap;
m.plot_fns{1}.pretty_name = 'All FIR Coefficients (Heat map)';
m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, m.time, m.output);
m.plot_fns{2}.pretty_name = 'FIR Response vs Time';
m.plot_fns{3}.fn = @do_plot_all_filtered_channels;
m.plot_fns{3}.pretty_name = 'All Filtered Channels';
m.plot_fns{4}.fn = @do_plot_single_filtered_channel;
m.plot_fns{4}.pretty_name = 'Single Filtered Channel';
m.plot_fns{5}.fn = @do_plot_fir_coefs;
m.plot_fns{5}.pretty_name = 'FIR Coefficients (Stem)';
m.plot_gui_create_fn = @(hh, stack, xxx) create_chan_selector_gui(hh, stack, xxx(1:end-1), m.input);

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_dims m.num_coefs], size(m.fircoefs))
    m.fircoefs = zeros(m.num_dims, m.num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function mm = auto_init_fir_filter(stack, xxx)
    % NOTE: Unlike most plot functions, auto_init functions get a 
    % STACK and XXX which do not yet have this module or module's data
    % added to them. 
    mm = m;
    
    if ~isfield(mm, 'fit_fields') 
        return
    end
    
    % Count the number of elements in the training set
    sfs = fieldnames(xxx{end}.dat);
    n_sfs = length(sfs);
    
    % TODO: Loop through the XXX data structure RESPFILES, passing each of
    % them to get_coefs() and monitor how many elements there are in that
    % set?
    
    % TODO: Init the fircoefs to be the right dimensionality
    
    % TODO: Make all other fields 'match' the previous module output
    
end

function x = do_fir_filtering(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf=fns{ii};
        
        % Figure out the particular parameter set to use
        params = mdl.get_parameter_set(sf);
        mmm = mdl.mdl(params);
        
        
        % Now use mmm instead of mdl to compute crap
        
        
        % 
        
        
    end       
    
    % Old way below the line
    % -----------------------------------------------
    
    % Apply the FIR filter across every stimfile
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
         sf=fns{ii};
         
         coefs = mdl.get_coefs(mdl, sf);
    
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
         
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = squeeze(sum(tmp, 3)) + mdl.baseline; 
    end
end

function do_plot_all_filtered_channels(stack, xxx)
    mdl = stack{end};
    xold = xxx{end-1}; % To print the inputs, you need to go up one
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
    
    [T, S, C] = size(xold.dat.(sf).(mdl.input));
    
    coefs = mdl.get_coefs(mdl, sf);
    
    hold on;
    for c = 1:C       
        [~, Zf] = filter(coefs(c,:)', [1], ...
                         ones(length(coefs(c,:)') * 2, 1) .* xold.dat.(sf).(mdl.input)(1, stim_idx, c));
                              
        plot(xold.dat.(sf).(mdl.time), ...
             filter(coefs(c, :), [1], ...
                    xold.dat.(sf).(mdl.input)(:, stim_idx, c), Zf), ...
             pickcolor(c));
    end
    axis tight;
    hold off;
end

function do_plot_single_filtered_channel(stack, xxx)
    mdl = stack{end};
    xold = xxx{end-1}; % To print the inputs, you need to go up one
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
    chan_idx = popup2num(mdl.plot_gui.selected_chan_popup);
    dat = x.dat.(sf);

    coefs = mdl.get_coefs(mdl, sf);

    [~, Zf] = filter(coefs(chan_idx,:)', [1], ...
                     ones(length(coefs(chan_idx,:)') * 2, 1) .* xold.dat.(sf).(mdl.input)(1, stim_idx, chan_idx));

    plot(xold.dat.(sf).(mdl.time), ...
         filter(coefs(chan_idx, :), [1], ...
                    xold.dat.(sf).(mdl.input)(:, stim_idx, chan_idx), Zf), ...
          pickcolor(chan_idx));
    axis tight;
end

function do_plot_fir_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    hold on;
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    coefs = mdl.get_coefs(mdl, sf);
    for c = 1:(mdl.num_dims)
        stem(1:mdl.num_coefs, coefs(c, :), pickcolor(c));
    end
    hold off;
    axis tight;
end

function do_plot_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    coefs = mdl.get_coefs(mdl, sf);
    
    mm=max(abs(coefs(:)));
    imagesc(coefs,[-mm mm]);
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    axis tight;
    textLoc(sprintf('Sparsity: %f\nSmoothness: %f', ...
        sparsity_metric(coefs), smoothness_metric(coefs)), 'NorthWest');
    
end

function do_plot_all_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    mm=max(abs(mdl.fircoefs(:)));
    imagesc(mdl.fircoefs,[-mm mm]);
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    axis tight;
    textLoc(sprintf('Sparsity: %f\nSmoothness: %f', ...
        sparsity_metric(mdl.fircoefs), smoothness_metric(mdl.fircoefs)), 'NorthWest');
        
end

end