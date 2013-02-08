function m = fir_filter(args)
% A Single N-dimensional FIR that spans the input space. 
%
% The total number of filter coefficients = num_coefs * num_dims
% 
% NUM_DIMS should always equal the size of the 'channel' input dimension.
%
% If THE_XXX is defined, NUM_DIMS and COEFS will be initialized from that,
% regardless of what you put into ARGS. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @fir_filter;
m.name = 'fir_filter';
m.fn = @do_fir_filtering;
m.pretty_name = 'FIR Filter';
m.editable_fields = {'num_coefs', 'num_dims', 'coefs', 'baseline', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_dims = 2;
m.baseline = 0;
m.coefs = zeros(m.num_coefs, m.num_dims);
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';
m.init_fit_sig = 'respavg'; % For initializing coefficients only

% Optional fields
m.auto_init = @auto_init_fir_filter;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{1}.pretty_name = 'FIR Coefficients (Heat map)';
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
if ~isequal([m.num_dims m.num_coefs], size(m.coefs))
    m.coefs = zeros(m.num_dims, m.num_coefs);
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
    
    % Initialize coefs automatically if it's in fit_fields
    if any(strcmp('coefs', mm.fit_fields))
        stim = [];
        resp = [];
        for ii = 1:length(xxx{end}.training_set),
            f = xxx{end}.training_set{ii};
            stim = cat(1, stim, xxx{end}.dat.(f).(mm.input));        
            resp = cat(1, resp, xxx{end}.dat.(f).(mm.init_fit_sig));
        end
        resp = resp(:);
        [Tx,Sx] = size(resp);
        stim=reshape(permute(stim, [1 3 2]), Tx, numel(stim) / Tx);
        params = [];
        params.altcore     = 'xccorefet';  % Either 'cdcore' or 'xccorefet'
        params.maxlag      = mm.num_coefs - 1;
        params.resampcount = 12;
        params.sfscount    = 10;
        params.sfsstep     = 3;
        strf = cellxcdataloaded(stim, resp, params);
        mm.coefs = strf(1).h;
        mm.num_dims = size(mm.coefs, 1); 
    end
end

function x = do_fir_filtering(stack, xxx)
    mdl = stack{end};
    x = xxx{end};   % Unfortunately, this doesn't seem to copy deeply
    
    % Apply the FIR filter across every stimfile
    for sf = fieldnames(x.dat)', sf=sf{1};

        % ---------------------------
        % TODO: PLACEHOLDER FOR ARBITRARY DIMENSION FILTERING         
%         % Build up the input matrix
%         M = x.dat.(sf).(mdl.input);
%         M_selcols = x.dat.(sf).([mdl.input '_selcols']);
%         
%         chan_idx = 1;
%         
%         % The fir filter acts on every channel
%         for chan_idx = 1:mdl.num_dims
%             cols = M_selcols('chan', chan_idx);
%             X = M(:, cols);
%             [T, N] = size(X);
%             tmp = zeros(T, N); 
%             % Apply the filter to every other column, one at a time
%             for d = 1:N,
%                 tmp(:, d) = filter(mdl.coefs(:, chan_idx), [1], X(:, d));
%             end
%             % Store the results back in the data structure.
%             x.dat.(sf).(mdl.output)(:, cols) = tmp; 
%         end
%        % ------------------------
        
        % Compute the size of the filter
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         
         if ~isequal(C, mdl.num_dims)
            error('Dimensions of (mdl.input) don''t match channel count.');
         end

         tmp = zeros(T, S, C);        
         for s = 1:S
             for c = 1:C,
                 % Find proper initial conditions for the filter
                 [~, Zf] = filter(mdl.coefs(c,:)', [1], ...
                                  ones(length(mdl.coefs(c,:)') * 2, 1) .* x.dat.(sf).(mdl.input)(1, s, c));
                              
                 tmp(:, s, c) = filter(mdl.coefs(c,:)', [1], ...
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
    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    
    [T, S, C] = size(xold.dat.(sf).(mdl.input));
    
    hold on;
    for c = 1:C       
        [~, Zf] = filter(mdl.coefs(c,:)', [1], ...
                         ones(length(mdl.coefs(c,:)') * 2, 1) .* xold.dat.(sf).(mdl.input)(1, stim_idx, c));
                              
        plot(xold.dat.(sf).(mdl.time), ...
             filter(mdl.coefs(c, :) , [1], ...
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
    [sf, stim_idx, baphy_chan_idx] = get_baphy_plot_controls(stack);
    chan_idx = popup2num(mdl.plot_gui.selected_chan_popup);
    dat = x.dat.(sf);

    [~, Zf] = filter(mdl.coefs(c,:)', [1], ...
                     ones(length(mdl.coefs(c,:)') * 2, 1) .* xold.dat.(sf).(mdl.input)(1, stim_idx, chan_idx));

    plot(xold.dat.(sf).(mdl.time), ...
         filter(mdl.coefs(chan_idx, :) , [1], ...
                    xold.dat.(sf).(mdl.input)(:, stim_idx, chan_idx)), Zf);
    axis tight;
end

function do_plot_fir_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    hold on;
    for c = 1:(mdl.num_dims)
        stem([1:mdl.num_coefs], mdl.coefs(c, :), pickcolor(c));
    end
    hold off;
    axis tight;
end

function do_plot_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    mm=max(abs(mdl.coefs(:)));
    imagesc(mdl.coefs,[-mm mm]);
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    axis tight;
end

end