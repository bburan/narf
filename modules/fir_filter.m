function m = fir_filter(args)
% A Single N-dimensional FIR that spans the input space. 
%
% The total number of filter coefficients = num_coefs * num_dims
% 
% num_dims should always equal the size of the 'chan' dimension.
%
% Module fields that must ALWAYS be defined
m = [];
m.mdl = @fir_filter;
m.name = 'fir_filter';
m.fn = @do_fir_filtering;
m.pretty_name = 'FIR Filter';
m.editable_fields = {'num_coefs', 'num_dims', 'coefs', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_dims = 2;
m.coefs = zeros(m.num_coefs, m.num_dims);
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_;
m.plot_fns{1}.pretty_name = 'FIR Response vs Time';
m.plot_fns{2}.fn = @do_plot_filtered_channels;
m.plot_fns{2}.pretty_name = 'Filtered Channels vs Time';
m.plot_fns{3}.fn = @do_plot_fir_coefs;
m.plot_fns{3}.pretty_name = 'FIR Coefficients (Stem)';
m.plot_fns{4}.fn = @do_plot_fir_coefs_as_heatmap;
m.plot_fns{4}.pretty_name = 'FIR Coefficients (Heat map)';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_dims m.num_coefs], size(m.coefs))
    m.coefs = zeros(m.num_dims, m.num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
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
                 tmp(:, s, c) = filter(mdl.coefs(:, c), [1], ...
                                       x.dat.(sf).(mdl.input)(:, s, c));
             end
         end
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = squeeze(sum(tmp, 3)); 
    end
end

function do_plot_filter_output(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);
    
    [T, S, C] = size(x.dat.(sf).(mdl.output));
    
    plot(dat.(mdl.time), squeeze(dat.(mdl.output)(:, stim_idx)), 'k-');
    axis tight;
end

function do_plot_filtered_channels(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);
    
    [T, S, C] = size(x.dat.(sf).(mdl.input));
    
    hold on;
    for c = 1:C
        plot(dat.(mdl.time), ...
             filter(mdl.coefs(:,c) , [1], dat.(mdl.input)(:, stim_idx, c)), ...
             pickcolor(c));
    end
    axis tight;
    hold off;
end

function do_plot_fir_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    hold on;
    for c = 1:(mdl.num_dims)
        stem([1:mdl.num_coefs], mdl.coefs(:, c), pickcolor(c));
    end
    hold off;
    axis tight;
end

function do_plot_fir_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % TODO: This way of scaling the image intensity is specific to what is
    % selected and therefore probably wrong in general. Replace with proper
    % filtering
    %     tmp = mdl.coefs;
    %     [M, N] = size(tmp);
    %     for ii = 1:M
    %         tmp(ii,:) = tmp(ii,:) * abs(mean(squeeze(dat.(mdl.output)(stim_idx, :, ii))));
    %     end
    
    imagesc(mdl.coefs);
    set(gca,'YDir','normal');
    axis tight;
end

end