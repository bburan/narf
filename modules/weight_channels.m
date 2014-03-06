function m = weight_channels(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @weight_channels;
m.name = 'weight_channels';
m.fn = @do_weight_channels;
m.pretty_name = 'Weight Channels';
m.editable_fields = {'weights', 'y_offset', 'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.weights = [1]; % Each column weights several channels to produce
m.y_offset = [0]; % A column of y-offsets to be added to each output chan. 
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_wc_weights_as_heatmap;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_wc_weights_as_heatmap;
m.plot_fns{1}.pretty_name = 'Heat Map';
m.plot_fns{2}.fn = @do_plot_all_default_outputs;
m.plot_fns{2}.pretty_name = 'Output Channels (All)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_weight_channels(mdl, x, stack, xxx)   
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
         sf=fns{ii};
         
         % Check dimensions are OK
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         if ~isequal(C, size(mdl.weights, 1))
            error('Dimensions of (mdl.input) don''t match weights.');
         end

         % Sum the weighted values of each stream to create new channels
         tmp = zeros(T, S, size(mdl.weights, 2));
         for s = 1:S
             in = squeeze(x.dat.(sf).(mdl.input)(:, s, :));
             tmp(:,s,:) = bsxfun(@plus, (in * mdl.weights)', mdl.y_offset)';
         end         
         x.dat.(sf).(mdl.output) = tmp;
    end
end

% ------------------------------------------------------------------------
% Plot methods


function do_plot_wc_weights_as_heatmap(sel, stack, xxx)
    mdls = stack{end};

    % Find the min and max values so colors are scaled appropriately
    c_max = 0;
    for ii = 1:length(mdls)
        c_max = max(c_max, max(abs(mdls{ii}.weights(:))));
    end
    
    % Plot all parameter sets' coefficients. Separate them by white pixels.
    xpos = 1;
    hold on;
    for ii = 1:length(mdls)
        weights = mdls{ii}.weights';
        [w, h] = size(weights');
        imagesc([xpos xpos+w-1], [1, h], weights, [-c_max-eps c_max+eps]);
        text(xpos, 1, sprintf('Sparsity: %f\nSmoothness: %f', ...
             sparsity_metric(weights), ...
             smoothness_metric(weights)));
        xpos = xpos + 1 + size(mdls{ii}.weights, 2);
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



end