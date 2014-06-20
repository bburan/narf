function m = weight_channels(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @weight_channels;
m.name = 'weight_channels';
m.fn = @do_weight_channels;
m.pretty_name = 'Weight Channels';
m.editable_fields = {'weights', 'y_offset', 'force_positive', ...
                    'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.weights = [1]; % Each column weights several channels to produce
m.y_offset = [0]; % A column of y-offsets to be added to each output chan. 
m.force_positive = 0;  % if 1, all values of weight matrix
                       % rectified before projection
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields

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

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_weight_channels(mdl, x)   
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
         sf=fns{ii};
         
         % Check dimensions are OK
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         if ~isequal(C, size(mdl.weights, 1))
            error('Dimensions of (mdl.input) don''t match weights.');
         end
         
         if isfield(mdl,'force_positive') && mdl.force_positive,
             W=abs(mdl.weights);
         else
             W=mdl.weights;
         end
         D = size(W, 2);
         
         % Rewritten to avoid squeeze and to have improved performance
         % Simple matrix multiplies, reshapes, and bsxfuns are fastest
         d = reshape(x.dat.(sf).(mdl.input), T*S, C);
         tmp = d*W;
         tmp2 = reshape(tmp, T,S,D);
         x.dat.(sf).(mdl.output) = bsxfun(@plus, tmp2, reshape(mdl.y_offset, 1,1,D));

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
    wmat=[];
    for ii = 1:length(mdls)
        wts=mdls{ii}.weights';
        if isfield(mdls{ii},'force_positive') && mdls{ii}.force_positive,
            wts=abs(wts);
        end
        ws=sqrt(mean(wts.^2,2));
        ws(ws<eps)=1;
        wts=wts./repmat(ws,[1 size(wts,2)]);
        [w, h] = size(wts');
        wmat=cat(2,wmat,wts);
    end
    imagesc(wmat);
    hold on;
    for ii=1:(length(mdls)-1);
        plot([1 1].*w.*ii+0.5,[0.5 h+0.5],'w-','LineWidth',2);
    end
    hold off;
    set(gca,'YDir','normal');
    
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    axis tight;    
    do_xlabel('Coef Channel Index');
    do_ylabel('Channel Index');
end



end
