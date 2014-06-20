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
                     'sum_channels',...
                    'input', 'filtered_input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_dims = 2;
m.baseline = 0;
m.sum_channels=1;
m.coefs = zeros(m.num_dims, m.num_coefs);
m.input =  'stim';
m.filtered_input = 'stim_filtered';
m.time =   'stim_time';
m.output = 'stim';
m.init_fit_sig = 'respavg';

% Optional fields
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

% Optimize this module for tree traversal  
m.required = {m.input, m.time, m.init_fit_sig};  % Signal dependencies
m.modifies = {m.output, m.filtered_input};       % These signals are modified

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
        mdl = m;
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

function x = do_fir_filtering(mdl, x)   
    % Apply the FIR filter across every stimfile
    if ~isfield(mdl, 'sum_channels') 
        mdl.sum_channels=1;
    end
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
         
         zf = {};
         for c = 1:C
            % Find proper initial conditions for the filter
            [~, zftmp] = filter(coefs(c,:)', [1], ...
                     init_data_space .* x.dat.(sf).(mdl.input)(1, 1, c));
            zf{c} = zftmp;
         end
                     
%          % Old way of filtering
%          fprintf('Old Way:');
%          tmp = zeros(T, S, C);     
%          for s = 1:S
%              for c = 1:C,                 
%                  
%                  [tmp(:, s, c), zftmp] = filter(coefs(c,:)', [1], ...
%                      x.dat.(sf).(mdl.input)(:, s, c), ...
%                      zf{c});
%                  zf{c} = zftmp;      
%              end
%          end        

         % Tests revealed that below is 2-5x faster than the above: 
         tmp = zeros(T*S, C);   
         dd = reshape(x.dat.(sf).(mdl.input), T*S, C); % 210000x3   
         for c = 1:C, 
             tmp(:, c) = filter(coefs(c,:)', [1], dd(:,c), zf{c}, 1);
         end                
         tmp = reshape(tmp, T, S, C);
         x.dat.(sf).(mdl.filtered_input) = tmp;
         
         % The output is the sum of the filtered channels
         if mdl.sum_channels
             x.dat.(sf).(mdl.output) = sum(tmp, 3) + mdl.baseline; 
             % tmp2 = sum(tmp, 3);
         else
             x.dat.(sf).(mdl.output) = tmp + mdl.baseline; 
         end
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
    
    %  Plot all parameter sets' coefficients. Separate them by white pixels.
    wmat=[];
    
    weight_mdl=find_modules(stack,'weight_channels');
    if ~isempty(weight_mdl),
        weight_mdl=weight_mdl{end};
        
        for ii = 1:length(mdls)
            wts=weight_mdl{ii}.weights;
            if isfield(weight_mdl{ii},'force_positive') && ...
                    weight_mdl{ii}.force_positive,
                wts=abs(wts);
            end
            sw=std(wts,0,1);
            sw(sw<eps)=1;
            wts=wts./repmat(sw,[size(wts,1) 1]);
            wcoefs = mdls{ii}.coefs;
            if size(wts,2)==size(wcoefs,1),
                coefs = wts*wcoefs;
                wcoefs=wcoefs./max(abs(wcoefs(:))).*max(abs(coefs(:)));
            else
                coefs=[];
            end
            [wc, hc] = size(coefs');
            [w, h] = size([coefs;wcoefs]');
            wmat=cat(2,wmat,[coefs;wcoefs]);
        end
    else
        for ii = 1:length(mdls)
            coefs = mdls{ii}.coefs;
            [wc, hc] = size(coefs');
            [w, h] = size(coefs');
            wmat=cat(2,wmat,coefs);
        end
    end
    imagesc(wmat);
    
    hold on;
    for ii=1:length(mdls);
        if ii<length(mdls),
            plot([1 1].*w.*ii+0.5,[0.5 h+0.5],'w-','LineWidth',2);
        end
        if hc<h,
            plot([0.5 size(wmat,2)+0.5],[hc+0.5 hc+0.5],'w-','LineWidth',2);
        end
        coefs = mdls{ii}.coefs;
        text(w.*ii, 1, sprintf('Sp: %f\nSm: %f', ...
             sparsity_metric(coefs),smoothness_metric(coefs)),...
             'VerticalAlign','bottom','HorizontalAlign','right');
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
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.filtered_input, ...
            sel, 'Time [s]', 'Filtered Channel [-]');
end

function do_plot_single_filtered_channel(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Filtered Channel [-]');
end

function do_plot_filter_output(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
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
