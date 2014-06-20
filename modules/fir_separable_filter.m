function m = fir_separable_filter(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @fir_separable_filter;
m.name = 'fir_separable_filter';
m.fn = @do_fir_separable_filter;
m.pretty_name = 'Separable FIR';
m.editable_fields = {'spec_weights','time_weights', ...
    'num_coefs', 'num_dims', 'v', ...
    'filtered_input', 'sum_channels', ... %'coefs', ...
    'baseline', 'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.spec_weights = 0;
m.time_weights = 0;
m.v=0;
m.num_coefs = 20;
m.num_dims = 2;
m.coefs = zeros(m.num_dims, m.num_coefs);
m.baseline = 0;
m.filtered_input = 'stim_filtered';
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';
m.sum_channels=1;
m.init_fit_sig = 'respavg';
% m.weights_old = m.weights; % Stored copy
% m.y_offset_old = m.y_offset; % Stored copy

% Optional fields
m.is_splittable = true;
m.plot_gui_create_fn = @create_chan_selector_gui;
m.auto_plot = @do_plot_fir_separable_filter_as_heatmap;
m.auto_init = @auto_init_fir_separable_filter;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_fir_separable_filter_as_heatmap;
m.plot_fns{1}.pretty_name = 'Heat Map';
m.plot_fns{2}.fn = @do_plot_all_default_outputs;
m.plot_fns{2}.pretty_name = 'Output Channels (All)';

% Overwrite the default module fields with arguments
if nargin > 0
    m = merge_structs(m, args);
end

% Reset the FIR filter coefficients if its size doesn't match num_coefs
if m.num_dims ~= length(m.spec_weights)
    m.spec_weights = zeros(1, m.num_dims);
end
if m.num_coefs ~= length(m.time_weights)
    m.time_weights = zeros(1, m.num_coefs);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

    function mdl = auto_init_fir_separable_filter(stack, xxx)
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
        
        mdl.time_weights = ones(1, mdl.num_coefs);
        mdl.spec_weights = ones(1, mdl.num_dims);
        mdl.coefs = ones(mdl.num_dims, mdl.num_coefs);
        
%         mdl.time_weights = normrnd(0, 1, [1 mdl.num_coefs]);
%         mdl.spec_weights = normrnd(0, 1, [1 mdl.num_dims]); % random initialization!
%         mdl.coefs = mdl.spec_weights' * mdl.time_weights;
%         [xi, yi] = meshgrid(1:mdl.num_coefs, 1:mdl.num_dims);
%         mdl.coefs = interp2(mdl.coefs, xi, yi-repmat((0:(mdl.num_coefs-1))*v,mdl.num_dims,1), 'linear', 0);
    end

    function x = do_fir_separable_filter(mdl, x, stack, xxx)
        % Apply the FIR filter across every stimfile
        if ~isfield(mdl, 'sum_channels')
            mdl.sum_channels=1;
        end
        fns = fieldnames(x.dat);
        for ii = 1:length(fns)
            sf=fns{ii};
            
            mdl.coefs = mdl.spec_weights' * mdl.time_weights;
            [xi, yi] = meshgrid(1:mdl.num_coefs, 1:mdl.num_dims);
            mdl.coefs = interp2(mdl.coefs, xi, yi-repmat((0:(mdl.num_coefs-1))*mdl.v,mdl.num_dims,1), 'linear', 0);
            
            % Have a space allocated for computing initial filter conditions
            init_data_space = ones(size(mdl.coefs, 2) * 2, 1);
            
            % Compute the size of the filter
            [T, S, C] = size(x.dat.(sf).(mdl.input));
            
            if ~isequal(C, length(mdl.spec_weights))
                error('Dimensions of (mdl.input) don''t match channel count.');
            end
            
            zf = {};
            for c = 1:C
                % Find proper initial conditions for the filter
                [~, zftmp] = filter(mdl.coefs(c,:)', [1], ...
                    init_data_space .* x.dat.(sf).(mdl.input)(1, 1, c));
                zf{c} = zftmp;
            end
            
            tmp = zeros(T*S, C);
            dd = reshape(x.dat.(sf).(mdl.input), T*S, C);
            for c = 1:C,
                tmp(:, c) = filter(mdl.coefs(c,:)', [1], dd(:,c), zf{c}, 1);
            end
            tmp = reshape(tmp, T, S, C);
            x.dat.(sf).(mdl.filtered_input) = tmp;
            
            % The output is the sum of the filtered channels
            if mdl.sum_channels
                x.dat.(sf).(mdl.output) = sum(tmp, 3) + mdl.baseline;
            else
                x.dat.(sf).(mdl.output) = tmp + mdl.baseline;
            end
        end
    end


% ------------------------------------------------------------------------
% Plot methods

    function do_plot_fir_separable_filter_as_heatmap(sel, stack, xxx)
        mdls = stack{end};
        
        % Find the min and max values so colors are scaled appropriately
        c_max = 0;
        for ii = 1:length(mdls)
            coefs = mdls{ii}.spec_weights' * mdls{ii}.time_weights;
            [xi, yi] = meshgrid(1:mdls{ii}.num_coefs, 1:mdls{ii}.num_dims);
            coefs = interp2(coefs, xi, yi-repmat((0:(mdls{ii}.num_coefs-1))*mdls{ii}.v,mdls{ii}.num_dims,1), 'linear', 0);
            c_max = max(c_max, max(abs(coefs(:))));
        end
        
        %  Plot all parameter sets' coefficients. Separate them by white pixels.
        wmat=[];
        
        weight_mdl=find_modules(stack,'weight_channels');
        if ~isempty(weight_mdl),
            weight_mdl=weight_mdl{end};
            
            for ii = 1:length(mdls)
                wts=weight_mdl{ii}.weights;
                sw=std(wts,0,1);
                sw(sw<eps)=1;
                wts=wts./repmat(sw,[size(wts,1) 1]);
                wcoefs = mdls{ii}.spec_weights' * mdls{ii}.time_weights;
                [xi, yi] = meshgrid(1:mdls{ii}.num_coefs, 1:mdls{ii}.num_dims);
                wcoefs = interp2(wcoefs, xi, yi-repmat((0:(mdls{ii}.num_coefs-1))*mdls{ii}.v,mdls{ii}.num_dims,1), 'linear', 0);
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
                coefs = mdls{ii}.spec_weights' * mdls{ii}.time_weights;
                [xi, yi] = meshgrid(1:mdls{ii}.num_coefs, 1:mdls{ii}.num_dims);
                coefs = interp2(coefs, xi, yi-repmat((0:(mdls{ii}.num_coefs-1))*mdls{ii}.v,mdls{ii}.num_dims,1), 'linear', 0);
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
            coefs = mdls{ii}.spec_weights' * mdls{ii}.time_weights;
            [xi, yi] = meshgrid(1:mdls{ii}.num_coefs, 1:mdls{ii}.num_dims);
            coefs = interp2(coefs, xi, yi-repmat((0:(mdls{ii}.num_coefs-1))*mdls{ii}.v,mdls{ii}.num_dims,1), 'linear', 0);
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

end