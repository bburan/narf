function m = lindeberg_spectral(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @lindeberg_spectral;
m.name = 'lindeberg_spectral';
m.fn = @do_lindeberg_spectral;
m.pretty_name = 'Lindeberg spectral filter';
m.editable_fields = {'bf', 's', 'add_factor', 'norm_factor', 'order', 'num_channels', ...
    'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.bf  = 0.5; % best freq in units of CHANNELS
m.s   = 0.1; % sigma param in units of CHANNELS
% m.s_optim = {'upper', 30, ...
%              'lower', 0};

m.add_factor = 0;
m.norm_factor = 1;
m.num_channels = 10;
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_lindeberg_spectral_as_heatmap;
m.auto_init = @auto_init_lindeberg_spectral;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_lindeberg_spectral_as_heatmap;
m.plot_fns{1}.pretty_name = 'Heat Map';
m.plot_fns{2}.fn = @do_plot_all_default_outputs;
m.plot_fns{2}.pretty_name = 'Output Channels (All)';

% Optimize this module for tree traversal  
m.required = {m.input, m.time};  % Signal dependencies
m.modifies = {m.output};       % These signals are modified

% Overwrite the default module fields with arguments
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

% % % Old non-scaled version
%     function weights = noncausal(mdl)
%         xx = (1:mdl.num_channels) - mdl.bf;
%         if (mdl.order==0)
%             weights = exp(-((xx).^2)/2/mdl.s)/(2*pi*mdl.s);
%         elseif (mdl.order==1)
%             weights = (-xx).*exp(-(xx.^2)/2/mdl.s)/(2*pi * mdl.s^2);
%         elseif (mdl.order==2)
%             weights = ((-xx).^2-mdl.s).*exp(-(xx.^2)/2/mdl.s)/(2*pi * mdl.s^3);
%         elseif (mdl.order==3)
%             weights = (-xx).*(xx.^2-3*mdl.s).*exp(-(xx.^2)/2/mdl.s)/(2*pi * mdl.s^4);
%         else
%             error('Spectral order not handled (must be an integer in [0,3])')
%         end
%         weights = mdl.add_factor + weights * mdl.norm_factor;
%     end

% % % % New scaled version (v1)
%     function weights = noncausal(mdl)
%         s = mdl.s * mdl.num_channels;
%         bf = mdl.bf * mdl.num_channels;
%         xx = (1:mdl.num_channels) - bf;
%         if (mdl.order==0)
%             weights = exp(-((xx).^2)/2/s)/(2*pi*s);
%         elseif (mdl.order==1)
%             weights = (-xx).*exp(-(xx.^2)/2/s)/(2*pi * s^2);
%         elseif (mdl.order==2)
%             weights = ((-xx).^2-s).*exp(-(xx.^2)/2/s)/(2*pi * s^3);
%         elseif (mdl.order==3)
%             weights = (-xx).*(xx.^2-3*s).*exp(-(xx.^2)/2/s)/(2*pi * s^4);
%         else
%             error('Spectral order not handled (must be an integer in [0,3])')
%         end
%         weights = mdl.add_factor + weights * mdl.norm_factor;
%     end


% % % New scaled version (v2)
    function weights = noncausal(mdl)
        s = mdl.s;
        bf = mdl.bf;
        xx = (1:mdl.num_channels)/mdl.num_channels - bf;
        if (mdl.order==0)
            weights = exp(-((xx).^2)/2/s)/(2*pi*s);
        elseif (mdl.order==1)
            weights = (-xx).*exp(-(xx.^2)/2/s)/(2*pi * s^2);
        elseif (mdl.order==2)
            weights = ((-xx).^2-s).*exp(-(xx.^2)/2/s)/(2*pi * s^3);
        elseif (mdl.order==3)
            weights = (-xx).*(xx.^2-3*s).*exp(-(xx.^2)/2/s)/(2*pi * s^4);
        else
            error('Spectral order not handled (must be an integer in [0,3])')
        end
        weights = mdl.add_factor + weights * mdl.norm_factor;
    end




    function mdl = auto_init_lindeberg_spectral(stack, xxx)
        % NOTE: Unlike most plot functions, auto_init functions get a
        % STACK and XXX which do not yet have this module or module's data
        % added to them.
        if ~isfield(m, 'fit_fields')
            return
        end
        
        % Init num_channels to have the right dimensionality
        mdl = m;
        x = xxx{end};
        fns = fieldnames(x.dat);
        sf = fns{1};
        [~, ~, C] = size(x.dat.(sf).(mdl.input));
        mdl.num_channels = C;
    end


    function x = do_lindeberg_spectral(mdl, x, stack, xxx)
        fns = fieldnames(x.dat);
        for ii = 1:length(fns)
            sf=fns{ii};
            
            weights = noncausal(mdl)';
            
            % Check dimensions are OK
            [T, S, C] = size(x.dat.(sf).(mdl.input));
            if ~isequal(C, size(weights, 1))
                error('Dimensions of (mdl.input) don''t match weights.');
            end
            
            % Rewritten to avoid squeeze and to have improved performance
            % Simple matrix multiplies, reshapes, and bsxfuns are fastest
            d = reshape(x.dat.(sf).(mdl.input), T*S, C);
            %             x.dat.(sf).(mdl.output) = d*weights;
            x.dat.(sf).(mdl.output) = reshape(d*weights, T,S,1);
            
            % % not needed anymore: the baseline is included in the spectral filter
            %             tmp = d*weights;
            %             tmp2 = reshape(tmp, T,S,D);
            %             D = size(weights, 2);
            %             x.dat.(sf).(mdl.output) = bsxfun(@plus, tmp2, reshape(mdl.baseline, 1,1,D));
            
        end
        
    end

% ------------------------------------------------------------------------
% Plot methods

    function do_plot_lindeberg_spectral_as_heatmap(sel, stack, xxx)
        mdls = stack{end};
        
        % Find the min and max values so colors are scaled appropriately
        c_max = 0;
        for ii = 1:length(mdls)
            weights = noncausal(mdls{ii})';
            c_max = max(c_max, max(abs(weights)));
        end
        
        % Plot all parameter sets' coefficients. Separate them by white pixels.
%         xpos = 1;
%         wmat=[];
%         for ii = 1:length(mdls)
%             wts=noncausal(mdls{ii})';
%             ws=sqrt(mean(wts.^2,2));
%             ws(ws<eps)=1;
%             wts=wts./repmat(ws,[1 size(wts,2)]);
%             [w, h] = size(wts');
%             wmat=cat(2,wmat,wts);
%         end
%         imagesc(wmat);
        imagesc(noncausal(mdls{1})');
        colorbar('East','YColor','white');
        
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