function m = lindeberg_filter(args)
% March 2014 - lienard
%
% modification of the fir_filter file, to constrain a s
%
% for now, this is mostly a test to try to get to a symmetrical STRF
% -- not a real implementation of Linderberg's family of functions

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @lindeberg_filter;
m.name = 'lindeberg_filter';
m.fn = @do_lindeberg_filtering;
m.pretty_name = 'Tony Lindeberg Filter';
m.editable_fields = {'lincoefs', 'coefs', 'num_coefs', 'num_dims', 'baseline', ...
    'sum_channels', ...
    'input', 'filtered_input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.order_x = 0;
m.order_t = 0;
m.type = 0; % switch between non-causal (type==0) and causal (type==1) kernels
m.num_coefs = 20;
m.num_dims = 2;
m.baseline = 0;
m.sum_channels=1;
m.lincoefs = [9 0.5 1 2 0 1 0];
m.coefs = zeros(m.num_dims, m.num_coefs);
m.input =  'stim';
m.filtered_input = 'stim_filtered';
m.time =   'stim_time';
m.output = 'stim';
m.init_fit_sig = 'respavg';

% Optional fields
m.is_splittable = true;
m.plot_gui_create_fn = @create_chan_selector_gui;
m.auto_plot = @do_plot_lindeberg_coefs_as_heatmap;
m.auto_init = @auto_init_lindeberg_filter;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_lindeberg_coefs_as_heatmap;
m.plot_fns{1}.pretty_name = 'FIR Coefs (Heat map)';
m.plot_fns{2}.fn = @do_plot_fir_coefs;
m.plot_fns{2}.pretty_name = 'FIR Coefficients (Stem)';
m.plot_fns{3}.fn = @do_plot_all_filtered_channels;
m.plot_fns{3}.pretty_name = 'Filtered Channels (All)';
m.plot_fns{4}.fn = @do_plot_single_filtered_channel;
m.plot_fns{4}.pretty_name = 'Filtered Channels (Single)';
m.plot_fns{5}.fn = @do_plot_filter_output;
m.plot_fns{5}.pretty_name = 'FIR Output';
m.plot_fns{6}.fn = @do_plot_lindeberg_coefs_as_heatmap2;
m.plot_fns{6}.pretty_name = 'FIR Coefs (Heat map x2)';

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

%     function h = noncausal(x, t, s, tau, v, order_x, order_t)
%         if (~exist('s','var'))
%             s = 1;
%         end
%         if (~exist('tau','var'))
%             tau = 2;
%         end
%         if (~exist('v','var'))
%             v = 0;
%         end
%         if (~exist('order_x','var'))
%             order_x = 0;
%         end
%         if (~exist('order_t','var'))
%             order_t = 0;
%         end
%
%         if (t<0)
%             h = 0;
%         else
%             if (order_x==0)
%                 term1 = exp(-((x-v*t)^2)/2/s)/(2*pi*s);
%             elseif (order_x==1)
%                 term1 = (t*v-x)*exp(-((x-v*t)^2)/2/s)/(2*pi * s^2);
%             elseif (order_x==2)
%                 term1 = ((t*v-x)^2-s)*exp(-((x-v*t)^2)/2/s)/(2*pi * s^3);
%             elseif (order_x==3)
%                 term1 = (t*v-x)*((x-t*v)^2-3*s)*exp(-((x-v*t)^2)/2/s)/(2*pi * s^4);
%             end
%
%             if (order_t==0)
%                 term2 = tau*exp(-(tau^2)/2/t) / (sqrt(2*pi) * t^(3/2));
%             elseif (order_t==1)
%                 term2 = (tau^3-3*t*tau)*exp(-(tau^2)/2/t) / (2*sqrt(2*pi) * t^(7/2));
%             elseif (order_t==2)
%                 term2 = tau*(15*t^2-10*t*tau^2+tau^4)*exp(-(tau^2)/2/t) / (4*sqrt(2*pi) * t^(11/2));
%             elseif (order_t==3)
%                 term2 = tau*(-105*t^3*tau+105*t^2*tau^3-21*t*tau^5+tau^7)*exp(-(tau^2)/2/t) / (8*sqrt(2*pi) * t^(15/2));
%             end
%
%             h = term1*term2;
%         end
%     end

%     function h = timecausal(x, t, s, tau, v, order_x, order_t)
%         if (~exist('s','var'))
%             s = 1;
%         end
%         if (~exist('tau','var'))
%             tau = 2;
%         end
%         if (~exist('v','var'))
%             v = 0;
%         end
%         if (~exist('order_x','var'))
%             order_x = 0;
%         end
%         if (~exist('order_t','var'))
%             order_t = 0;
%         end
%
%         if (t<0)
%             h = 0;
%         else
%             if (order_x==0)
%                 term1 = exp(-((x-v*t)^2)/2/s)/(2*pi*s);
%             elseif (order_x==1)
%                 term1 = (t*v-x)*exp(-((x-v*t)^2)/2/s)/(2*pi * s^2);
%             elseif (order_x==2)
%                 term1 = ((t*v-x)^2-s)*exp(-((x-v*t)^2)/2/s)/(2*pi * s^3);
%             elseif (order_x==3)
%                 term1 = (t*v-x)*((x-t*v)^2-3*s)*exp(-((x-v*t)^2)/2/s)/(2*pi * s^4);
%             end
%
%             k=4;
%             if (order_t==0)
%                 term2 = t^(k-1)*exp(-t/tau)/(tau^k*factorial(k));
%             elseif (order_t==1)
%                 term2 = tau^(-k-1)*t^(k-2)*( (k-1)*tau-t )/(factorial(k))*exp(-t/tau);
%             elseif (order_t==2)
%                 term2 = tau^(-k-2)*t^(k-3)*( (k^2-3*k+2)*tau^2-2*(k-1)*t*tau+t^2 )/(factorial(k))*exp(-t/tau);
%             end
%
%
%             h = term1*term2;
%         end
%     end


    function h = timecausal_vectorized(x, t, s, tau, v, order_x, order_t)
        if (~exist('s','var'))
            s = 1;
        end
        if (~exist('tau','var'))
            tau = 2;
        end
        if (~exist('v','var'))
            v = 0;
        end
        if (~exist('order_x','var'))
            order_x = 0;
        end
        if (~exist('order_t','var'))
            order_t = 0;
        end
        
        xlen = length(x);
        x = repmat(x,length(t),1)';
        t = repmat(t,xlen,1);
        
        if (order_x==0)
            term1 = exp(-((x-v*t).^2)/2/s)/(2*pi*s);
        elseif (order_x==1)
            term1 = (t*v-x).*exp(-((x-v*t).^2)/2/s)/(2*pi * s^2);
        elseif (order_x==2)
            term1 = ((t*v-x).^2-s).*exp(-((x-v*t).^2)/2/s)/(2*pi * s^3);
        elseif (order_x==3)
            term1 = (t*v-x).*((x-t*v).^2-3*s).*exp(-((x-v*t).^2)/2/s)/(2*pi * s^4);
        end
        
        
        k=4;
        if (order_t==0)
            term2 = t.^(k-1).*exp(-t/tau)/(tau^k*factorial(k));
            %             term2 = tau*exp(-(tau^2)./t/2) ./ (sqrt(2*pi) * t.^(3/2));
        elseif (order_t==1)
            term2 = tau^(-k-1)*t.^(k-2).*( (k-1)*tau-t )/(factorial(k)).*exp(-t/tau);
            %             term2 = (tau^3-3*t*tau).*exp(-(tau^2)./t/2) ./ (2*sqrt(2*pi) * t.^(7/2));
        elseif (order_t==2)
            term2 = tau^(-k-2)*t.^(k-3).*( (k^2-3*k+2)*tau^2-2*(k-1)*t*tau+t.^2 )/(factorial(k)).*exp(-t/tau);
            %             term2 = tau*(15*t.^2-10*t*tau^2+tau^4).*exp(-(tau^2)./t/2) ./ (4*sqrt(2*pi) * t.^(11/2));
        end
        
        h = term1.*term2;
        h(t<0) = 0;
    end



    function h = noncausal_vectorized(x, t, s, tau, v, order_x, order_t)
        if (~exist('s','var'))
            s = 1;
        end
        if (~exist('tau','var'))
            tau = 2;
        end
        if (~exist('v','var'))
            v = 0;
        end
        if (~exist('order_x','var'))
            order_x = 0;
        end
        if (~exist('order_t','var'))
            order_t = 0;
        end
        
        xlen = length(x);
        x = repmat(x,length(t),1)';
        t = repmat(t,xlen,1);
        
        if (order_x==0)
            term1 = exp(-((x-v*t).^2)/2/s)/(2*pi*s);
        elseif (order_x==1)
            term1 = (t*v-x).*exp(-((x-v*t).^2)/2/s)/(2*pi * s^2);
        elseif (order_x==2)
            term1 = ((t*v-x).^2-s).*exp(-((x-v*t).^2)/2/s)/(2*pi * s^3);
        elseif (order_x==3)
            term1 = (t*v-x).*((x-t*v).^2-3*s).*exp(-((x-v*t).^2)/2/s)/(2*pi * s^4);
        end
        
        if (order_t==0)
            term2 = tau*exp(-(tau^2)./t/2) ./ (sqrt(2*pi) * t.^(3/2));
        elseif (order_t==1)
            term2 = (tau^3-3*t*tau).*exp(-(tau^2)./t/2) ./ (2*sqrt(2*pi) * t.^(7/2));
        elseif (order_t==2)
            term2 = tau*(15*t.^2-10*t*tau^2+tau^4).*exp(-(tau^2)./t/2) ./ (4*sqrt(2*pi) * t.^(11/2));
        elseif (order_t==3)
            term2 = tau*(-105*t.^3*tau+105*t.^2*tau^3-21.*t*tau^5+tau^7).*exp(-(tau^2)./t/2) ./ (8*sqrt(2*pi) * t.^(15/2));
        end
        
        h = term1.*term2;
        h(t<0) = 0;
    end



    function mdl = constrain_lincoefs(mdl)
        
        %xshift
        mdl.lincoefs(1) = min(mdl.lincoefs(1), mdl.num_dims);
        mdl.lincoefs(1) = max(mdl.lincoefs(1), 0);
        
        %tshift
        mdl.lincoefs(2) = min(mdl.lincoefs(2), mdl.num_coefs/4);
        mdl.lincoefs(2) = max(mdl.lincoefs(2), 0);
        
        % Here are some MINIMAL constraints
        %s
        mdl.lincoefs(3) = max(mdl.lincoefs(3), 0.01);
        
        % tau
        mdl.lincoefs(4) = max(mdl.lincoefs(4), 0.1);
        
        % v
        mdl.lincoefs(5) = min(mdl.lincoefs(5), 10);
        mdl.lincoefs(5) = max(mdl.lincoefs(5), -10);
        
        % Here are some reasonable constraints... they are maybe too
        % restrictive
%         %s
%         mdl.lincoefs(3) = min(mdl.lincoefs(3), mdl.num_dims/2);
%         mdl.lincoefs(3) = max(mdl.lincoefs(3), 0.01);
%         
%         % tau
%         mdl.lincoefs(4) = min(mdl.lincoefs(4), mdl.num_coefs/3);
%         mdl.lincoefs(4) = max(mdl.lincoefs(4), 0.1);
%         
%         % v
%         mdl.lincoefs(5) = min(mdl.lincoefs(5), 5);
%         mdl.lincoefs(5) = max(mdl.lincoefs(5), -5);
%         
%         % norm_factor
%         mdl.lincoefs(6) = max(mdl.lincoefs(6), -100);
%         mdl.lincoefs(6) = min(mdl.lincoefs(6), 100);
        
        % add_factor is mdl.lincoefs(7) => should we constrain it too?
        
    end


    function coefs = initialize_coefs(num_dims, num_coefs, xshift, tshift, s, tau, v, order_x, order_t, norm_factor, add_factor, type)
        
        if (~exist('norm_factor','var'))
            norm_factor = 1;
        end
        if (~exist('add_factor','var'))
            add_factor = 0;
        end
        
        %         xshift = min(xshift, num_dims); xshift = max(xshift, 0);
        %         tshift = min(tshift, num_coefs/4); tshift = max(tshift, 0);
        %         s = min(s,num_dims/2); s = max(s,0.1);
        %         tau = min(tau,num_coefs/3); tau = max(tau,0.5);
        %         v = min(v,0.75); v = max(v,-0.75);
        %         norm_factor = max(norm_factor,0);
        
        if (type==0)
            coefs = noncausal_vectorized((1:num_dims)-xshift, (1:num_coefs)-tshift, s, tau, v, order_x, order_t);
        elseif (type==1)
            coefs = timecausal_vectorized((1:num_dims)-xshift, (1:num_coefs)-tshift, s, tau, v, order_x, order_t);
        end
        
        %         coefs = coefs / max(max(coefs));
        coefs = add_factor + coefs * norm_factor;
    end

    function mdl = auto_init_lindeberg_filter(stack, xxx)
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
        [~, ~, C] = size(x.dat.(sf).(mdl.input));
        mdl.num_dims = C;
        
        mdl = constrain_lincoefs(mdl);
        
        xshift = mdl.lincoefs(1);
        tshift = mdl.lincoefs(2);
        s = mdl.lincoefs(3);
        tau = mdl.lincoefs(4);
        v = mdl.lincoefs(5);
        norm_factor = mdl.lincoefs(6);
        add_factor = mdl.lincoefs(7);
        if ~isfield(mdl,'type'),
            type = 0;
        else
            type = mdl.type;
        end
        mdl.coefs = initialize_coefs(mdl.num_dims, mdl.num_coefs, xshift, tshift, s, tau, v, mdl.order_x, mdl.order_t, norm_factor, add_factor, type);
    end

    function x = do_lindeberg_filtering(mdl, x, stack, xxx)
        % Apply the FIR filter across every stimfile
        if ~isfield(mdl, 'sum_channels')
            mdl.sum_channels=1;
        end
        fns = fieldnames(x.dat);
        for ii = 1:length(fns)
            sf=fns{ii};
            
            mdl = constrain_lincoefs(mdl);
            
            xshift = mdl.lincoefs(1); %xshift = min(xshift,10);
            tshift = mdl.lincoefs(2); %tshift = min(tshift,10);
            s = mdl.lincoefs(3);
            tau = mdl.lincoefs(4);
            v = mdl.lincoefs(5);
            norm_factor = mdl.lincoefs(6);
            add_factor = mdl.lincoefs(7);
            if ~isfield(mdl,'type'),
                type = 0;
            else
                type = mdl.type;
            end
            coefs = initialize_coefs(mdl.num_dims, mdl.num_coefs, xshift, tshift, s, tau, v, mdl.order_x, mdl.order_t, norm_factor, add_factor, type);
            
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
            if mdl.sum_channels
                x.dat.(sf).(mdl.output) = squeeze(sum(tmp, 3)) + mdl.baseline;
            else
                x.dat.(sf).(mdl.output) = tmp + mdl.baseline;
            end
        end
    end

% ------------------------------------------------------------------------
% Plot methods

    function do_plot_lindeberg_coefs_as_heatmap(sel, stack, xxx)
        mdls = stack{end};
        
        % Find the min and max values so colors are scaled appropriately
        c_max = 0;
        for ii = 1:length(mdls)
            
            mdls{ii} = constrain_lincoefs(mdls{ii});
            
            xshift = mdls{ii}.lincoefs(1); %xshift = min(xshift,10);
            tshift = mdls{ii}.lincoefs(2); %tshift = min(tshift,10);
            s = mdls{ii}.lincoefs(3);
            tau = mdls{ii}.lincoefs(4);
            v = mdls{ii}.lincoefs(5);
            norm_factor = mdls{ii}.lincoefs(6);
            add_factor = mdls{ii}.lincoefs(7);
            if ~isfield(mdls{ii},'type'),
                type = 0;
            else
                type = mdls{ii}.type;
            end
            coefs = initialize_coefs(mdls{ii}.num_dims, mdls{ii}.num_coefs, xshift, tshift, s, tau, v, mdls{ii}.order_x, mdls{ii}.order_t, norm_factor, add_factor, type);
            c_max = max(c_max, max(abs(coefs(:))));
        end
        
        % Plot all parameter sets' coefficients. Separate them by white pixels.
        xpos = 1;
        hold on;
        for ii = 1:length(mdls)
            
            
            mdls{ii} = constrain_lincoefs(mdls{ii});
            
            xshift = mdls{ii}.lincoefs(1); %xshift = min(xshift,10);
            tshift = mdls{ii}.lincoefs(2); %tshift = min(tshift,10);
            s = mdls{ii}.lincoefs(3);
            tau = mdls{ii}.lincoefs(4);
            v = mdls{ii}.lincoefs(5);
            norm_factor = mdls{ii}.lincoefs(6);
            add_factor = mdls{ii}.lincoefs(7);
            if ~isfield(mdls{ii},'type'),
                type = 0;
            else
                type = mdls{ii}.type;
            end
            coefs = initialize_coefs(mdls{ii}.num_dims, mdls{ii}.num_coefs, xshift, tshift, s, tau, v, mdls{ii}.order_x, mdls{ii}.order_t, norm_factor, add_factor, type);
            
            %             coefs = mdls{ii}.coefs;
            %             % JL hackish
            %             for (i=1:mdls{ii}.num_dims)
            %                 for (j=1:mdls{ii}.num_coefs)
            %                     % gauss([1,1],[1,0;0,1],[1,1])
            %                     mdls{ii}.coefs(i,j) = mdls{ii}.lincoefs(4)*gauss([mdls{ii}.lincoefs(1)*mdls{ii}.num_dims,mdls{ii}.lincoefs(2)*mdls{ii}.num_coefs],[mdls{ii}.lincoefs(3),0;0,mdls{ii}.lincoefs(3)],[i,j]);
            %                 ends
            %             end
            %             coefs = mdls{ii}.coefs;
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

    function do_plot_lindeberg_coefs_as_heatmap2(sel, stack, xxx)
        mdls = stack{end};
        for ii = 1:length(mdls)
            mdls{ii}.lincoefs(1) = mdls{ii}.lincoefs(1) *2;
            mdls{ii}.num_dims = mdls{ii}.num_dims*2;
        end
        do_plot_lindeberg_coefs_as_heatmap(sel, mdls, xxx);
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

    function do_plot_lindeberg_coefs(sel, stack, xxx)
        mdls = stack{end};
        
        % Plot all parameter sets' coefficients. Separate them by white pixels
        hold on;
        handles = [];
        names = {};
        n_mdls = length(mdls)
        for ii = 1:n_mdls
            %             mdl = mdls{ii};
            
            mdls{ii} = constrain_lincoefs(mdls{ii});
            
            
            xshift = mdls{ii}.lincoefs(1); %xshift = min(xshift,10);
            tshift = mdls{ii}.lincoefs(2); %tshift = min(tshift,10);
            s = mdls{ii}.lincoefs(3);
            tau = mdls{ii}.lincoefs(4);
            v = mdls{ii}.lincoefs(5);
            norm_factor = mdls{ii}.lincoefs(6);
            add_factor = mdls{ii}.lincoefs(7);
            type = mdls{ii}.type;
            coefs = initialize_coefs(mdls{ii}.num_dims, mdls{ii}.num_coefs, xshift, tshift, s, tau, v, mdls{ii}.order_x, mdls{ii}.order_t, norm_factor, add_factor, type);
            
            %             for (i=1:mdl.num_dims)
            %                 for (j=1:mdl.num_coefs)
            %                     % gauss([1,1],[1,0;0,1],[1,1])
            %                     mdl.coefs(i,j) = mdl.lincoefs(4)*gauss([mdl.lincoefs(1)*mdl.num_dims,mdl.lincoefs(2)*mdl.num_coefs],[mdl.lincoefs(3),0;0,mdl.lincoefs(3)],[i,j]);
            %                 end
            %             end
            %             coefs = mdl.coefs;
            [w, h] = size(coefs);
            for c = 1:w
                handles(end+1) = stem((0.3*(ii/n_mdls))+(1:mdls{ii}.num_coefs), coefs(c, :), 'Color', pickcolor(c), 'LineStyle', pickline(ii));
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