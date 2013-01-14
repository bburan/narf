function m = exp_filter(args)
% A Single N-dimensional FIR that spans the input space constrained
% to have an exponetial decay timecourse
%
% The total number of filter coefficients = 3 * num_dims
%  3 parameters are lag, amplitude, and decay time const in units
%  of time bins
% 
% num_dims should always equal the size of the 'chan' dimension.
%
% Module fields that must ALWAYS be defined
m = [];
m.mdl = @exp_filter;
m.name = 'exp_filter';
m.fn = @do_exp_filtering;
m.pretty_name = 'EXP Filter';
m.editable_fields = {'num_coefs', 'num_dims', 'coefs', 'baseline', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_coefs = 20;
m.num_parms = 3;
m.num_dims = 2;
m.coefs = repmat([1.5 0.1 1.1],[m.num_dims 1]);
m.baseline = 0;
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_exp_coefs_as_heatmap;
m.plot_fns{1}.pretty_name = 'EXP Coefficients (Heat map)';
m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, m.time, m.output);
m.plot_fns{2}.pretty_name = 'EXP Response vs Time';
m.plot_fns{3}.fn = @do_plot_all_filtered_channels;
m.plot_fns{3}.pretty_name = 'All Filtered Channels';
m.plot_fns{4}.fn = @do_plot_single_filtered_channel;
m.plot_fns{4}.pretty_name = 'Single Filtered Channel';
m.plot_fns{5}.fn = @do_plot_exp_coefs;
m.plot_fns{5}.pretty_name = 'EXP Coefficients (Stem)';
m.plot_gui_create_fn = @(hh, stack, xxx) create_chan_selector_gui(hh, stack, xxx(1:end-1), m.input);

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

% Reset the EXP filter coefficients if its size doesn't match num_coefs
if ~isequal([m.num_dims m.num_parms], size(m.coefs))
    m.coefs = repmat([1.5 0.1 1.5],[m.num_dims 1]);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS
function f=gen_exp_filter(parms,num_coefs);
    
    latency=parms(1);
    amplitude=parms(2);
    tau=parms(3);
    
    x=(0:(num_coefs-1))';
    
    %if latency<0 || tau<1,
    %    f=zeros(size(x));
    %    return;
    %end
    
    f=amplitude.*exp(-(x-latency)./tau);
    f(x<latency-1)=0;
    
    ff=max(find(x<latency));
    if ~isempty(ff),
        f(ff)=f(ff)*(1-latency+x(ff));
    end
end

function x = do_exp_filtering(stack, xxx)
    mdl = stack{end};
    x = xxx{end};   % Unfortunately, this doesn't seem to copy deeply
    
    % Apply the EXP filter across every stimfile
    for sf = fieldnames(x.dat)', sf=sf{1};

        % Compute the size of the filter
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        
        if ~isequal(C, mdl.num_dims)
            error('Dimensions of (mdl.input) don''t match channel count.');
        end
        
        tmp = zeros(T, S, C);
        for s = 1:S
            for c = 1:C,
                f=gen_exp_filter(mdl.coefs(c,:),mdl.num_coefs);
                tmp(:, s, c) = filter(f, [1], ...
                                      x.dat.(sf).(mdl.input)(:, s, c));
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
        f=gen_exp_filter(mdl.coefs(c,:),mdl.num_coefs);
        plot(xold.dat.(sf).(mdl.time), ...
             filter(f, [1], xold.dat.(sf).(mdl.input)(:, stim_idx, c)), ...
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
    
    f=gen_exp_filter(mdl.coefs(chan_idx,:),mdl.num_coefs);
    plot(xold.dat.(sf).(mdl.time), ...
         filter(f, [1], xold.dat.(sf).(mdl.input)(:, stim_idx, chan_idx)));
    axis tight;
end

function do_plot_exp_coefs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    hold on;
    for c = 1:(mdl.num_dims)
        f=gen_exp_filter(mdl.coefs(c,:),mdl.num_coefs);
        stem([1:mdl.num_coefs], f, pickcolor(c));
    end
    hold off;
    axis tight;
end

function do_plot_exp_coefs_as_heatmap(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    fmatrix=zeros(mdl.num_coefs,mdl.num_dims);
    for c = 1:(mdl.num_dims)
        fmatrix(:,c)=gen_exp_filter(mdl.coefs(c,:),mdl.num_coefs);
    end
    mm=max(abs(fmatrix(:)));
    if mm
        imagesc(fmatrix',[-mm mm]);
    else
        imagesc(fmatrix');
    end
    set(gca,'YDir','normal');
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    axis tight;
end

end