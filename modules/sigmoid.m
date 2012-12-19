function m = sigmoid(args)
% Applies a gaussian sigmoid to the desired field 
% The curve is actually just the CDF of a Gaussian.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @sigmoid;
m.name = 'sigmoid';
m.fn = @do_sigmoid;
m.pretty_name = 'Gaussian Sigmoid';
m.editable_fields = {'mu', 'sigma', 'amp', ...
                     'input', 'input_time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time = 'ds_stim_time';
m.output = 'stim';
m.mu = 0.5;
m.sigma = 1/2*pi;
m.amp = 1.0;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @plot_output_vs_time;
m.plot_fns{1}.pretty_name = 'Output vs Time';

m.plot_fns{2}.fn = @plot_nonlinearity;
m.plot_fns{2}.pretty_name = 'Nonlinearity';

m.plot_fns{3}.fn = @plot_nonlinearity_and_histogram;
m.plot_fns{3}.pretty_name = 'Nonlinearity + Histogram';

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_sigmoid(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    fn = @(x) 0.5 * (1 + erf((x - mdl.mu) / sqrt(2 * mdl.sigma^2)));
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, N, F] = size(x.dat.(sf).(mdl.input));
        x.dat.(sf).(mdl.output) = zeros(S, N, F);
        for s = 1:S
            for f = 1:F
                x.dat.(sf).(mdl.output)(s,:,f) = ...
                    fn(x.dat.(sf).(mdl.input)(s,:,f));
            end    
        end
    end
end

function plot_output_vs_time(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    do_plot_time_series(stack, xxx, mdl.input_time, mdl.output);
    
end

function plot_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    [bins, centers] = hist(dat.(mdl.input)(stim_idx,:), 50);
    
    xs = linspace(centers(1), centers(end), 200);
    
    fn = @(x) 0.5 * (1 + erf((x - mdl.mu) / sqrt(2 * mdl.sigma^2)));
    
    plot(xs, fn(xs), 'k-');
    axis tight;
end

function plot_nonlinearity_and_histogram(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    dat = x.dat.(sf);
    
    [bins, centers] = hist(dat.(mdl.input)(stim_idx,:), 50);
    xs = linspace(centers(1), centers(end), 200);
    
    fn = @(x) 0.5 * (1 + erf((x - mdl.mu) / sqrt(2 * mdl.sigma^2)));
    
    %hold on;
    [AX, H1, H2] = plotyy(centers, bins, xs, fn(xs), 'bar', 'plot');
    axis tight;
    %hold off;
       
    
end

end