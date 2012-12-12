function m = exponential(args)
% Applies a scaled exponential function to an input.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @exponential;
m.name = 'exponential';
m.fn = @do_exponential;
m.pretty_name = 'Exponential';
m.editable_fields = {'input', 'time', 'output', ...
                     'offset', 'gain'};
m.isready_pred = @isready_general_purpose;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.input_time = 'stim_time'; % TODO: Let's just name these time
m.output = 'stim';
m.offset = 0;
m.gain = 1;

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

function x = do_exponential(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    for sf = fieldnames(x.dat)', sf=sf{1};
        [S, N, F] = size(x.dat.(sf).(mdl.input));
        x.dat.(sf).(mdl.output) = zeros(S, N, F);
        fn = @(z) mdl.offset + exp(mdl.gain * z); % TODO: Refactor to remove the duplication here and below
        
        for s = 1:S
            for f = 1:F
                x.dat.(sf).(mdl.output)(s,:,f) = ...
                   fn(x.dat.(sf).(mdl.input)(s,:,f));
            end    
        end
    end
end


function isready = isready_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    if all(isfield(x, {'dat'}))
        isready = true;
        for sf = fieldnames(x.dat)', sf=sf{1};
            isready = isready && ...
                      all(isfield(x.dat.(sf), {'pp_stim', ...
                                               'raw_stim',...
                                               'raw_stim_time', ...
                                               'raw_stim_fs'}));
        end
    else
        isready = false;
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
    
    fn = @(z) mdl.offset + exp(mdl.gain * z);
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
        
    fn = @(z) mdl.offset + exp(mdl.gain * z);
    [AX, H1, H2] = plotyy(centers, bins, xs, fn(xs), 'bar', 'plot');
    axis tight;
    
end

end