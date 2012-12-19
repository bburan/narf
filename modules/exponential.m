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
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time = 'stim_time';
m.output = 'stim';
m.offset = 0;
m.gain = 1;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_output_vs_time(stack, xxx, m.time, m.output);
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
        y = zeros(S, N, F);
        fn = @(z) exp((mdl.gain * z) - mdl.offset);  % TODO: Remove spot violation here and below
        for s = 1:S
            for f = 1:F
                y(s,:,f) = fn(x.dat.(sf).(mdl.input)(s,:,f));
            end    
        end
    end
    
    x.dat.(sf).(mdl.output) = y;
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

function plot_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    sf = popup2str(baphy_mod.plot_gui.selected_stimfile_popup);
    stim_idx = popup2num(baphy_mod.plot_gui.selected_stim_idx_popup);
    dat = x.dat.(sf);   
    
    [bins, centers] = hist(dat.(mdl.input)(:, stim_idx), 50);
    
    xs = linspace(centers(1), centers(end), 200);
    
    fn = @(z) exp((mdl.gain * z) - mdl.offset); 
    
    plot(xs, fn(xs), 'k-');
    axis tight;
end

function plot_nonlinearity_and_histogram(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    sf = popup2str(baphy_mod.plot_gui.selected_stimfile_popup);
    stim_idx = popup2num(baphy_mod.plot_gui.selected_stim_idx_popup);
    dat = x.dat.(sf);
    
    
    
    [bins, centers] = hist(dat.(mdl.input)(:, stim_idx), 50);
    xs = linspace(centers(1), centers(end), 200);
        
    fn = @(z) exp((mdl.gain * z) - mdl.offset); 
    
    [AX, H1, H2] = plotyy(centers, bins, xs, fn(xs), 'bar', 'plot');
    axis tight;
    
end

end