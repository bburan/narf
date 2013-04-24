function m = elliptic_bandpass_filter_bank(args)
% An elliptic bandpass filter bank module creation function.
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @elliptic_bandpass_filter_bank;
m.name = 'elliptic_bandpass_filter_bank';
m.fn = @do_elliptic_filter;
m.pretty_name = 'Elliptic Bandpass Filter Bank';
m.editable_fields = {'low_freqs', 'high_freqs', 'order', ...
                     'sampfs', 'stop_dB', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.low_freqs = [1000 4000];   % Bottom frequencies of bandpass filters
m.high_freqs = [4000 8000];  % Top frequencies of bandpass filters
m.order = 4;                 % What order should the filter be?
m.sampfs = 100000;           % TODO: Rename this
m.stop_dB = 50;              % Ratio of passband/stopband attenuation
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim;
m.plot_fns{1}.pretty_name = 'Filtered Stimulus vs Time';
m.plot_fns{2}.fn = @do_plot_filtered_spectrogram;
m.plot_fns{2}.pretty_name = 'Filtered Stimulus Spectrogram';
m.plot_fns{3}.fn = @do_plot_elliptic_bandpass_filter_bank_frq_resp;
m.plot_fns{3}.pretty_name = 'Filter Frequency Response';
m.plot_gui_create_fn = @(h, stk, xx) create_filter_selector_gui(h, stk, xx, length(m.low_freqs));

% Values computed from here on are not directly editable, but are based on
% the above values in a deterministic manner.
m.coefs={}; 
for i = 1:length(m.low_freqs);
    [B,A] = ellip(m.order, 0.5, m.stop_dB, ...
                  [(m.low_freqs(i)/m.sampfs)*2,...
                   (m.high_freqs(i)/m.sampfs)*2]);   
    m.coefs{i} = {B,A};
end 

% Finally, define the 'methods' of this module, as if it were a class
function x = do_elliptic_filter(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % For each stimulus file 
    for sf = fieldnames(x.dat)', sf = sf{1};
        % Accumulate a matrix to store all the filter responses...
        filtered_x = [];
        for idx = 1:length(mdl.low_freqs)
            tmp = filter(mdl.coefs{idx}{1}, mdl.coefs{idx}{2}, x.dat.(sf).(mdl.input));
            filtered_x = cat(4, filtered_x, tmp); 
        end
        % Store that matrix in our data structure
        x.dat.(sf).(mdl.output) = filtered_x;
    end
   
    
end

% Plot the filter responses
function do_plot_filtered_stim(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    filt_idx = get(mdl.plot_gui.selected_filter_popup, 'Value');
    dat = x.dat.(sf);
    
    plot(dat.(mdl.time), ...
         dat.(mdl.output)(:, stim_idx, chan_idx, filt_idx), ...
         pickcolor(filt_idx));
    axis tight;
end

function do_plot_filtered_spectrogram(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
    [baphy_mod, ~] = find_modules(stack, 'load_stim_resps_from_baphy', true);
    dat = x.dat.(sf);
    
    logfsgram(dat.(mdl.output)(:, stim_idx, chan_idx)', 4048, baphy_mod.raw_stim_fs, [], [], 500, 12); 
    caxis([-20,40]);
end


function do_plot_elliptic_bandpass_filter_bank_frq_resp(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    hold on;
    for filt_idx = 1:length(mdl.low_freqs)
        ww = 0:(pi/1000):pi;
        H = freqz(mdl.coefs{filt_idx}{1}, mdl.coefs{filt_idx}{2}, ww);
        loglog(ww, abs(H), pickcolor(filt_idx));
        setAxisLabelCallback('X', @(f) (f*mdl.sampfs/(3.14*2)));
        axis tight;
    end 
    hold off;
end


end