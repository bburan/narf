function m = bandpass_filter_bank(args)
% An bandpass filter bank module creation function.
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @bandpass_filter_bank;
m.name = 'bandpass_filter_bank';
m.fn = @do_bandpass_filter;
m.pretty_name = 'Bandpass Filter Bank';
m.editable_fields = {'low_freqs', 'high_freqs', 'order', ...
                     'sampfs', 'stop_dB', ...
                     'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.low_freqs = [1 4];   % Bottom frequencies of bandpass filters in kHz
m.high_freqs = [4 8];  % Top frequencies of bandpass filters in kHz
m.order = 4;           % Filter order
m.sampfs = 50000;           
m.stop_dB = 50;        % Ratio of passband/stopband attenuation, if applicable
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';
m.function = @ellip;   % @butter or @ellip for now

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

% Optional fields
m.auto_plot = @do_plot_bandpass_filter_bank_frq_resp;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_bandpass_filter_bank_frq_resp;
m.plot_fns{1}.pretty_name = 'Filter Frequency Response';
m.plot_fns{2}.fn = @do_plot_bandpass_filtered_channels;
m.plot_fns{2}.pretty_name = 'Filtered Channels';
% m.plot_fns{2}.fn = @do_plot_filtered_spectrogram;
% m.plot_fns{2}.pretty_name = 'Input Spectrogram';
% m.plot_fns{2}.fn = @do_plot_filtered_spectrogram;
% m.plot_fns{2}.pretty_name = 'Output Spectrogram';

% Helper function
function [mdl, coefs] = build_coefs(mdl)              
    % Make sure the values are ordered as we expect
    mins = min(abs(mdl.low_freqs), abs(mdl.high_freqs));
    maxs = max(abs(mdl.low_freqs), abs(mdl.high_freqs));
    mdl.low_freqs = mins;
    mdl.high_freqs = maxs;
        
    % Now bound the values (and perhaps swap them again if needed)
    
    % Bound the values that may be fit.
    mdl.high_freqs = min(mdl.high_freqs, 0.99999999*0.001*(mdl.sampfs*0.5)*ones(size(mdl.high_freqs))); % Max is barely less than 1/2 sampfs
    mdl.low_freqs = min(max(mdl.low_freqs, 0.01), 0.9999*mdl.high_freqs); % Minimum is 10 Hz, max is 99.99% of the max value    
    mdl.stop_dB = max(1, min(100, abs(mdl.stop_dB)));   
    
    % Build the coefficients for the filters
    for i = 1:length(mdl.low_freqs);
        if strcmp(func2str(mdl.function), 'ellip')
            [B,A] = ellip(mdl.order, 0.5, mdl.stop_dB, ...
                          [(1000*mdl.low_freqs(i)/mdl.sampfs)*2,...
                           (1000*mdl.high_freqs(i)/mdl.sampfs)*2]);   
        elseif strcmp(func2str(mdl.function), 'butter')            
            [B,A] = butter(mdl.order, [(1000*mdl.low_freqs(i)/mdl.sampfs)*2,...
                                       (1000*mdl.high_freqs(i)/mdl.sampfs)*2]);  
        else
            error('Invalid or unknown filtering function.');
        end
        coefs{i} = {B,A};
    end 
end

% Finally, define the 'methods' of this module, as if it were a class
function x = do_bandpass_filter(mdl, x)
    
    [mdl, coefs] = build_coefs(mdl);
    
    % For each stimulus file 
    for sf = fieldnames(x.dat)', sf = sf{1};
        % Accumulate a matrix to store all the filter responses...
        filtered_x = [];
        for idx = 1:length(mdl.low_freqs)
            tmp = filter(coefs{idx}{1}, coefs{idx}{2}, x.dat.(sf).(mdl.input));
            tmp2 = abs(hilbert(tmp));
            filtered_x = cat(3, filtered_x, tmp2); 
        end
        % Store that matrix in our data structure
        x.dat.(sf).(mdl.output) = filtered_x;
    end      
end

% Plot the filter responses
function do_plot_bandpass_filtered_channels(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', 'Channel [-]');
end

function do_plot_bandpass_filter_bank_frq_resp(sel, stack, xxx)
    mdl = stack{end}{1};
    x = xxx{end};

    [mdl, coefs] = build_coefs(mdl);
    
    N = 10000;
    hold on;
    for filt_idx = 1:length(mdl.low_freqs)
        ww = 0:(pi/N):pi;
        H = freqz(coefs{filt_idx}{1}, coefs{filt_idx}{2}, ww);
        H = abs(H)';
        loglog(ww', H, 'Color', pickcolor(filt_idx));
        setAxisLabelCallback('X', @(f) (f*mdl.sampfs/(pi*2)));
        setAxisLabelCallback('Y', @(f) f);
        axis tight;
        
        % Add info about the CF and BW
        idx = find(H == max(H(:)), 1, 'first');
        CF = N * ww(idx);
        powerband = (H >= H(idx) * 0.5);
        bot = find(powerband, 1, 'first');
        top = find(powerband, 1, 'last');
        BW = log2(ww(top) / ww(bot)); % In octaves   
        hh = text(CF*2*pi/mdl.sampfs, max(H)-0.12*filt_idx, sprintf('          CF: %0.1f Hz, BW: %0.3f Oct\n', CF, BW));
        set(hh, 'color', pickcolor(filt_idx));
    end 
    hold off;
    
    do_xlabel('Frequency');
    do_ylabel('Amplitude');
        
end

% function do_plot_filtered_spectrogram(sel, stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
% 
%     [sf, stim_idx, chan_idx] = get_baphy_plot_controls(stack);
%     [baphy_mod, ~] = find_modules(stack, 'load_stim_resps_from_baphy', true);
%     dat = x.dat.(sf);
% 
%     logfsgram(dat.(mdl.output)(:, stim_idx, chan_idx)', 4048, baphy_mod.raw_stim_fs, [], [], 500, 12);
%     caxis([-20,40]);
% end

end
