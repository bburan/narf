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
m.editable_fields = {'low_freqs', 'high_freqs', 'order', 'sampfs', 'stop_dB'};
m.plot_fns = {'Freq. Response', @do_plot_elliptic_bandpass_filter_bank_frq_resp};
m.isready_pred = @elliptic_bandpass_isready;

% Module fields that are specific to THIS MODULE
m.low_freqs = [2000 20000];  % Bottom frequencies of bandpass filters
m.high_freqs = [4000 27000]; % Top frequencies of bandpass filters
m.order = 4;                 % What order should the filter be?
m.sampfs = 100000;           % What is the sampling frequency?
m.stop_dB = 50;              % Ratio of passband/stopband attenuation

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

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
function [new_stack, new_x] = do_elliptic_filter(stack, x)
    mdl = stack{end};
    % TODO: Check the structure of x to ensure its validity
    filtered_x=[];
    for idx = 1:length(mdl.low_freqs)
        tmp = filter(mdl.coefs{idx}{1}, mdl.coefs{idx}{2}, x,[],2);
        filtered_x = cat(3, filtered_x, tmp); 
    end
    new_stack = stack;
    new_stack{end+1} = mdl;
end

function do_plot_elliptic_bandpass_filter_bank_frq_resp(mdl)    
    for filt_idx = 1:length(mdl.low_freqs)
        ww = 0:(pi/1000):pi;
        H = freqz(mdl.coefs{filt_idx}{1}, mdl.coefs{filt_idx}{2}, ww);
        loglog(ww, abs(H), pickcolor(filt_idx));
        setAxisLabelCallback(gca, @(f) (f*mdl.sampfs/(3.14*2)), 'X');
        axis tight;
    end 
end

function isready = elliptic_bandpass_isready(mdl, x)
    isready = all(isfield(mdl, {'low_freqs', 'high_freqs', 'order', 'sampfs', 'stop_dB'})) & ...
              isvector(mdl.low_freqs) & ...
              isvector(mdl.high_freqs) & ...
              length(mdl.low_freqs) == length(m.high_freqs);
end

end