function m = gammatone_filter(args)
% A module for a small number of gammatones.
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @gammatone_filter;
m.name = 'gammatone_filter';
m.fn = @do_gammatone_filter;
m.pretty_name = 'Gammatone Filter';
m.editable_fields = {'center_freq_khz', 'bandwidth', 'align', 'use_env',... 
                     'input', 'time', 'output' };
m.isready_pred = @isready_always;

m.center_freq_khz = 2; % In kHz. Can be a vector.
m.bandwidth       = 1; % in octaves. Can be a vector.
m.sampfs          = 50000;   
m.align           = false;
m.use_env         = true;
m.input           = 'stim';
m.time            = 'stim_time';
m.output          = 'stim';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_filtered_stim; 
m.plot_fns{1}.pretty_name = 'Output Channels vs Time';

function x = do_gammatone_filter(mdl, x, stack, xxx)

    if length(mdl.center_freq_khz) ~= length(mdl.bandwidth)
        error('Must be same number of bandwidth and CF params');
    end
    
    mdl.center_freq_khz = max(min(mdl.center_freq_khz, 0.001*mdl.sampfs/2), 0.01);
    mdl.bandwidth = min(abs(mdl.bandwidth), 20); % No more than 10 octaves
    % Convert the "octaves" now into a frequency
    bw = 0.107 * 1000 * mdl.center_freq_khz;
    for sf = fieldnames(x.dat)', sf = sf{1};
        [T, S, ~] = size(x.dat.(sf).(mdl.input));
        C = length(mdl.center_freq_khz);
        ret = zeros(T, S, C);
        for s = 1:S
            for c = 1:C                
                [BM, ENV] = gammatone(x.dat.(sf).(mdl.input)(:,s,c), ...
                                mdl.sampfs, 1000*mdl.center_freq_khz(c), ...
                                bw(c), mdl.align);
                if (mdl.use_env)
                    ret(:, s, c) = ENV;
                else
                    ret(:, s, c) = abs(hilbert(BM)); % Indistinguishable from ENV?
                end                
            end
        end
        x.dat.(sf).(mdl.output) = ret; 
    end
end

function do_plot_filtered_stim(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));  
    if mdls{1}.use_env
        ylab = 'Basilar Envelope';
    else
        ylab = 'Basilar Motion';
    end
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, mdls{1}.output, ...
            sel, 'Time [s]', ylab);
end

end
