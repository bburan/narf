function m = neural_statistics(args)
% A module which computes neural statistics for the RESP signal:
%   Total number of spikes
%   Spike distance between each trial (Determinism)
%
% Doesn't compare the neuron to anything in particular
% 
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @neural_statistics;
m.name = 'neural_statistics';
m.fn = @do_neural_statistics;
m.pretty_name = 'Neural Statistics';
m.editable_fields = {'resp', 'respavg', 'time', 'spike_count'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.resp    = 'resp';
m.respavg = 'respavg';
m.time    = 'resp_time';
m.spike_count = 'spike_count';
m.tau_range   = logspace(-3,-0.5,100); % From 0.001 s to 0.3 sec time constants
m.similarity  = 'metric_self_similarity';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_similarity;
m.plot_fns{1}.pretty_name = 'Inputs, Error vs Time';

function x = do_neural_statistics(mdl, x, stack, xxx)   
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        resp = x.dat.(sf).(mdl.resp);
        
        spike_count = sum(resp(:));
        x.dat.(sf).(mdl.spike_count) = spike_count;
    
        % convolve RESP with exponential kernels of size tau
        pairings = combnk(1:size(resp,3), 2);        
        avgdist = nan(size(mdl.tau_range));
        for jj = 1:length(mdl.tau_range)
            tau = mdl.tau_range(jj);
            sample_time_delta = diff(x.dat.(sf).(mdl.time)(1:2));
            % nsamps = ceil(-log(10^-6)*tau/sample_time_delta); % The efficient, proper way to do this
            nsamps = 10000; % Trading extra computation to ensure max precision
            zz = 0:min(nsamps-1, length(x.dat.(sf).(mdl.time)) - 1);
            kernel = exp((-zz.*sample_time_delta)./tau);    
            fresp = filter(kernel, 1, resp);
            distances = nan(size(pairings,1), 1);
            
            for kk = 1:size(pairings,1)
                v1 = fresp(:, :, pairings(kk,1));
                v2 = fresp(:, :, pairings(kk,2));                
                distances(kk) = nanmean((v1(:) - v2(:)).^2);
            end
            avgdist(jj) = 1000 * nanmean(distances(:)) / (sum(kernel) * spike_count);
        end
    
        x.dat.(sf).tau_range = mdl.tau_range(:);
        x.dat.(sf).(mdl.similarity) = avgdist(:);
        [themin, ind] = min(avgdist(:));
        x.dat.(sf).best_tau = mdl.tau_range(ind);
        x.dat.(sf).min_dist = themin;
    end
end

function do_plot_similarity(sel, stack, xxx)    
    mdls = stack{end};
    xouts = xxx(end);
    hold on;
    do_plot(xouts, 'tau_range', mdls{1}.similarity, ...
            sel, 'Tau', 'Average Spike Distance');    
        
    textLoc(sprintf('Best Tau: %f\nMin Dist: %f', ...
                    xouts{1}.dat.(sel.stimfile).best_tau, ...
                    xouts{1}.dat.(sel.stimfile).min_dist), 'SouthWest');
    hold off;
end

end