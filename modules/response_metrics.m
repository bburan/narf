function m = response_metrics(args)
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
m.mdl = @response_metrics;
m.name = 'response_metrics';
m.fn = @do_response_metrics;
m.pretty_name = 'Response Metrics';
m.editable_fields = {'resp', 'respavg', 'time', 'tau'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.resp    = 'resp';
m.respavg = 'respavg';
m.time    = 'resp_time';
m.tau     = 0.9;            % Time constant

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_similarity;
m.plot_fns{1}.pretty_name = 'Similarity Raster';

function x = do_response_metrics(mdl, x, stack, xxx)   
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        resp = x.dat.(sf).(mdl.resp);
        n_reps = size(resp,3);
        pairings = combnk(1:n_reps, 2);
        sample_time_delta = diff(x.dat.(sf).(mdl.time)(1:2));
        expfilt = 
        distances = nan(n_reps, n_reps);
        
        for jj = 1:length(pairings)
            fresp = filter(kernel, 1, resp);            
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