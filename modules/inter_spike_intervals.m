function m = inter_spike_intervals(args)
% Compute the inter-spike intervals (ISIs) of a neural response.
% Or maybe more precisely, the ISIs returned are actually just the
%  "time since the previous spike"

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @inter_spike_intervals;
m.name = 'inter_spike_intervals';
m.fn = @do_inter_spike_intervals;
m.pretty_name = 'Inter Spike Intervals';
m.editable_fields = {'input', 'n_bins', 'time', 'output', 'output_time'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input  = 'resp';
m.n_bins = 500;
m.time   = 'resp_time';
m.output = 'resp_ISIs'; % ISIs are measured from one spike to its previous
m.output_time = 'resp_spiketimes';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_isis;
m.plot_fns{1}.pretty_name = 'Inter-spike Intervals';
m.plot_fns{2}.fn = @do_plot_autocorrelation;
m.plot_fns{2}.pretty_name = 'Auto-Correlation';
m.plot_fns{3}.fn = @do_plot_raw_spikes;
m.plot_fns{3}.pretty_name = 'Spikes vs. Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output, m.output_time};   % These signals are modified

function x = do_inter_spike_intervals(mdl, x)     
   fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        in = x.dat.(sf).(mdl.input);
        time =  x.dat.(sf).(mdl.time);                
        [ti, si, ri] = size(in);
        x.dat.(sf).(mdl.output) = zeros(ti, si, ri);
        x.dat.(sf).(mdl.output_time) = zeros(ti, si, ri);
        for s = 1:si
            for r = 1:ri
                % If there are no NaNs, we can do this the fast way        
                if ~any(isnan(in(:,s,r)))
                    sidxs = in(:,s,r) > 0;
                    spiketimes = time(sidxs);
                    x.dat.(sf).(mdl.output_time)(sidxs,s,r) = spiketimes;
                    prevs = [spiketimes(:); 0] - [0; spiketimes(:)];
                    x.dat.(sf).(mdl.output)(sidxs,s,r) = prevs(1:end-1);
                else
                    % Otherwise, I think we must use for loops and reset
                    % every time we have NaNs
                    tprev = time(1);
                    for t = 1:ti
                        if isnan(in(t,s,r))
                            tprev = time(t);
                        end
                        if in(t,s,r) > 0
                            x.dat.(sf).(mdl.output_time)(t,s,r) = time(t);
                            x.dat.(sf).(mdl.output)(t,s,r) = time(t) - tprev;
                            tprev = time(t);
                        end
                    end 
                end                              
            end
        end
    end
end

function do_plot_raw_spikes(sel, stack, xxx)
    mdl = stack{end}{1};
    sel.chan_idx = []; 
    do_plot(xxx(end), mdl.time, mdl.input, ...
            sel, 'Time [s]', 'Spike Events [-]');
end

function do_plot_isis(sel, stack, xxx)
    % [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    xout = xxx{end};
    mdl = stack{end}{1};
     
    sidxs = (xout.dat.(sel.stimfile).(mdl.input) > 0);
    hist(xout.dat.(sel.stimfile).(mdl.output)(sidxs), mdl.n_bins);
 
    do_xlabel('Raw Inter-Spike Intervals [s]');
    do_ylabel('# of neurons');
    
end

function do_plot_autocorrelation(sel, stack, xxx)
    % [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    xout = xxx{end};
    mdl = stack{end}{1};
    
    sidxs = (xout.dat.(sel.stimfile).(mdl.input) > 0);
    isis = xout.dat.(sel.stimfile).(mdl.output)(sidxs);

    plot(isis(2:end), isis(1:end-1), 'k.');
 
    R = corrcoef(isis(2:end), isis(1:end-1));
    textLoc(sprintf('Spike Count:%d\nSelected Stimfile Corr=%f', ...
                    nnz(sidxs), R(2,1)), 'NorthEast');
    do_xlabel('ISI at t(n) [s]');
    do_ylabel('ISI at t(n-1) [s]');
end

end
