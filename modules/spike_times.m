function m = spike_times(args)
% Compute the times of each spike

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @spike_times;
m.name = 'spike_times';
m.fn = @do_spike_times;
m.pretty_name = 'Spike Times';
m.editable_fields = {'input', 'n_bins', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input  = 'resp10000';
m.n_bins = 500;
m.time   = 'resp10000time';
m.output = 'resp_spiketimes';

% Optional fields
m.plot_fns = {};

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};   % These signals are modified

function x = do_spike_times(mdl, x)     
   fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};        
        resp = x.dat.(sf).(mdl.input);
        dt =  x.dat.(sf).(mdl.time)(2) - x.dat.(sf).(mdl.time)(1);        
        [ti, si, ri] = size(resp);
        resp = reshape(resp, [], ri); % Concat all stimuli
        time = cumsum(dt * ones(ti * si, 1));
        x.dat.(sf).(mdl.output) = cell(ri, 1);
                        
        for r = 1:ri
            % Skip a trial if it is completely NAN
            if all(isnan(resp(:,r)))
                x.dat.(sf).(mdl.output){r} = [];
                continue;
            end
            if any(isnan(resp(:,r)))
                fprintf('Warning: Mixed NANs in resp is incompatible with spike_times.m\n');
                keyboard;
            end
            x.dat.(sf).(mdl.output){r} = time(resp(:,r) > 0);
        end
    end
end

end
