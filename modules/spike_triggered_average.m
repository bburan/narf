function m = spike_triggered_average(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @spike_triggered_average;
m.name = 'spike_triggered_average';
m.fn = @do_spike_triggered_average;
m.pretty_name = 'Spike Triggered Average';
m.editable_fields = {'stim', 'resp', 'pretrigger_time', 'n_stim_bins', 'stim_time', 'resp_time'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.stim  = 'stim';
m.resp  = 'resp'; 
m.pretrigger_time = 0.2;     % Amount of history to record in seconds
m.n_stim_bins = 30;          % How many stimulus level bins to use
m.stim_time = 'stim_time';
m.resp_time = 'resp_time';
m.st_matrix = 'st_matrix'; % Spike triggered activity matrix

% Optional fields
% m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_triggered_heatmap;
m.plot_fns{1}.pretty_name = 'Triggered Heatmap';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.stim, m.resp, m.stim_time, m.resp_time};   % Signal dependencies
m.modifies = {m.st_matrix};         % These signals are modified

% Initialize a HISTORY
% For each repetition
%   for every instant in time
%       If time bin contains a spike 
%           Add the previous 500 bin values to each of the histories
%           Do this N times if the bin contains N spikes
function x = do_spike_triggered_average(mdl, x)
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        stim = x.dat.(sf).(mdl.stim);
        resp = x.dat.(sf).(mdl.resp);
        stim_time =  x.dat.(sf).(mdl.stim_time);  
        resp_time =  x.dat.(sf).(mdl.resp_time);  
        
        [ti1, si1, ci1] = size(stim);
        [ti2, si2, ri2] = size(resp);
        if (ti1 ~= ti2 || si1 ~= si2)
            error('Stim must be the same size as resp!');
            % Throw an error if stim and RESP aren't the same resolution & size? 
            % TODO: do the "smart thing" and interpret as best we can
        end
        
        % Calculate the vertical bin sizes
        top = max(stim(:));
        bot = min(stim(:));
        vbins = linspace(bot, top, mdl.n_stim_bins);
       
        % Raster the stimulus into bins
        [~, stim_binned] = histc(stim, vbins);               
        
        % Using stim_time, estimate the sampling frequency
        stim_fs = 1 / (stim_time(2) - stim_time(1));       
        n_time_bins = ceil(stim_fs * mdl.pretrigger_time);
        
        x.dat.(sf).(mdl.st_matrix) = zeros(ci1, mdl.n_stim_bins, n_time_bins);

        ctr = 1;
        ctr_end = ci1 * si2 * ri2;
        
        for c = 1:ci1 % Make a different STA for each channel
            for s = 1:si2
                for r = 1:ri2
                    trial = resp(:,s,r);
                    indexes_of_spikes = find(trial > 0);
                    
                    fprintf('%d/%d\n', ctr, ctr_end); ctr = ctr+1;
                    for zz = 1:length(indexes_of_spikes)
                        
                        % This loop accounts for when two spikes occur
                        % in the same raster period. 
                        idx = indexes_of_spikes(zz);
                        for yy = 1:resp(idx, s,r)

                            spike_time = resp_time(idx);
                            
                            % Build a vector of the pre-spike stimulus
                            idxs = (spike_time - mdl.pretrigger_time < stim_time & stim_time < spike_time);
                            prespike_stim_binned = flipud(stim_binned(idxs));
                            % prespike_stim_times = stim_time(idxs);
                            
                            if (length(prespike_stim_binned) > n_time_bins)
                                keyboard;                                
                            end
                            
                            for pp = 1:length(prespike_stim_binned),
                                xcoord = pp;
                                ycoord = prespike_stim_binned(pp);
                                %fprintf('[%d, %d]\n', xcoord, ycoord);
                                x.dat.(sf).(mdl.st_matrix)(c, ycoord, xcoord) = 1 + x.dat.(sf).(mdl.st_matrix)(c, ycoord, xcoord);    
                            end
                        end
                    end
                end
            end
            
            % Once the STA has been built, go back and scale each column
            % by the number of events added to each column
            hmap = squeeze(x.dat.(sf).(mdl.st_matrix)(c,:,:));
            for aa = 1:size(hmap,2)
                x.dat.(sf).(mdl.st_matrix)(c, :, :) = hmap(:, aa) ./ sum(hmap(:,aa));
            end
        end
    end
end

function do_plot_triggered_heatmap(sel, stack, xxx)
    x = xxx{end};
    mdl = stack{end}{1};

    hmap = squeeze(x.dat.(sel.stimfile).(mdl.st_matrix)(sel.chan_idx, :,:)); 
      
    stim = x.dat.(sel.stimfile).(mdl.stim);
	resp = x.dat.(sel.stimfile).(mdl.resp);
    top = max(stim(:));
    bot = min(stim(:));
    vbins = linspace(bot, top, mdl.n_stim_bins);
        
    % Plot a MEAN CURVE that is twice as hot as any other color. 
    % vec = 1:size(hmap,2);
    % HOT=[];
    % for ii = vec,
    %     % TODO: Plot a very hot MEAN curve (which should basically be the 
    %     hmap(ii) = hmap(:, ii) .* vec;
    % end
    
    % ERASE LAST LINE and 0 stimulus line for better visibility
    %hmap(:, end) = 0;
    %hmap(1,:) = 0;
    
    % figure;
    imagesc(hmap);
    % hold on;    
    % plot(1:size(hmap,2), mdl.n_stim_bins * mdl.n_stim_bins * mean(hmap), 'Color','r','LineWidth',2);  
    % hold off;
    
    do_xlabel('Time Before Spike [s]');
    do_ylabel('Stimulus Intensity [-]');
end

end
