function plot_polezeros_vs_fir(batch, cellids, modelnames)
% THIS IS A VERY HACKY FUNCTION
% I AM JUST USING IT TEMPORARILY
% YOU MUST SELECT ONLY A SINGLE CELLID

global XXX STACK META;

if length(cellids) ~= 1
    error('Only one cellid can be selected.');
end

n_subplots = length(modelnames);
idx_subplot = 1;
fig = figure('Name', cellids{1}, ...
             'NumberTitle', 'off', 'Position', [20 50 900 1000]);
         
    function w = pack_polezero_impulses(unused_stack)        
        firfn = [];
        ok_load = true;
        ok_ds = true;
        ok_wc = true;
        sr = 200;
        imp = {};        
        weights = [];
        fprintf('%s\n', XXX{1}.cellid);
        
        for ii = 1:length(STACK)
            mm = STACK{ii};
            nsplits = length(mm);
            if nsplits > 1
                error('I cannot support split/unify things yet.');
            end
            m = mm{1};
            
            if strcmp(m.name, 'load_stim_resps_from_baphy')
                if ok_load
                    fprintf('LSRFB: Changing sr to %f\n', m.raw_stim_fs);
                    sr = m.raw_stim_fs;
                    ok_load = false;
                else
                    error('too many Loadstimsfrombaphy');
                end
            end
            
            if strcmp(m.name, 'downsample_signal')
                if ok_ds
                    fprintf('DS: Changing sr to %f\n', m.output_freq);
                    sr = m.output_freq;
                    ok_ds = false;
                else
                    error('too many downsamples');
                end
            end
            
            t = 0:0.0001:0.1;
            
            if strcmp(m.name, 'pole_zeros')
                p = {};
                z = {};
                for ii = 1:m.n_inputs
                    p{ii} = m.poles(ii, :);
                    z{ii} = m.zeros(ii, :);
                end
                
                sys = zpk(z, p, m.gains);
                sys.InputDelay = abs(m.delays) / 1000; % (milliseconds)
                imp{end+1} = impulse(sys, t);
                
                % Find out the total area under the curve and normalize
                tt = 0:0.00001:10; % A ridiculous number of samples
                mm = nanmean(impulse(sys, tt));
                imp{end} = imp{end} / (10*mm);
                
            end
            
            if strcmp(m.name, 'weight_channels')
                if ok_wc
                    weights = m.weights;
                    ok_wc = false;
                else
                    error('Too many weight_channels');
                end
            end
            
            if strcmp(m.name, 'fir_filter')
                subplot(n_subplots, 1, idx_subplot);
                idx_subplot = idx_subplot + 1;
                m.plot_fns{1}.fn({}, STACK(1:ii), {}); 
                w = [];
                hl = xlabel(META.modelname); set(hl,'interpreter','none');
                return;
            end
        end
               
        if length(weights) ~= length(imp);
            if length(imp) == 1 && isempty(weights)
                weights(1) = 1; % Single filter case
            else
                error('Number of weights and PZ filters don''t match!?');
            end
        end        
            
        total = zeros(size(imp{1}));
        for ii = 1:length(weights)
            total = total + weights(ii) * imp{ii};
        end
        subplot(n_subplots, 1, idx_subplot);
        idx_subplot = idx_subplot + 1;
        imagesc(total');
        set(gca,'YDir','normal');
        ca = caxis;
        lim = max(abs(ca));
        caxis([-lim, +lim]);
        hl = xlabel(META.modelname); set(hl,'interpreter','none');
        w = [];
        
    end

stack_extractor = @pack_polezero_impulses;

[params, ~, ~] = load_model_batch(batch, cellids, modelnames, ...       
                                  stack_extractor);

%hold on;
%plot(allzx, allzy, 'r.');
%plot(allpx, allpy, 'b.');
%set(gca,'XScale','log');
%set(gca,'YScale','log');
%axis tight;
%xlabel('Latency');
%ylabel('Pole or Zero');
%hold off;
%
%figure('Name', 'Poles/Zeros Hist', 'NumberTitle', 'off', ...
%       'Position', [20 50 900 900]);
%data = log10(allpy(:));
%data(:,2) = NaN;
%data(1:length(allzy),2) = log10(allzy);
%[y,b] = hist(data, 100);
%bar(b,y, 'grouped', 'BarWidth', 1);
%title('Pole (blue) and Zero (red) Histogram');
%
%xlabel('Log10(pole) or log10(zero)');
%ylabel('# of poles/zeros');

end