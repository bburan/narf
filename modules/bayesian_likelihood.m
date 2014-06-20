function m = bayesian_likelihood(args)
% For displaying the bayesian likelihood of a noise distribution
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @bayesian_likelihood;
m.name = 'bayesian_likelihood';
m.fn = @do_bayesian_likelihood;
m.pretty_name = 'Bayesian Likelihood';
m.editable_fields = {'stim', 'n_bins', 'raw_ISIs', 'scaled_ISIs', ...
                     'time', 'train_bic', 'test_bic', 'train_nlogl', 'test_nlogl', 'probdist', 'probcutoff'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.stim = 'stim';
m.n_bins = 500;
m.raw_ISIs = 'resp_ISIs';
m.raw_spiketimes = 'resp_spiketimes';
m.scaled_ISIs = 'scaled_resp_ISIs';
m.time = 'stim_time';
m.train_aic = 'score_train_aic';
m.test_aic = 'score_test_aic';
m.train_bic  = 'score_train_bic';
m.test_bic  = 'score_test_bic';
m.train_nlogl = 'score_train_nlogl';
m.test_nlogl = 'score_test_nlogl';
m.train_autocorr = 'score_train_autocorr';
m.test_autocorr = 'score_test_autocorr';
m.probdist = 'exponential'; % Also try 'inversegaussian', 'gamma'
m.probcutoff = 0.001; % Cutoff in seconds

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.stim, m.time, m.raw_ISIs, m.raw_spiketimes};   % Signal dependencies
m.modifies = {m.train_nlogl, m.train_bic, m.train_aic, m.train_autocorr, m.test_nlogl, m.test_bic m.test_aic m.test_autocorr};          % These signals are modified

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_scaled_isis;
m.plot_fns{1}.fn = @do_plot_scaled_isis;
m.plot_fns{1}.pretty_name = 'Scaled ISI Distribution';
m.plot_fns{2}.fn = @do_plot_scaled_autocorr;
m.plot_fns{2}.pretty_name = 'Scaled ISI AutoCorrelation';
% m.plot_fns{5}.fn = @do_plot_kolmorgorov;
% m.plot_fns{5}.pretty_name = 'Kolmogorov-Smirnov Plot';

function [nlogl, scaled_ISIs] = helper_fn(mdl, x, stimfiles)       
    global STACK;
    
    scaled_ISIs = [];
    total_time = 0;
    
    if isempty(stimfiles)
        nlogl = 10^10;
        return
    end
    
    for ii = 1:length(stimfiles)
        sf = stimfiles{ii};
        stim = x.dat.(sf).(mdl.stim);
        time = x.dat.(sf).(mdl.time);
        
        
        [ti, si, ci] = size(stim);
        if ci > 1
            error('LL expected only a single channel STIM');
        end
        
        if any(isnan(stim(:)))
            error('LL cannot handle NaNs in STIM yet.');
        end
        
        if any(stim(:) < 0)
            %fprintf('WARNING: LL detected negative numbers! Rectifying....\n');
            stim(stim<0) = abs(stim(stim<0)); % Is this the right thing to do?
        end
        
        if all(stim(:) == 0)
            fprintf('WARNING: Skipping NLOGL calculation because stim is all zero...\n');  
            nlogl = 10^10;
            return
        end
        
        [ssi, ri] = size(x.dat.(sf).(mdl.raw_spiketimes));
        if ssi ~= si
            keyboard;
            error('Incorrectly sized raw_spiketimes cell array.');
        end
        
        % Find the total absolute stim integral (on a per-stimulus basis)
        for s = 1:si
            % Total predicted spiking activity of model, unnormalized
        	CDF = cumsum(stim(:,s));
           
            % Build up the scaled inter-spike intervals distribution
           for r = 1:ri
                spiketimes = x.dat.(sf).(mdl.raw_spiketimes){s,r};                                                
                sISIs = diff(interp1q([0; time], [0; CDF], [0; spiketimes]));

                % FIXME: Sometimes scaled ISIs of 0 occur
                % I'm not sure what to do about this.
                % Right now I'll just 'modify' them to be slightly nonzero
                sISIs(sISIs < 0) = 10^-9; % FIXME
    
                if any(isnan(sISIs))
                    keyboard;
                    error('How did a scaled ISI become NaN?!');
                end
               
                scaled_ISIs = cat(1, scaled_ISIs, sISIs);
           end
        end    
        
         % We compare to RESPAVG, so we don't need to multiply by ri
        total_time = total_time + si * (time(end));
        
    end
    
    % Anything less than the temporal cutoff should be NaN'd
    scaled_ISIs(scaled_ISIs <= mdl.probcutoff) = NaN;
    scaled_ISIs = excise(scaled_ISIs);
    
    % Everything should also be shifted left to the cutoff
    scaled_ISIs(:) = scaled_ISIs(:) - mdl.probcutoff;
    
    n_spikes = length(scaled_ISIs);
    
    % Average lambda for the stimulus-scaled ISIs
	l_avg = 1 / nanmean(scaled_ISIs); 
    
    % Average lambda we should expect, from looking at RESP
    lambda = total_time / n_spikes;
    
	% Map scaled ISIs (x-axis) to the uniform distribution (Y-axis)
    %z1 = 1 - exp(- l_avg * scaled_ISIs); % Old way normalizes to stim
    z2 = 1 - exp(- lambda * scaled_ISIs); % Normalizing to respavg
    z2 = sort(z2);
    
    % Expected: uniform distribution 
    uni = linspace(0, 1, n_spikes)';     
        
    % Pick L0, L1, L2 norm here       
    %     global HID;
    %     if isempty(HID)
    %         HID = figure();    
    %     else
    %         figure(HID);        
    %         plot([0,1], [0,1], 'k--', z2, uni, 'b-');
    %     end
    
    err = sum((z2 - uni).^2); % Sum of squares error works
    % err = sum(abs(z2 - uni)) / n_spikes; % L1 error worked pretty well
    %err = max(abs(z2 - uni)); % L0 error may also work (KS plot-style)
    nlogl = err;
    
    % Negative log likelihood would be the sum of all the ISIs
    %nlogl = - sum(log(l_tot .* exp(-l_tot .* scaled_ISIs)));
    
    % Mathematically it's impossible to have z become less than zero, but
    % occasionally numerical noise is bringing us there... I think.     
    %z(z < 0) = 10^-9; 
        
    % PD = fitdist(scaled_ISIs, mdl.probdist);    
    % PD = fitdist(z2, mdl.probdist);    
    % k = length(pack_fittables(STACK));
    % n = numel(scaled_ISIs);
    
    %nlogl = PD.NLogL;
    %bic = -2*(-PD.NLogL) + k*log(n);
    %aic = -2*(-PD.NLogL) + 2*k;
    
    %tmp = excise(scaled_ISIs);
    %R = corrcoef(tmp(2:end), tmp(1:end-1));
    %autocorr = R(2,1);    
end


function x = do_bayesian_likelihood(mdl, x)
    [train_nlogl, training_scaled_ISIs] = helper_fn(mdl, x, x.training_set);
    [test_nlogl, test_scaled_ISIs] = helper_fn(mdl, x, x.test_set);    
    
    x.(mdl.train_nlogl) = train_nlogl;    
    x.(mdl.test_nlogl) = test_nlogl;
           
    x.training_scaled_ISIs = training_scaled_ISIs;
    x.test_scaled_ISIs = test_scaled_ISIs;

end

function do_plot_scaled_isis(sel, stack, xxx)
    x = xxx{end};
    mdl = stack{end}{1};
    
    % Choose which ISI set to use
    if any(strcmp(sel.stimfile, x.training_set))
        isis = x.training_scaled_ISIs;  
    else    
        isis = x.test_scaled_ISIs;
    end
    
    histfit(isis, mdl.n_bins, mdl.probdist)        
    do_xlabel('Scaled Inter-Spike Intervals [s]');
    do_ylabel('# of neurons');
    
    textLoc(sprintf('train nlogl:%f\ntest nlogl:%f', ...
                x.(mdl.train_nlogl), x.(mdl.test_nlogl)), 'NorthEast');
    
end

function do_plot_scaled_autocorr(sel, stack, xxx)
    x = xxx{end};
    mdl = stack{end}{1};
    
    % Choose which ISI set to use
    if any(strcmp(sel.stimfile, x.training_set))
        isis = x.training_scaled_ISIs;  
    else    
        isis = x.test_scaled_ISIs;
    end
        
    %train_autocorr = x.(mdl.train_autocorr);
    %test_autocorr = x.(mdl.test_autocorr);
    plot(isis(2:end), isis(1:end-1), 'k.');
     
    %textLoc(sprintf('Spike Count:%d\nTrainAutoCorr=%f\nTestAutoCorr=%f,', ...
    %                length(isis), train_autocorr, test_autocorr), 'NorthEast');
    
    do_xlabel('ISI at t(n) [s]');
    do_ylabel('ISI at t(n-1) [s]');
end

% function do_plot_raw_ISIs(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     [sf, ~, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);  
%     
%     hist(dat.(mdl.raw_ISIs)(:), nbins);
%     %xlabel('Inter Spike Interval [s]');
%     %ylabel('N');
% end
% 
% function do_plot_scaled_ISIs(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     [sf, ~, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);  
%     
%     hold on;
%     hist(dat.(mdl.scaled_ISIs)(:), nbins);
%     PDTMP = fitdist(isis_scaled', 'exponential');
%     lam = 1 / PDTMP.mu;
%     plot(ts, nbins * lam * exp(-lam*ts), 'r-');
%     title('Time-Scaled ISI and Scaled Unit Exponential');
%     %xlabel('Scaled Inter-Spike Interval [-]');
%     %ylabel('N');
%     hold off;    
% end
% 
% function do_plot_raw_corr(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     [sf, ~, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);  
%     
%     z = dat.(mdl.raw_ISIs)(:);
%     plot(z(2:end), z(1:end-1), 'k.');
% %    R = corrcoef(z(2:end), z(1:end-1));
% %     title(sprintf('ISI Correlation: %f', R(1,2)));
% %     xlabel('z_i');
% %     ylabel('z_{i-1}');
% %     hold off;
% end
% 
% function do_plot_scaled_corr(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     [sf, ~, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);  
%     
%     z = dat.(mdl.scaled_ISIs)(:);  
%     plot(z(2:end), z(1:end-1), 'k.');
% %    R = corrcoef(z(2:end), z(1:end-1));
% %     title(sprintf('ISI Correlation: %f', R(1,2)));
% %     xlabel('z_i');
% %     ylabel('z_{i-1}');
% %    hold off;
% end
% 
% function do_plot_kolmorgorov(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);  
%     
%     z = dat.(mdl.scaled_ISIs)(:);  
%     hold on;
%     plot([0,1],[0,1], 'k-', [0,0.95], [0.05,1], 'r-', [0.05,1], [0, 0.95], 'r-');
%     plot([0:1/(length(z)-1):1], sort(z), 'b-');
% %     title('Kolmogorov-Smirnov Plot');
% %     xlabel('Empirical CDF');
% %     ylabel('Spiking Stochasticity CDF');    
%     text(0, 0.9, sprintf(' BIC: %f', x.(mdl.score)), ...
%         'VerticalAlignment','top',...
%         'HorizontalAlignment','left');
%     axis tight;
%     hold off;
% end

end
