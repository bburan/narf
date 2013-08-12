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
                     'time', 'train_bic', 'test_bic', 'train_nlogl', 'test_nlogl', 'cstim', 'probdist', 'probcutoff'};
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
m.cstim = 'cstim'; % Cumulative stimulus

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_scaled_isis;
m.plot_fns{1}.fn = @do_plot_scaled_isis;
m.plot_fns{1}.pretty_name = 'Scaled ISI Distribution';
m.plot_fns{2}.fn = @do_plot_scaled_autocorr;
m.plot_fns{2}.pretty_name = 'Scaled ISI AutoCorrelation';
% m.plot_fns{5}.fn = @do_plot_kolmorgorov;
% m.plot_fns{5}.pretty_name = 'Kolmogorov-Smirnov Plot';

function [x, nlogl, bic, aic, autocorr] = helper_fn(mdl, x, stimfiles)       
    global STACK;
    scaled_ISIs = [];
    for ii = 1:length(stimfiles)
        sf = stimfiles{ii};
        stim = x.dat.(sf).(mdl.stim);
        time = x.dat.(sf).(mdl.time);
        
        [ti, si, ci] = size(stim);
        if ci > 1
            error('BIC/LL expected only a single channel STIM');
        end
        
        if any(isnan(stim(:)))
            error('BIC/LL cannot handle NaNs in STIM yet.');
        end
        
        if any(stim(:) < 0)
            % fprintf('WARNING: BIC/NLOGL detected negative numbers! Rectifying....\n');
            stim(stim<0) = 0;
        end
        
        if all(stim(:) == 0)
            fprintf('WARNING: Skipping BIC/NLOGL calculation because stim is all zero...\n');  
            x.dat.(sf).(mdl.scaled_ISIs) = [];
            nlogl = inf;
            bic = 0; 
            return
        end
        
        [~, ~, ri] = size(x.dat.(sf).(mdl.raw_ISIs));
        
        % Find the total absolute stim integral (on a per-stimulus basis)
        for s = 1:si
        	tmp = cumsum(abs(stim(:,s)));
            CDF = tmp ./ tmp(end);
            x.dat.(sf).(mdl.cstim)(:,s) = CDF;
           
            % Build up the scaled inter-spike intervals distribution
           for r = 1:ri
                v = x.dat.(sf).(mdl.raw_ISIs)(:,s,r);
                idxs = v > 0;
                spiketimes = x.dat.(sf).(mdl.raw_spiketimes)(idxs,s,r);
                sISIs = diff(interp1([0; time], [0; CDF], [0;spiketimes]));

                % FIXME: Sometimes scaled ISIs of 0 occur
                % I'm not sure what to do about this.
                % Right now I'll just 'modify' them to be slightly nonzero
                sISIs(sISIs < 0) = 10^-9; % FIXME
    
                if any(isnan(sISIs))
                    error('How did a scaled ISI become NaN?!');
                end
                
                scaled_ISIs = cat(1, scaled_ISIs, sISIs);
           end
        end
    end
    
    % Anything less than the cutoff should be NaN'd
    scaled_ISIs(scaled_ISIs <= mdl.probcutoff) = NaN;   
    
    % Everything should also be shifted left to the cutoff
    scaled_ISIs(:) = scaled_ISIs(:) - mdl.probcutoff;
    
    % Average lambda for the scaled ISIs
	l_avg = 1 / nanmean(scaled_ISIs);
    
	% Where scaled ISIs fall on the unit poisson's cumulative density function
    z = 1 - exp(- l_avg * scaled_ISIs);      
    
    % Mathematically it's impossible to have z become less than zero, but
    % occasionally numerical noise is bringing us there... I think.     
    z(z < 0) = 10^-9; 
    
    PD = fitdist(z, mdl.probdist);    
    k = length(pack_fittables(STACK));
    n = numel(scaled_ISIs);
    
    x.dat.(sf).(mdl.scaled_ISIs) = scaled_ISIs;        
    nlogl = PD.NLogL;    
    bic = -2*(-PD.NLogL) + k*log(n);
    aic = -2*(-PD.NLogL) + 2*k;
    
    tmp = excise(scaled_ISIs);
    R = corrcoef(tmp(2:end), tmp(1:end-1));
    autocorr = R(2,1);    
    
end


function x = do_bayesian_likelihood(mdl, x, stack, xxx)
    
    [x, train_nlogl, train_bic, train_aic, train_autocorr] = helper_fn(mdl, x, x.training_set);
    [x, test_nlogl, test_bic, test_aic, test_autocorr] = helper_fn(mdl, x, x.test_set);    
    
    x.(mdl.train_nlogl) = train_nlogl;
    x.(mdl.train_bic) = train_bic;       
    x.(mdl.train_aic) = train_aic;
    x.(mdl.train_autocorr) = train_autocorr;
    
    x.(mdl.test_nlogl) = test_nlogl;
    x.(mdl.test_bic) = test_bic;   
    x.(mdl.test_aic) = test_aic;
    x.(mdl.test_autocorr) = test_autocorr;
       
end

function do_plot_scaled_isis(sel, stack, xxx)
    % [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    xout = xxx{end};
    mdl = stack{end}{1};
     
    %sidxs = (xout.dat.(sel.stimfile).(mdl.scaled_ISIs) > 0);
    tmp = xout.dat.(sel.stimfile).(mdl.scaled_ISIs);
    tmp(tmp < mdl.probcutoff) = NaN;
    tmp(:) = tmp(:) - mdl.probcutoff;
    %hold on;
    %hist(tmp, mdl.n_bins);
    histfit(tmp, mdl.n_bins, mdl.probdist)
    %PDTMP = fitdist(xout.(mdl.scaled_ISIs), mdl.probdist);
    %lam = 1 / PDTMP.mu;
    %ts = linspace(0, max(xout.(mdl.scaled_ISIs)(:)), 500);
    %plot(ts, length(tmp)/mdl.n_bins * lam * exp(-lam*ts), 'r-');
    %hold off;

    do_xlabel('Scaled Inter-Spike Intervals [s]');
    do_ylabel('# of neurons');
    
    textLoc(sprintf('train nlogl:%f\ntest nlogl:%f', ...
                xout.(mdl.train_nlogl), xout.(mdl.test_nlogl)), 'NorthEast');
    
end

function do_plot_scaled_autocorr(sel, stack, xxx)
    % [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    xout = xxx{end};
    mdl = stack{end}{1};
    
    isis = xout.dat.(sel.stimfile).(mdl.scaled_ISIs);
    train_autocorr = xout.(mdl.train_autocorr);
    test_autocorr = xout.(mdl.test_autocorr);
    plot(isis(2:end), isis(1:end-1), 'k.');
     
    textLoc(sprintf('Spike Count:%d\nTrainAutoCorr=%f\nTestAutoCorr=%f,', ...
                    length(isis), train_autocorr, test_autocorr), 'NorthEast');
    
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