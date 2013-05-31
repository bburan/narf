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
                     'time', 'train_BIC', 'test_BIC', 'train_ll', 'test_ll', 'cumstim'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.stim = 'stim';
m.n_bins = 200;
m.raw_ISIs = 'resp_ISIs';
m.scaled_ISIs = 'scaled_resp_ISIs';
m.time = 'stim_time';
m.train_BIC  = 'score_test_bic';
m.test_BIC  = 'score_test_bic';
m.train_ll = 'score_train_ll';
m.test_ll = 'score_train_ll';
m.cumstim = 'cumstim'; % Cumulative stimulus

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
% m.plot_fns{1}.fn = @do_plot_raw_ISIs;
% m.plot_fns{1}.pretty_name = 'Raw ISI Distribution';
% m.plot_fns{2}.fn = @do_plot_scaled_ISIs;
% m.plot_fns{2}.pretty_name = 'Scaled ISI Distribution';
% m.plot_fns{3}.fn = @do_plot_raw_corr;
% m.plot_fns{3}.pretty_name = 'Raw ISI AutoCorrelation';
% m.plot_fns{4}.fn = @do_plot_scaled_corr;
% m.plot_fns{4}.pretty_name = 'Time-Scaled ISI AutoCorrelation';
% m.plot_fns{5}.fn = @do_plot_kolmorgorov;
% m.plot_fns{5}.pretty_name = 'Kolmogorov-Smirnov Plot';

function x = do_bayesian_likelihood(mdl, x, stack, xxx)
    
    scaled_ISIs = [];
    
    fns = x.training_set;
    for ii = 1:length(fns)
        sf = x.training_set{ii};
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
            fprintf('WARNING: BIC detected negative numbers!\n');
        end
        
        [~, ~, ri] = size(x.dat.(sf).(mdl.raw_ISIs));
        
        % Find the total absolute stim integral (on a per-stimulus basis)
        for s = 1:si
        	tmp = cumsum(abs(stim(:,s)));
            CDF = tmp ./ tmp(end);
            x.dat.(sf).(mdl.cumstim)(:,s) = CDF;
            
            % Concatenate on the scaled inter-spike intervals
            scaled_ISIs = cat(1, scaled_ISIs, ...
                                 diff(interp1(time, CDF, x.dat.(sf).(mdl.raw_ISIs)(:,s))));
        end        
        
    end
    
    % Average lambda for the scaled ISIs
	l_avg = 1 / mean(scaled_ISIs);
    
	% Where scaled ISIs fall on the unit poisson's cumulative density function
    z = 1 - exp(- l_avg * scaled_ISIs);
    
    PD = fitdist(z','exponential');  % TODO: Also try gamma, inverse gaussian
    k = 1;
    n = length(isi_h); 
       
    x.(mdl.scaled_ISIs) = isi_h;
    x.(mdl.BIC) = -2*(-PD.NLogL) + 2*k*log(n);
    x.(mdl.loglikelihood) = PD.NLogL; % The Negative Log Likelihood
    
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