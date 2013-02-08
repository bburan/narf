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
                     'time',  'BIC', 'loglikelihood'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.stim = 'stim';
m.n_bins = 50;
m.raw_ISIs = 'resp_ISIs';
m.scaled_ISIs = 'scaled_resp_ISIs';
m.time = 'stim_time';
m.BIC  = 'score_bic';
m.loglikelihood = 'score_ll';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_raw_ISIs;
m.plot_fns{1}.pretty_name = 'Raw ISI Distribution';
m.plot_fns{2}.fn = @do_plot_scaled_ISIs;
m.plot_fns{2}.pretty_name = 'Time-Scaled ISI Distribution';
m.plot_fns{3}.fn = @do_plot_raw_corr;
m.plot_fns{3}.pretty_name = 'Raw ISI AutoCorrelation';
m.plot_fns{4}.fn = @do_plot_scaled_corr;
m.plot_fns{4}.pretty_name = 'Time-Scaled ISI AutoCorrelation';
m.plot_fns{5}.fn = @do_plot_kolmorgorov;
m.plot_fns{5}.pretty_name = 'Kolmogorov-Smirnov Plot';

function x = do_bayesian_likelihood(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Hypothesis CDF is the numerical integral of the stim signal
    hypothesis_CDF = @(t) (123*t + (123/(2*pi*l))*sin(2*pi*l*t));
  
    % Compute the log likelihood
    for ii = 1:length(x.training_set) 
        sf = x.training_set{ii};
        [T S C] = size(x.dat.(sf).(mdl.input1));
        for s = 1:S, 
            % DO stuff
        end
    end
        
    % Time-Scale the spikes according to our CDF
    isi_h = diff(hypothesis_CDF([0,spikes]));

    % Average lambda for the scaled ISIs
    l_avg = 1 / mean(isi_h);

    % Where scaled ISIs fall on the unit poisson's cumulative density function
    z = 1 - exp(- l_avg * isi_h);
    
    % The all-important scaled CDF
    z = sort(z);
    
    PD = fitdist(z','exponential');
    k = 1;
    n = length(isi_h); 
       
    x.(mdl.scaled_ISIs) = isi_h;
    x.(mdl.BIC) = -2*(-PD.NLogL) + 2*k*log(n);
    x.(mdl.loglikelihood) = PD.NLogL; % The Negative Log Likelihood
    
end

function do_plot_raw_ISIs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    hist(dat.(mdl.raw_ISIs)(:), nbins);
    %xlabel('Inter Spike Interval [s]');
    %ylabel('N');
end

function do_plot_scaled_ISIs(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    hold on;
    hist(dat.(mdl.scaled_ISIs)(:), nbins);
    PDTMP = fitdist(isis_scaled', 'exponential');
    lam = 1 / PDTMP.mu;
    plot(ts, nbins * lam * exp(-lam*ts), 'r-');
    title('Time-Scaled ISI and Scaled Unit Exponential');
    %xlabel('Scaled Inter-Spike Interval [-]');
    %ylabel('N');
    hold off;    
end

function do_plot_raw_corr(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    z = dat.(mdl.raw_ISIs)(:);
    plot(z(2:end), z(1:end-1), 'k.');
%    R = corrcoef(z(2:end), z(1:end-1));
%     title(sprintf('ISI Correlation: %f', R(1,2)));
%     xlabel('z_i');
%     ylabel('z_{i-1}');
%     hold off;
end

function do_plot_scaled_corr(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    z = dat.(mdl.scaled_ISIs)(:);  
    plot(z(2:end), z(1:end-1), 'k.');
%    R = corrcoef(z(2:end), z(1:end-1));
%     title(sprintf('ISI Correlation: %f', R(1,2)));
%     xlabel('z_i');
%     ylabel('z_{i-1}');
%    hold off;
end

function do_plot_kolmorgorov(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
    dat = x.dat.(sf);  
    
    z = dat.(mdl.scaled_ISIs)(:);  
    hold on;
    plot([0,1],[0,1], 'k-', [0,0.95], [0.05,1], 'r-', [0.05,1], [0, 0.95], 'r-');
    plot([0:1/(length(z)-1):1], sort(z), 'b-');
%     title('Kolmogorov-Smirnov Plot');
%     xlabel('Empirical CDF');
%     ylabel('Spiking Stochasticity CDF');    
    text(0, 0.9, sprintf(' BIC: %f', x.(mdl.score)), ...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');
    axis tight;
    hold off;
end

end