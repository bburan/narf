function m = bernoulli_trials(args)
% A module that maps stimulus intensity to RESPAVG as Bernoulli trials.
% For best accuracy, use at a high enough sampling rate that there is 
% generally only one spike per bin. 
% 
% Basically, a Bayesian version of the Spike Triggered Average? Suffers
% from the same problem of autocorrelation of stimulus #($*&ing everything
% up. 
%
% Idea: Try to remove autocorrelation somehow? Or differentiate stimulus
% before doing STA? I'm not sure what to do here. 
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @bernoulli_trials;
m.name = 'bernoulli_trials';
m.fn = @do_bernoulli_trials;
m.pretty_name = 'Bernoulli Trial';
m.editable_fields = {'n_stimbins', 'n_respbins', 'delay_bins', ...
                     'stim', 'respavg'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.n_stimbins = 20;  % How many bins to bin the stimulus into
m.n_respbins = 200; % Vertical "probabilistic" bins for bernoulli
m.delay_bins = 0; % From the first bin
m.stim = 'stim';
m.respavg = 'respavg';
m.bernoulli_grid = [];
m.is_perf_metric = true;

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_bernoulli_grid;
m.plot_fns{1}.fn = @do_plot_bernoulli_grid;
m.plot_fns{1}.pretty_name = 'Bernoulli Grid';

function [posterior, z, N] = bernoulli_bayes(theta, prior, data)
    z = sum(data == 1);
    N = length(data);
    % likelihood = theta.^z .* (1-theta).^(N-z);   % p(D|theta)
    % Use log likelihood because likelihood gets too small
    loglikelihood = z.*log(theta) + (N-z).*log(1-theta);
    likelihood = exp(loglikelihood); % Still has some problems
    evidence = sum(likelihood .* prior);   % Marginalizes out    
    posterior = likelihood .* prior ./ evidence; % Bayes Theorem
    % BROKEN ATTEMPT TO USE LOGS all the way through.
    %logprior = log(prior);
    %posterior = logprior + loglikelihood - sum(loglikelihood + logprior);
    %prob = z/N;
end

function x = do_bernoulli_trials(mdl, x, stack, xxx)    
    % TODO: Use test data as well
    stim = flatten_field(x.dat, x.training_set, mdl.stim);
    resp = flatten_field(x.dat, x.training_set, mdl.respavg);      
    
    coefs = zeros(500,3);
    
    for ii = 1:size(coefs,1)
        fprintf('Find quadratic approximation for delay = %d \n', ii);
        mdl.delay_bins = ii - 1;
        %ii = mdl.delay_bins + 1;
        
        % Delay the stimulus by some amount, truncate the response
        s = stim(1+mdl.delay_bins:end);
        r = resp(1:end-mdl.delay_bins);
        D = [s r]; 
    
        % Set up probabilistic analysis
        stim_bin_edges  = linspace(min(s), max(s), mdl.n_stimbins + 1);
        grid = linspace(0, 1, mdl.n_respbins); % Bernoulli variable theta
    
        bernoulli_grid = nan(mdl.n_stimbins, mdl.n_respbins);
        % raw_prob = [];
    
        % Compute posterior density of theta for a given intensity level    
        for b = 1:mdl.n_stimbins 
            mask = stim_bin_edges(b) <= s & s <= stim_bin_edges(b+1);
            relevant_data = D(mask, 2);
            relevant_data(relevant_data > 0) = 1; 
            relevant_data(relevant_data <= 0) = 0;  % Make it 1 or 0
    
            prior = ones(size(grid))/length(grid);
            [posterior, z, N] = bernoulli_bayes(grid, prior, relevant_data);
            bernoulli_grid(b, :) = posterior;
        end
    
        % SMOOTH the distribution (ignoring relative frequencies of stim data)
        sigma = 0.02*mdl.n_stimbins; % Sigma is the "neighbor smoothing" value. 
        wid = ceil(sigma*3);  % Number of bins wide for discrete filter
        ww = linspace(-wid/2, wid/2, wid);
        gaussFilter = exp(-ww .^ 2 / (2 * sigma ^ 2));
        gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    
        yfilt = conv2(bernoulli_grid, gaussFilter', 'same'); % MATLAB is Broken
        % Following workaround not needed: CONV2 does the "right thing"?
        % yfilt = nan(size(bernoulli_grid));
        % for ii = 1:mdl.n_respbins
        %     yfilt(:, ii) = iconv(bernoulli_grid(:, ii), gaussFilter);
        % end        
    
        x.bernoulli_grid = yfilt;
        %x.bernoulli_grid = bernoulli_grid;
    
        % TODO: This "smoothed" 2D probability is no longer normalized
        % properly. Vertical columns may no longer sum up to 1.  
    
        % Pull out the MAP model from the bernoulli grid, and scale to get prob
        [~, MAP_indexes] = max(yfilt, [], 2);
        MAP_model = MAP_indexes / mdl.n_respbins;
    
        %figure; 
        
        xx = linspace(min(s), max(s), mdl.n_stimbins)';
        %plot(xx, MAP_model, 'k.');
        %hold on;
        % TODO: Technically this is the WRONG thing to do since it weights all
        % MAP points equally, which is not true. MAP points should be weighted
        % by their posterior probabilities during "linear" fits, etc. 
    
        % Generate a linear fit
        %P = polyfit(xx, MAP_model, 1); 
        %plot(xx, P(1)*xx + P(2), 'r-');    
        %coefs(ii, :) = [0 P];
        
        % Generate a quadratic fit
        P = polyfit(xx, MAP_model, 2);
        coefs(ii, :) = P;
        %plot(xx, P(1)*xx.^2 + P(2)*xx + P(3), 'b-');
    
        % CUBIC
        % P = polyfit(xx, MAP_model, 3);
        % plot(xx, P(1)*xx.^3 + P(2)*xx.^2 + P(3)*xx + P(4), 'g-');
    
        % Fit a SIGMOID to the damn thing
        %g = fittype('base + (peak - base) / (1 + exp( gain * (- x + thresh)))');
        %p0 = [];
        %p0.baserate = 0;
        %p0.peakrate = 0.5;
        %p0.gain = 1;
        %p0.thresh = 0.0;
        %myFit = fit(xx, MAP_model, g, 'Startpoint', [0.05 2 0.4 2]);
        %plot(xx, myFit.base + (myFit.peak - myFit.base) / (1 + exp(myFit.gain * (- xx + myFit.thresh))), 'g-');
        % plot(xx, 0.05 + 1 / (1 - exp( - 0.1 * (xx + 2))), 'g-');
    
        %hold off;
    end
    
     yy = [];
     if true
         t = 0:0.01:3;
         figure;
         hold on;
         len = size(coefs,1);
         for ii = 1:len
             yy(ii,:) = coefs(ii, 1)*t.^2 + coefs(ii, 2)*t + coefs(ii, 3);
             plot(t, yy(ii),...
                  'Color', [ii/len 0 (len-ii)/len]);
         end
         hold off;
     end
     figure;
     mesh(yy);        
     
    %     Extract the marginalized probability of a spike at once bin latency
    %     Multiply, sum, or average with (not sure yet) other bins' predictions
    %     That should give you the "spike probability distribution" at each
    %       moment in time
    
%     fprintf('Computing predictions based on firing probabilities...\n');
%     for sf = fieldnames(x.dat)', sf=sf{1};
%         in = x.dat.(sf).(mdl.stim);
%         [T, S, C] = size(in);
%         if all(isnan(in(:)))
%             yout = in(:);
%         else
%             for tt = 1:T
%                 
%                 yout = npfnl(in(:));
%             end
%         end
%         
%         x.dat.(sf).(mdl.output) = reshape(yout,[T,S,C]);
%     end
    
    % OPTIONAL: Plot a heat map of the probability of a given response
    
    % OPTIONAL: Compute the overall likelihood of the data given that
    %           prediction. 
    
end


function do_plot_bernoulli_grid(sel, stack, xxx)    
    x = xxx{end};
    mdl = stack{end}{1};    
    
    %x.bernoulli_grid
    %figure;
    imagesc(x.bernoulli_grid(:, 1:end)');
    ca = caxis;
    lim = max(abs(ca));
    caxis([0, +lim]);
    
    set(gca,'YDir','normal');
    do_xlabel('STIM Intensity [-]');
    do_ylabel('Bernoulli Variable (\theta)[-]');    
    %title('Color is posterior belief of theta.');
    
    
end

end
