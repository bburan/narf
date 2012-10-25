% Nonstationary Poisson problem: 
% Given that rate Lambda changes with time, what is the best estimate?

% A stationary poisson process 
duration = 3; % seconds
lambda =  100; % Spikes/sec
Fs = 500;    % Sampling frequency
ts = [1:Fs*duration]*(1/Fs); % Time series

% Time 'til Next Spike: (constant lambda ONLY)
ttns_stationary = @(lam) (- log(rand(1,1)) / lam);

% Generate spikes
t=0.0;
spikes = zeros(1,Fs*duration);
rawspiketimes = [];
while (t < duration)
    t = t + ttns_stationary(lambda);
    i = ceil(t*Fs);        % Samples appear at the END of a measurement time
    if (i < duration*Fs)
        spikes(i) = spikes(i) + 1;
    else
        break;
    end
end

% Infer the spike rate
lambda_est = sum(spikes)/duration;

% Plot spike events
figure;
hold on;
bar(ts, spikes, 0.01,'k');
plot(ts, (lambda/Fs)*ones(1,duration*Fs), 'r-');
plot(ts, (lambda_est/Fs)*ones(1,duration*Fs), 'g-');
title('Stationary Poisson: histogram (black), true (red), inferred (green)');
hold off;

% Plot the interspike intervals (ISIs)
figure;
hist(diff(spikes));


% -----------------------------------------------------------------------
% A nonstationary (aka nonhomogeneous) poisson process

% Lambda becomes a time-varying function
lambda = @(t) 100*(1 + cos(2*pi*t));
    
% A plot of lambda vs time
lambdas = lambda([1/Fs:1/Fs:duration]);

% To compute the time until the next spike we need to find 
% the Inverse of "Complementary Cumulative Distribution Function".
CCDF = @(t) (100*t + (100/(2*pi))*sin(2*pi*t));

% Generate spikes using forward integration method
t = 0.0;
spikes = zeros(1,Fs*duration);
    
while (t < duration)
    % Pick a random point on the uniform distribution
    z = rand(1,1);
       
    % Define the integrating estimate of the CCDF
    %F = @(t_next) exp(-quad(lambda, t, max(t_next, t))); % Slow 
    F = @(t_next) exp(-(CCDF(max(t,t_next)) - CCDF(t))); % Faster 
        
    % Invert F numerically to compute the next time
    t_next = fzero(@(t)(z - F(t)), 0.5); 
    
    t = t_next;
    i = ceil(t*Fs);

    if (i < duration*Fs)
        spikes(i) = spikes(i) + 1;
    else
        break;
    end
end

% Estimate as if it were stationary
lambda_est = sum(spikes) / duration;

% Estimate with the Kalman filter
[km, kt, kb] = spikerate.KalmanPSTH(spikes, Fs);

% Plot spike events
figure;
hold on;
bar(ts, spikes, 0.01,'k');
plot(ts, (1/Fs)*lambdas, 'r-');
plot(ts, (lambda_est/Fs)*ones(1,duration*Fs), 'g-');
plot(ts, [km; kt; kb], 'b-');
title('Non-stationary Poisson: histogram (black), true (red), stationary inferred (green), Kalman Filter (blue)');
hold off;

% ------------------------------------------------------------------------
% And last but not least, inference!

% What is the likelihood of the model which produced the spiking data?
% Since each point
%      P(lambda_ab, t_a, t_b) = lambda_ab * exp ^ -lambda_ab * (t_b - t_a)
%
% We should be able to find the log-probability of lambda(t) simply by
% taking the sum of the above.
% Compute:
%     lambda_ab(t_a, t_b) =   integral of model lambda, from t_a to t_b

% P(DATA|MODEL) = P(MODEL|DATA) P(DATA) / P(MODEL)

% ------------------------------------------------------------------------

% Unrasterize the spike timings from a histogram (hopefully very little
% binning occured, but maybe there are 2 or 3 spikes per bin, and we must
% get rid of that)

k = 1;
dt = 1/Fs;
spike_times=[]; 
for i = 1:length(spikes)
    if spikes(i) ~= 0
        % Divide up the spikes times evenly for that bin
        for j = 1:spikes(i);            
            spike_times(k) = ts(i) - dt + dt*(j / spikes(i));
            k=k+1;
        end
    end
end

% TODO TEST: Rasterizing and unrasterizing should yield exactly the same data!
        
spike_gaps = [spike_times(1), diff(spike_times)];
i=1;
mylamb = [];
likelihoods = [];
% Exhaustive search:
for l = 0.1:0.01:2
    
    % Create the hypothesis
    % hypo_lambda = @(t) 100*(1 + cos(2*pi*l*t));
  
    % Using this avoids integration to estimate lambda_ab
    hypo_CDF = @(t) (100*t + (100/(2*pi*l))*sin(2*pi*l*t));
  
    % Loop through the spikes
    logprob = 0.0;
    for j = 1:length(spike_times)
        if (i == 1)
            t_a = 0.0;
        else 
            t_a = spike_times(i-1);
        end
        
        t_b = spike_times(i); %t_a + spike_gaps(i);

        % Estimate average lambda over that period? 
        t = t_b - t_a;
        
        % lamb = (hypo_CDF(t_b) - hypo_CDF(t_a)) / ((t_b - t_a) * Fs);
        % lamb = hypo_CDF(t_b) - hypo_CDF(t_a);
        
        % Compute the log probability of this observation
        logprob = logprob + log((t*lamb^2 / 1) * exp(- lambda * t));
    end
    
    likelihoods(i) = logprob;
    
    i=i+1;
end

plot([0.1:0.01:2], likelihoods, 'k-');

% Look at lambda
l = 1.8;
hypo_CDF = @(t) (100*t + (100/(2*pi*l))*sin(2*pi*l*t));
mylamb = [hypo_CDF(spike_times(1)), ...
          hypo_CDF(spike_times(2:end)) - hypo_CDF(spike_times(1:end-1))] ./ (spike_gaps * Fs) ;

figure;
hold on;
bar(ts, spikes, 0.01,'k');
plot(ts, (1/Fs)*lambdas, 'r-');
plot(ts, (lambda_est/Fs)*ones(1,duration*Fs), 'g-');
plot(ts, [km; kt; kb], 'b-');
plot(spike_times, mylamb, 'y-');
title('Non-stationary Poisson: histogram (black), true (red), stationary inferred (green), Kalman Filter (blue)');
hold off;

% Try:
% 1. Polynomials
% 2. Parameterizing with gaussian mixture models and 1D k-means as a 
%   quick and easy way to get a smooth function
%
% What is the most likely sinusoid in light of the spiking?

% Graphing some basic curves because I don't understand:

t = [0:0.1:10]
lambda = .9;
%alpha = 1;
beta = lambda;
P1 = @(t) t*lambda^2 .* exp(-lambda * t);
P2 = @(alpha,t) (beta^alpha / gamma(alpha)) .* t.^(alpha-1) .* exp(-beta * t);
P3 = @(t) t.^2*lambda^3/2 .* exp(-lambda * t);
P4 = @(t) t.^3*lambda^4/6 .* exp(-lambda * t);

plot(t, P1(t), 'bo', ...
     t, P2(1,t), 'r-', ...
     t, P2(2,t), 'g-', ...
     t, P2(3,t), 'y-', ...
     t, P2(4,t), 'y-', ...
     t, P4(t), 'gx', ...
     t, P3(t), 'bx');


% Plot the interspike intervals


