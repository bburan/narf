% Time-scaling method for generating spikes

clear all;
close all;

duration = 3;                % Duration in seconds
Fs = 1000;                   % Sampling frequency
ts = [1:Fs*duration]*(1/Fs); % Time series
secret = 0.9;                % This is the number we will try to infer

% lambda becomes a time-varying function
lambda = @(t) 123*(1 + cos(2*pi*secret*t));

% Integral of lambda
Lambda = @(t) (123*t + (123/(2*pi))*sin(2*pi*secret*t));

% The last time
T_end = Lambda(ts(end));

% How many times to sample?  TODO: This is probably wrong!
samples_to_take = duration * T_end;

% Generate spikes using random sampling
unsorted_spike_times = zeros(1,samples_to_take);

for i = 1:samples_to_take
    % Pick a random point on the uniform distribution between 0 and T_scale
    z = T_end * rand(1,1);
              
    % Invert Lambda to find the point in time which corresponds to z
    t_spike = fzero(@(t)(z - Lambda(t)), 0.5); 
    
    % Save that value
    unsorted_spike_times(i) = t_spike;   
end

spikes = sort(unsorted_spike_times);

% Inter Spike Intervals
isis_raw = diff([0, spikes]);

% Time scaled inter spike intervals: % TODO THIS IS NOT BETWEEN 0 and 1
isis_scaled = diff(Lambda([0,spikes]));

% Average lambda for the scaled ISIs
lambda_avg = 1 / mean(isis_scaled);

% Where scaled ISIs fall on the unit poisson's cumulative density function
z = 1 - exp(- lambda_avg * isis_scaled);

% How many bins for the histogram?
nbins = 50;

% Plot things prettily
figure;

subplot(3,2,1);
hold on;
bar(spikes, 50*ones(1,length(spikes)), 0.01,'k');
plot(ts, lambda(ts), 'r-');
title('Non-stationary Poisson: Histogram (black), PDF (red)');
xlabel('Time [s]');
ylabel('\lambda(t) [spikes/s]');
hold off;

subplot(3,2,2);
hist(isis_raw, nbins);
title('Raw ISI');
xlabel('Inter Spike Interval [s]');
ylabel('N');

subplot(3,2,3);
plot(ts, Lambda(ts));
xlabel('Time [s]');
ylabel('Cumulative Density');
title('{\Lambda}(t)');

subplot(3,2,4);
hold on;
hist(isis_scaled, nbins);
PDTMP = fitdist(isis_scaled', 'exponential');
lam = 1 / PDTMP.mu;
plot(ts, nbins * lam * exp(-lam*ts), 'r-');
title('Time-Scaled ISI and Scaled Unit Exponential');
xlabel('Scaled Inter-Spike Interval [-]');
ylabel('N');
hold off;

subplot(3,2,5);
hold on
plot([0,1],[0,1], 'k-', [0,0.95], [0.05,1], 'r-', [0.05,1], [0, 0.95], 'r-');
plot([0:1/(length(z)-1):1], sort(z), 'b-');
title('Kolmogorov-Smirnov Plot');
xlabel('Empirical CDF from Rescaled Time');
ylabel('Spiking Stochasticity CDF');
hold off;

subplot(3,2,6);
plot(z(2:end), z(1:end-1), 'k.');
R = corrcoef(z(2:end), z(1:end-1));
title(sprintf('Time-Scaled ISI Correlation: %f', R(1,2)));
xlabel('z_i');
ylabel('z_{i-1}');
hold off;

% -----------------------------------------

% Now check the likelihoods to see how inference would work

i=1;
search = 0.1:0.01:1.4;
likelihoods = [];

% Exhaustive search example
for l = search
    
    % Define the CDF of our hypothesis
    hypothesis_CDF = @(t) (123*t + (123/(2*pi*l))*sin(2*pi*l*t));
  
    % Time-Scale the spikes according to our CDF
    isi_h = diff(hypothesis_CDF([0,spikes]));

    % Average lambda for the scaled ISIs
    l_avg = 1 / mean(isi_h);

    % Where scaled ISIs fall on the unit poisson's cumulative density function
    z = 1 - exp(- l_avg * isi_h);
    z = sort(z);  % This is the all-important scaled CDF
    
    % METHOD #1: Least Squares deviation from 45 degree line on KS Plot
    %likelihoods(i) = - sum(([0:1/(length(z)-1):1] - z).^2); 
    
    % METHOD #2: Negative Log Likelihood
    PD = fitdist(z','exponential');
    k = 1;
    n = length(isi_h); 
    %likelihoods(i) = -2*(-PD.NLogL) + 2*k*log(n); % THIS IS THE BIC!
    likelihoods(i) = PD.NLogL; % The Negative Log Likelihood

       
    i=i+1;
end

figure; plot(search, likelihoods, 'k-');
title('Model Likelihood (Higher is better)');
xlabel('Model Parameter');
ylabel('Negative Log Likelihood');

% figure;
% hold on;
% bar(spikes, 50*ones(1,length(spikes)), 0.01,'k');
% plot(ts, lambda(ts), 'r-');
% l=0.99; hypo = @(t) 123*(1 + cos(2*pi*l*t));
% plot(ts, hypo(ts), 'g-');
% l=1.02; hypo = @(t) 123*(1 + cos(2*pi*l*t));
% plot(ts, hypo(ts), 'g-');
% xlabel('Time [s]');
% ylabel('\lambda(t) [spikes/s]');
% hold off;
% 
% 
% 
