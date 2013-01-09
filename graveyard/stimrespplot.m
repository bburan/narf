function [] = stimrespplot(stim, stim_fs, resp, resp_fs, stimenv)
% STIMRESPPLOT: Generates a stimulus/response plot for a neuron.
% 
% The stimulus/response plot has five subplots on it:
%   1. Raw signal vs time
%   2. Gammatonegram                     (TODO: vary nonlinear wavelets)
%   3. PSTH + Kalman Filtered Response   (TODO: vary smoothing)
%   4. Model-based MaxLikelihood input
%      Gammatone filtered sound intensity + PSTH correlation
%      Raw signal envelope + PSTH correlation
%   5. Inverted STRF                     (correlation + inverted wavelet)
% 
% FUNCTION ARGUMENTS:
%   stim     Stimulus wave form    [1xW vector]
%   stim_fs  Stimulus sample rate  [Hz]
%   resp     Response spike counts [MxN, M: sample time, N: spike count] 
%   resp_fs  Response sample rate  [Hz]
%
% REUTURNS:
%   fh       Figure handle
%   
% Notes
%   1. By chopping up output correlations into smaller pieces, 
%      we can work out what a linear STRF must have been.  
%   2. Modeling stimulus/response intensity logarithmically might help. 
%
% 2012-09-20, Ivar Thorson.

% TODO: Check that input arguments are of equal length
% TODO: Check that the proper number of arguments were provided

% -----------------------------------------------------------------------
% User-settable values. TODO: Make arguments later?
alignphase = false;   % Align the gamma filters to be in phase?
t_win = 0.025;        % Gamma filter window size [seconds]
t_hop = 1/resp_fs;    % Time between gamma window centers [seconds]
N_gfs = 64;           % Number of gamma tone filter channels
f_min = 100;          % Lowest frequency [Hz]
f_max = stim_fs/2;    % Highest (nyquist) frequency [Hz]
N_xticks = 10;        % Number of ticks to display on the x axis
N_yticks = 5;         % Number of ticks to display on the y axis

% -----------------------------------------------------------------------
% Values derived mathematically from the above. 
mktimeaxis = @(x, fs) ( [0:(1/fs):(length(x) - 1)/fs]);
t_stim = mktimeaxis(stim,stim_fs);     % Stimulus time axis
t_resp = mktimeaxis(resp,resp_fs);     % Response time axis
N_win = round(t_win * stim_fs);        % Number of samples in each window
N_hop = round(t_hop * stim_fs);        % Number of samples to hop center by
N_cols = 1 + floor((length(stim) - N_win) / N_hop); % Total number of windows
gamma_pow = zeros(N_gfs, N_cols);

% Filter the stimulus sound using a gammatone bank
[gamma_resp, gamma_envs, gamma_frqs] = ...
    gammatonebank(stim, f_min, f_max, N_gfs, stim_fs, alignphase);

% Rectify and smooth the gamma response for visualization
gamma_sqr = [gamma_resp.^2];                      
for i = 1:N_cols
    gamma_pow(:,i) = 20*log10(sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2)));
    %gamma_pow(:,i) = sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2));
end

% Run a Kalman filter on the peristimulus raw spiking data
[psth_avg, psth_top, psth_bot] = ...
    KalmanPSTH(resp);

% Let's plot the figures!
fh = figure;

% -----------------------------------------------------------------------
%subplot(511);
% The Stimulus! Just show the original sound
% Options: Play sound
subplot('Position', [0.05 0.8 0.9 0.18]);
plot(t_stim, stim, 'b-');
set(gca,'XLim',[0 length(stim)/stim_fs]);
setAxisLabelCallback(gca, @(x)(1000*x), 'X');
xlabel('Time [ms]');
ylabel('Sound Pressure');
h = title('Raw Stimulus Waveform');
P = get(h,'Position'); set(h,'Position',[P(1) P(2)-10 P(3)]);

% -----------------------------------------------------------------------
% subplot(512);
% Linear Wavelet filters:  (Shape of the cochlea, transduction)
%   Square Window (Pure FFT)
%   Gammatone
%   SinC (low pass)
%   Lanczos filter.
%   Mexican hat, 
%   Meyer (First orthogonal wavelet?)
subplot('Position', [0.05 0.6 0.9 0.18]);
imagesc(gamma_pow); axis xy; 
% caxis(clipOutliers(gamma_pow,1));
% colorbar;
setAxisLabelCallback(gca, @(x)(1000*t_hop*x), 'X');
setAxisLabelCallback(gca, @(y)round(gamma_frqs(y)), 'Y');
% set(gca,'XLim',[0 length(gamma_pow)]);
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
if alignphase
    h=title('Gammatone Filter Response (Phase locked filters)');
else
    h=title('Gammatone Filter Response');
end
P = get(h,'Position'); set(h,'Position',[P(1) P(2)-10 P(3)]);

% -----------------------------------------------------------------------
% Nonlinear Model (Kalman Filtered)
%   Pre-synaptic depression (input state)
%   Logarithmic intensity compression (input filter)
%   Square root intensity compression (input filter)
%   Refractory periods (model state)
%   Observational model (Poisson)

% -----------------------------------------------------------------------
%subplot(513);
% The response!
% Options:
%   Binning
%   Convolution with Gaussian, Exponential or boxcar kernels
%   Model-based PSTH
subplot('Position', [0.05 0.4 0.9 0.18]);
hold on;
bar(t_resp, resp, 0.01,'k');
plot(t_resp, psth_avg,'b-');
plot(t_resp, psth_top,'r-');
plot(t_resp, psth_bot,'r-');
set(gca,'XLim',[0 length(resp)/resp_fs]);
setAxisLabelCallback(gca, @(x)(1000*x), 'X');
xlabel('Time [ms]');
ylabel('Spike Rate [-]');
h=title('PSTH + Kalman Filtered Poisson Model');
P = get(h,'Position'); set(h,'Position',[P(1) P(2)-10 P(3)]);
hold off;

% -----------------------------------------------------------------------
%subplot(514);
% Correlations with filtered values
% Options:
%    Whole-signal correlations
%    Wavelet correlations
%      Entire list of previous correlations
%      Symmetric or causal?
subplot('Position', [0.05 0.2 0.9 0.18]);
% The old way: 
% plot(t_resp, stimenv, 'k-');
% The newer, double axis way
hl1 = line(t_resp, stimenv, 'Color','r');
ax1 = gca;
set(ax1,'XColor','r','YColor','r')

% Adjust the axes as before
set(gca,'XLim',[0 length(resp)/resp_fs]);
setAxisLabelCallback(gca, @(x)(1000*x), 'X');
xlabel('Time [ms]');
ylabel('Stimulus Envelope');

ax2 = axes('Position', get(ax1, 'Position'), ... 
        'XAxisLocation', 'top', ...
        'YAxisLocation', 'right', ...
        'Color', 'none', ...           
        'XColor', 'k', ...
        'YColor', 'k');
     

thediff = length(resp) - length(gamma_pow);
halfdiff = floor(thediff/2);
thefront= ones(1, halfdiff) * gamma_pow(55,1);
theback = ones(1, thediff- halfdiff) * gamma_pow(55,1);
temppow = [thefront, gamma_pow(55,:), theback];

disp(size(resp));      %% 300x1
disp(size(psth_avg));  %% 1x300
disp(size(temppow)); %% 1x298
disp(size(stimenv));         %% 1x300

hl2 = line(t_resp, temppow, 'Color','k','Parent',ax2);
ylabel('Power of Best-Fitting Gammatone Filter');
% -----------------------------------------------------------------------
% subplot(515);
% Error
%    Predicted minus actual response. 
%    
subplot('Position', [0.05 0.0 0.9 0.18]);
% Inverted STRF
%plot((1/stim_fs)*1:length(gamma_envs(:,55)), gamma_pow(:,35), 'r-', ...
%    (1/stim_fs)*1:length(gamma_envs(:,55)), gamma_pow(:,25), 'b-', ...
%    (1/stim_fs)*1:length(gamma_envs(:,55)), gamma_pow(:,30), 'g-');
%ll = length(gamma_pow(55,:));
%plot((3/ll)*1:ll, gamma_pow(54:56,:));
%ch = get(gcf,'ch');
%set(ch(2:5),'xtick',[])
    

R=corrcoef(resp', stimenv);
disp(sprintf('Correlation psth-env:          %f', R(2,1)));

[r,d]=bestcorrcoef(resp', stimenv);
disp(sprintf('Correlation psth-env:          %f (delay: %d)', r, d));

R=corrcoef(psth_avg, stimenv);
disp(sprintf('Correlation kalman-env:        %f', R(2,1)));

[r,d]=bestcorrcoef(psth_avg, stimenv);
disp(sprintf('Correlation kalman-env:        %f (delay: %d)', r, d));

R=corrcoef(log(psth_avg), stimenv);
disp(sprintf('Correlation log(kalman)-env:   %f', R(2,1)));

[r,d]=bestcorrcoef(log(psth_avg), stimenv);
disp(sprintf('Correlation log(kalman)-env:   %f (delay: %d)', r, d));


R=corrcoef(resp', temppow);
disp(sprintf('Correlation psth-gamma:        %f', R(2,1)));

[r,d]=bestcorrcoef(resp', temppow);
disp(sprintf('Correlation pst-gamma:         %f (delay: %d)', r, d));

R=corrcoef(psth_avg, temppow);
disp(sprintf('Correlation kalman-gamma       %f', R(2,1)));

[r,d]=bestcorrcoef(psth_avg, temppow);
disp(sprintf('Correlation kalman-gamma:      %f (delay: %d)', r, d));

R=corrcoef(log(psth_avg), temppow);
disp(sprintf('Correlation log(kalman)-gamma: %f', R(2,1)));

[r,d]=bestcorrcoef(log(psth_avg), temppow);
disp(sprintf('Correlation log(kalman)-gamma: %f (delay: %d)', r, d));

R = corrcoef(stimenv, temppow);
disp(sprintf('Correlation env-gamma:         %f', R(2,1)))  ;


R = corrcoef(log(stimenv), temppow);
disp(sprintf('Correlation logenv-gamma:      %f', R(2,1)))  ;

%figure; plot(resp', stimenv, 'ko');
%figure; plot(psth_avg, stimenv, 'ro');
%figure; plot(resp', temppow, 'bo');
%figure; plot(psth_avg, temppow, 'go');
%   figure; plot(sqrt(psth_avg), temppow, 'go');
%figure; plot(log(psth_avg), temppow, 'go');

figure; corrplot([stimenv', temppow', resp, psth_avg']);
