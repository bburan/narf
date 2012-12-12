% This is just a quick test to figure out how to find an optimal impulse
% response from matlab!

[cfd, cellids, cellfileids] = dbgetscellfile('cellid','dai008a-c1','rawid', 64662);
index = 1;
channel = 1;
stimfs = 20000;      % Stimulus frequency
respfs = 200;        % Response histogram bin frequency (200Hz -> 5ms) 

options = []; 
options.includeprestim = 0;  % Include pre and post stimulation results
options.unit     = cfd(index).unit;
options.channel  = channel;
options.rasterfs = respfs; 

% Load and raster the auditory stimulus
stimfile = [cfd(index).stimpath cfd(index).stimfile];
stimenv  = loadstimfrombaphy(stimfile, [], [], 'envelope', options.rasterfs, 1, 0, options.includeprestim);

% Load and raster the neural spike response
respfile = [cfd(index).path cfd(index).respfile]; 
[resp, tags] = loadspikeraster(respfile, options); 

X = stimenv(1,:,1)';
Y = sum(resp(:,:,1),2);

mktimeaxis = @(x, fs) ( [0:(1/fs):(length(x) - 1)/fs]);
T = mktimeaxis(X, respfs);
	
% Well, that didn't work very well! I'll try again with the Kalman filtered
% data some time....    
[psth_avg, psth_top, psth_bot] = spikerate.KalmanPSTH(Y, respfs);

%g = filter([1 1 1 1], [1 0 0 0], psth_avg);
%g = smooth(psth_avg, 10);

%figure; plot(T,X, 'b-', T, Y, 'r-', T, psth_avg', 'g-');
figure; plot(T,X, 'b-', T, psth_avg'/6, 'g-');

rawdata = iddata(Y, X, 1/respfs);
kaldata = iddata(psth_avg', X, 1/respfs);
%figure; impulse(m,T);

% ir = cra(thedata);        
% m = arx(thedata,[2 2 1]);
% imp = [1;zeros(19,1)];
% irth = sim(m,imp);
% figure;
% subplot(211)
% plot(irth) 
% title('impulse responses')

% Technique #1: impulse(). Gives a big ugly, non-causal mess
impulse(kaldata);  

% Technique #2: correlation regression. 
%    Over the last 50 elements.
%    10th order model
%    1=plot
% CRA() filters the input to remove correlations and 'whiten' everything
% That should theoretically improve the least squares fit?
% Result is a another big mess, only slightly better
[imp_cra, R_cra, cl] = cra(kaldata, 50, 10, 1);

% Technique #3: Autoregressive model with external input (ARX)
%  y is input, 
%  u is output, 
%  filter order parameters specified with [na, nb, nk]
%   na     y's history of itself
%   nb     y's dependence on u's history
%   nk     dead time between k and 
na = 0;
nb = 20;
nk = delayest(kaldata);   % Estimates the dead time 
km = arx(kaldata, [na nb nk]);
kimp = [1;zeros(nb-1,1)];
imp_kal = sim(km, kimp);


% Now the same thing, but for the raw data
nk = delayest(rawdata);   % Estimates the dead time 
rm = arx(rawdata, [na nb nk]);
rimp = [1;zeros(nb-1,1)];
imp_raw = sim(rm, rimp);

ts=[1:nb];
plot(ts, imp_kal, 'g-', ts, imp_raw, 'r-');
title('Estimated impulse response for KalmanPSTH (green) and raw data (red), 5 params');


% Now get the predicted neural responses for each
rawpred = sim(rm, X);
kalpred = sim(km, X);
mypred  = filter([0 1 0 0], [1 0 0 0], psth_avg); % Just a simple delay

size(kalpred)
size(mypred)

% Evaluation time:
R = corrcoef([X,Y,psth_avg',rawpred, kalpred, mypred']);
disp('Correlations were:');
disp(R(1,:));
% X-Y      0.3711 
% X-kal    0.6029 
% X-linear 0.8642 
% X-linkal 0.9039   Wheee!
% X-shift  0.5722

% But those are probably waaay overfitting. A better metric would use
% cross-validation of the model. TODO.

% -------------------------------------------------------
% Plot correlations as a big matrix of graphs

corrplot([X,Y,psth_avg',rawpred, kalpred], {'Envelope', 'PSTH', 'Kal', 'PSTHLin', 'KalLin'});



% -------------------------------------------------------
% Plot the correlation matrix as a heat map?
imagesc(corrsort(R)); % plot the matrix
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', L); % set x-axis labels
set(gca, 'YTickLabel', L); % set y-axis labels
title('Your Title Here', 'FontSize', 14); % set title
colormap('jet'); % set the colorscheme
colorbar; % enable colorbar


% then to set the axis titles you'll have to use
% Please note the curly braces for the cell array
labelNames = {'X','PSTH','KalPSTH', 'PSTH-pred', 'Kal-Pred', 'mypred'};
set(gca,'XTickLabel',labelNames);   % gca gets the current axis
set(gca,'YTickLabel'labelNames);   % gca gets the current axis

% -------------------------------------------------------


% I think I need a jackknifed, leave-one-out correlation score too:
igma = 5;
y = normrnd(0,sigma,100,1);
m = jackknife(@var, y, 1);
n = length(y);
bias = -sigma^2 / n % known bias formula
jbias = (n - 1)*(mean(m)-var(y,1)) % jackknife bias estimate

% 
% % Now try to filter the input with that with kirth!
% Z = filter(kirth, 5, X);
% figure; plot(T,X, 'b-', T, Z, 'g-', T, psth_avg'/6, 'r-');
% xlabel('Time')
% ylabel('Intensity');
% title('Envelope (blue), linear filtered envelope (green), and kalman PSTH (red)')
% 
% R = corrcoef(Y,Z)
