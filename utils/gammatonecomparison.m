% Comparison of methods to create gamma tone spectrographs

%import gammatonegram.ERBSpace;
%import gammatonegram.ERBFilterBank;
%import gammatonegram.MakeERBFilters;
%import gammatonebank.gammatonebank;

% Load a waveform, calculate its gammatone spectrogram, then display:
[X,SR] = wavread('sa2.wav');

TWIN=0.025;  % Summation window size
THOP=0.010;  % Space between window centers (10ms)
N=128;       % Gamma tone channels
FMIN=100;    % Lowest frequency
FMAX=SR/2;   % Highest (nyquist) frequency
WIDTH=1.0;

% First, we try it with one filter bank
%
[fcoefs,F] = MakeERBFilters(SR, N, FMIN); % Make the filter bank 
fcoefs = flipud(fcoefs);                  % Flip coef matrix up/down
F = flipud(F);
XF = ERBFilterBank(X,fcoefs);             % Process the waveform

nwin = round(TWIN*SR);                    % The window size, in samples
XE = [XF.^2];                             % Square to get energy
hopsamps = round(THOP*SR);                %
ncols = 1 + floor((size(XE,2)-nwin)/hopsamps);
Y = zeros(N,ncols);

for i = 1:ncols
    Y(:,i) = sqrt(mean(XE(:,(i-1)*hopsamps + [1:nwin]),2));
end

% Then we try it with a second filter bank system

align=true; % Align the gamma filters?
ZF=gammatonebank(X,FMIN,FMAX,N,SR,align);
ZE = [ZF.^2];                             % Square to get energy
Z = zeros(N,ncols);
for i = 1:ncols
    Z(:,i) = sqrt(mean(ZE(:,(i-1)*hopsamps + [1:nwin]),2));
end

% Now plot everything!
subplot(311)
logfsgram(X, 1024, SR, nwin, round(nwin/2), FMIN, 12);
%spectrogram(X,nwin,nwin/2,SR,'yaxis');   % Linear scale for comparison
%caxis([-90 -30])                         % Manually limit color range
colorbar
set(gca,'XTickLabel', 1000*(get(gca,'XTick')));
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Square-window FFT');

subplot(312)
imagesc(20*log10(Y)); axis xy
%caxis([-90 -30])                        
colorbar
set(gca,'YTickLabel',round(F(get(gca,'YTick'))));
set(gca,'XTickLabel', 1000*THOP*(get(gca,'XTick')));
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Gammatonegram');

subplot(313)
imagesc(20*log10(Z)); axis xy
%caxis([-90 -30])
colorbar
set(gca,'YTickLabel',round(F(get(gca,'YTick'))));
set(gca,'XTickLabel', 1000*THOP*(get(gca,'XTick')));
xlabel('Time [ms]');
ylabel('Frequency [Hz]');
title('Another Gammatonegram (Phase locked filter bank)');

