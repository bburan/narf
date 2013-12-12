% Test of impulse responses of various filters
clear all;

% Three high-level parameters that (I hope) make life easy for fitters
CF   = 1000; % Peak transfer function frequency, in Hz. MUST BE > 0.
Q    = 10;   % Quality Factor. MUST BE > sqrt(2) or it will be unstable.
N    = 3;    % Order of the filter. MUST BE AN INTEGER.

% Calculate useful parameterizations
w_c = 2*pi*CF;
theta = asin(1-(1/(2*Q^2)));
w_n = w_c / sqrt(1-(1/(2*Q^2))); % Natural frequency in rad/sec

% Create a set of 2nd order poles that will be repeated N times
b   = w_n * cos(theta/2 + pi/4);
w_r = w_n * sin(theta/2 + pi/4);
poles = [-b + 1i*w_r; -b - 1i*w_r];

% Replicate the number of poles depending on the order of the filter
poles = repmat(poles, N, 1);

% Plotting limits in Hz
f_min = 100; 
f_max = 5000; 

% Create the APGF to have a gain of 1 at DC
APGF = zpk([], poles, 1);
APGF = zpk([], poles, 1/dcgain(APGF));

% Create the OZGF to have a gain that is same as APGF at CF
z = zeros(N, 1);
OZGF = zpk(z, poles, 1);
OZGF = zpk(z, poles, (1/cos(theta))/real(evalfr(OZGF, w_c*1j))); 

% Plot
h = bodeplot(APGF, OZGF);
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
legend('APGF', 'OZGF');

% TODO: If Q < 3, print a warning that the OZGF CF is probably off!

%figure; zplane(z, poles);
%figure; impulse(APGF);

%%%%%%%%%%%%%%%%%%%%%%%%

% Hmax = 20; % Peak transfer function gain at CF relative to 0Hz, in dB. MUST BE >0.
% Hmax = 10^(Hmax / 20);  % 40dB yields 100 since 20*log10(100)=40 
% Convert Hmax into non-decibel units

% Although I would love to have bandwidth be a high-level parameter, it's
% challenging to figure out how to solve the equations such that 
% the order N works out to be a reasonable integer value. 
% TODO: Revisit this problem later.
% BW   = 1;  % Bandwidth to -3dB power points, in octaves.
% theta = acos(1/(Hmax^(1/N)));

% Q = 1/(2*(1-sin(theta)));
