function [env]=ozgf_filterbank(input, lowcf, highcf, n_chans, inputfs, outputfs, align_peak)
% envs = ozgf_filterbank(input, lowcf, highcf, numchans, inputfs, outputfs,
% align_peak)
%
% input    Signal
% lowcf    Low cutoff frequency
% highcf   High cutoff frequency
% n_chans  Number of channels
% inputfs  Sample rate of input
% outputfs Sample rate of output
% align    Whether to align the peaks of the impulse responses

% Check for nyquist frequency problems
if inputfs < highcf * 2
    error('Cannot create filterbank with high cutoff frequency above nyquist maximum.');    
end
m = [];
m.zeros       = zeros(1,n_chans); % Picks a single zero that can be stacked by M_order
m.N_order     = 4;  % aka, number of pairs of poles
m.M_order     = 1;  % aka, number of stacked zeros
CFs = logspace(log10(lowcf/1000), log10(highcf/1000), n_chans);
m.center_freq_khz = CFs; % 0.2-20Khz
m.Q_factor    = ones(1, n_chans) * (1.0*CFs(1) ./ (CFs(2) - CFs(1))); % Actually the offset from 1/sqrt(2)
m.delayms     = zeros(1,n_chans); % Input delays in ms
m.time_align  = align_peak;

% Make the bank of filters
syses = {};
t_aligns = [];
for ii = 1:length(m.center_freq_khz)
    CF = abs(m.center_freq_khz(ii)) * 1000;
    Q  = 1/sqrt(2) + 0.00001 + abs(m.Q_factor(ii)); %.00001 is for numerical stability
    N  = ceil(abs(m.N_order));
    
    % Calculate useful parameterizations
    w_c = 2*pi*CF;
    theta = asin(1-(1/(2*Q^2)));
    w_n = w_c / sqrt(1-(1/(2*Q^2))); % Natural frequency in rad/sec
    
    % Create a set of 2nd order poles that will be repeated N times
    b   = w_c * cos(theta/2 + pi/4);
    w_r = w_c * sin(theta/2 + pi/4); % Using w_c instead of w_n...why?
    poles = [-b + 1i*w_r; -b - 1i*w_r];
    
    % Replicate the number of poles depending on the order of the filter
    poles = repmat(poles, N, 1);
    
    if ~isfield(m, 'zeros') || isempty(m.zeros)
        % Create the APGF to have a gain of 1 at DC
        sys = zpk([], poles, 1);
        sys = zpk([], poles, 1/dcgain(sys));
    else
        z = repmat(m.zeros(ii) * 1000, m.M_order, 1);
        sys = zpk(z, poles, 1);
        % Gain is 0 at DC for OZGF, so we just try to keep it "near" the
        % same levels by normalizing by the ringing frequency w_r. Why does
        % this get slightly erronious for very low Q factors?
        sys = zpk(z, poles, (1/cos(theta)^N) / abs(evalfr(sys, w_r*1j)));
    end
    
    % Correct for the time delay introduced by the wavelet
    if isfield(m, 'time_align') && m.time_align
        [imp, time] = impulse(sys);
        [~, idx] = max(abs(hilbert(imp)));
        t_align = time(idx);
    else
        t_align = 0;
    end
    
    sys.InputDelay = (abs(m.delayms) / 1000); % (milliseconds)
    syses{ii} = sys;
    t_aligns(ii) = t_align;
end


% Now perform the filtering

[T, S, C] = size(input);
if C ~= 1
    error('There must be only one input channel');
end

tmp = zeros(T,S,length(syses));

t = (1/inputfs).*(1:T);

for s = 1:S    
    u = squeeze(input(:, s, :));
    
    % If there are NANs in the input, treat them like 0's when
    % simulating with lsim, then NAN out the output later.
    % This is kind of an ugly hack.
    u_nan = isnan(u);
    u(u_nan) = 0;
    u(isinf(u)) = 10^6;
    warning off Control:analysis:LsimStartTime;
    warning off Control:analysis:LsimUndersampled;
    
    for ii = 1:length(syses)
        blah = lsim(syses{ii}, u, t);
        blah2 = hilbert(blah);
        tmp(:,s,ii) = abs(blah2);
        
        nanidxs = any(u_nan,2);
        tmp(nanidxs,s) = nan;
        
        % Shift everything by t_align bins
        if isfield(m, 'time_align') && m.time_align
            idx = find(t < t_aligns(ii), 1, 'last');
            tmp(:,s, ii) = [tmp(idx:end,s)' tmp(end)*ones(1, idx-1)];
        end
    end
    
    warning on Control:analysis:LsimUndersampled;
    warning on Control:analysis:LsimStartTime;
end

% Apply the post-filter function
env_inputfs = abs(tmp);

% Downsample
scale = inputfs / outputfs;

if ~(floor(scale) == scale)
    error('Input frequency must be an integer multiple of output frequency.');
end

% env = zeros(size(env_inputfs, 1)/scale, size(env_inputfs, 3));
% for ii = 1:size(env_inputfs, 3)
%     env(:,ii) = decimate(env_inputfs(:, 1, ii), scale);
% end

env_inputfs = squeeze(env_inputfs);
smfilt=ones(scale, 1)./ scale; % Low pass boxcar filter
tmp=conv2(smfilt, [1], env_inputfs, 'same');
env=tmp(scale/2:scale:size(tmp,1), :);

end