function demo_nlapgt(~, ~, ~)
%A demo of variable-Q ozgf filter properties

global XXX STACK META MODULES;

mdl = [];
mdl.Q_factor = 5;
mdl.N_order = 4;
mdl.center_freq_khz = 1;
mdl.time_align = false;
mdl.zeros = [0 0];
mdl.delayms = 0;

function sys = calcsys(mdl)

    CF = abs(mdl.center_freq_khz) * 1000;
    Q  = 1/sqrt(2) + 0.00001 + abs(mdl.Q_factor); %.00001 is for numerical stability
    N  = ceil(abs(mdl.N_order));
    
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
    
    if isempty(mdl.zeros)
        % Create the APGF to have a gain of 1 at DC
        sys = zpk([], poles, 1);
        sys = zpk([], poles, 1/dcgain(sys));
    else
        z = mdl.zeros * 1000;
        sys = zpk([], poles, 1);
        % Gain is 0 at DC for OZGF, so we just try to keep it "near" the
        % same levels by normalizing by the ringing frequency w_r. Why does
        % this get slightly erronious for very low Q factors?
        sys = zpk([], poles, (1/cos(theta)^N) / abs(evalfr(sys, w_r*1j)));         
    end

    % Correct for the time delay introduced by the wavelet
    if isfield(mdl, 'time_align') && mdl.time_align
        [imp, time] = impulse(sys);
        [~, idx] = max(abs(hilbert(imp)));
        t_align = time(idx);
    else
        t_align = 0;
    end
    
    sys.InputDelay = (abs(mdl.delayms) / 1000); % (milliseconds)
end
    
qq = 0.3:1:10;
len = length(qq);
samps = 1000;
w = logspace(3,4.5, samps);
matrix = zeros(samps, len);
ii = 1;
for q = qq
    mdl.Q_factor = q;
    sys = calcsys(mdl);
    [mag, phase] = bode(sys, w);
    matrix(:, ii) = log(squeeze(mag));
    ii = ii + 1;
end

figure;

%surf(matrix);
plot(matrix, 'Color', [0 0 0]);
xlabel('Log Frequency');
ylabel('Log Volume');
%zlabel('Log Neural Response');
title('Variable-Q APGT');    

% This code demos what happens if you try to apply log compression to one
% of the OZGF curves: It shifts it downward and makes the tails go to zero
% very quickly.
%
% figure();
% blah = matrix(:,end);
% ii = 1;
% for kk = linspace(-3, -1, len)   
%     matrix(:, ii) = nl_log([kk], blah);
%     ii = ii + 1;
% end
% plot(matrix, 'Color', [0 0 1]);
% xlabel('Log Volume');
% ylabel('Log Frequency');
% %zlabel('Log Neural Response');
% title('OZGF and then LOG COMPRESSION');    
% 
% 


end
